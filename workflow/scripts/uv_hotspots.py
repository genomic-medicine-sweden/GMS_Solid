"""
Count UV hotspot mutations at YYTTCCG ETS promoter sites in one DNA sample (experimental).

UV light forms pyrimidine dimers, and a cytosine in a dimer mutates C>T (CC>TT when
both are hit). At the YY dipyrimidine in front of ETS-bound TTCCG elements dimers form
far more often than elsewhere, so these positions act as sensitive UV mutation sensors.

The positions table (one row per position, hg38 1-based) has the columns
site, gene, rank, motif, motif_strand, cls, chr, pos, ref, alt, CpG, where cls is:
  Y1, Y2  cytosines in the YY dipyrimidine, the main hotspot
  C5      the cytosine after TT, a weaker hotspot, reported separately
  flank   C/G positions around the site, the in-sample C>T background

Bases are counted per fragment: overlapping mates count once when they agree and are
dropped when they disagree. Duplicate, secondary, supplementary, QC-failed and unmapped
reads, reads below min_mq and bases below min_bq are skipped. Counts are reported on
the pyrimidine (C) strand, so G>A on the reference is reported as C>T.

A Y1/Y2/C5 position is a hit with at least min_depth fragments, min_alt C>T fragments
and a C>T fraction of min_vaf; a hit at germline_vaf or above is flagged as a possible
germline variant. A flank position is left out of the background when its non-reference
fraction is at or above flank_max_vaf. Y1/Y2 are never CpG, so the hotspot rate is
compared with non-CpG flank positions.

Outputs: one row per position, one row per site, a one-row sample summary and a call
table (metric, value, cutoff, result) with provisional cut-offs, set from samples
without UV exposure only:

  UV_call                    not evaluable with fewer than min_sites_covered covered
                             sites; otherwise yes at uv_sites_hit or more sites with a
                             somatic (not possible germline) Y1/Y2 hit, possible at
                             possible_sites_hit or more or with an elevated
                             low_frequency_ratio, else no
  dnv_CCTT_fragments         one or more CC>TT fragments supports UV
  hotspot_alt_F1R2_fraction  outside orientation_bias..1-orientation_bias, with at least
                             orientation_min_alt C>T fragments, suggests FFPE artefacts
  low_frequency_ratio        C>T rate at Y1/Y2 positions without a hit against the non-CpG
                             flank rate; elevated at elevated_ratio or above with a one-sided
                             binomial p-value below elevated_p. It picks up hotspot
                             mutations below the hit threshold, e.g. at low tumour content
  manual_review              yes for a possible call (with the reason), possible germline
                             hits, or an orientation bias together with hits

All thresholds are in DEFAULTS; the ones in CONFIGURABLE are passed as rule params.
"""

import csv
import math
import statistics
from collections import defaultdict

import pysam

COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}
HOTSPOT_CLASSES = ("Y1", "Y2")
SCORED_CLASSES = ("Y1", "Y2", "C5")
POSITION_TABLE_COLUMNS = ["site", "gene", "rank", "motif", "motif_strand", "cls", "chr", "pos", "ref", "alt", "CpG"]
COUNT_COLUMNS = ["depth", "C", "C>T", "C>A", "C>G", "vaf_CT", "alt_F1R2", "alt_F2R1", "alt_unpaired", "discordant_pairs"]
POSITION_COLUMNS = POSITION_TABLE_COLUMNS + COUNT_COLUMNS + ["hit", "possible_germline", "flank_excluded"]
SITE_COLUMNS = [
    "site", "gene", "rank", "motif",
    "hotspot_min_depth", "hotspot_max_vaf", "hotspot_CT", "hotspot_CA", "hotspot_CG", "hotspot_depth_sum",
    "hit", "possible_germline", "dnv_depth", "dnv_CCTT",
    "C5_depth", "C5_CT", "C5_vaf", "C5_hit", "C5_possible_germline",
    "flank_CT", "flank_depth_sum", "flank_CT_rate", "flank_CpG_CT", "flank_CpG_depth_sum", "flank_CpG_CT_rate",
    "flank_excluded",
]
SUMMARY_COLUMNS = [
    "sample", "sites", "sites_covered", "median_hotspot_depth",
    "sites_hit", "sites_possible_germline", "dnv_CCTT_fragments",
    "hotspot_CT", "hotspot_depth", "hotspot_CT_rate", "hotspot_CA_rate", "hotspot_CG_rate", "hotspot_alt_F1R2_fraction",
    "flank_CT", "flank_depth", "flank_CT_rate", "flank_CpG_CT_rate", "flank_positions_excluded", "hotspot_vs_flank_CT",
    "C5_sites_hit", "C5_sites_possible_germline", "C5_CT", "C5_depth", "C5_CT_rate", "C5_vs_flank_CT",
]
CALL_COLUMNS = ["metric", "value", "cutoff", "result"]
DEFAULTS = {
    # Counting
    "min_bq": 30,
    "min_mq": 20,
    "flank_max_vaf": 0.1,
    # Hits
    "min_depth": 100,
    "min_alt": 3,
    "min_vaf": 0.02,
    "germline_vaf": 0.35,
    # Call
    "min_sites_covered": 80,
    "uv_sites_hit": 3,
    "possible_sites_hit": 1,
    "orientation_min_alt": 10,
    "orientation_bias": 0.1,
    "elevated_ratio": 1.5,
    "elevated_p": 0.01,
}
CONFIGURABLE = ("min_alt", "min_vaf", "uv_sites_hit", "possible_sites_hit", "elevated_ratio", "elevated_p")


def read_positions(path):
    """Positions grouped per site, in file order."""
    sites = defaultdict(list)
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            row["pos"] = int(row["pos"])
            row["CpG"] = row["CpG"] == "True"
            sites[row["site"]].append(row)
    return sites


def orientation(read):
    """F1R2 or F2R1 for paired reads, F or R for unpaired reads."""
    if not read.is_paired:
        return "R" if read.is_reverse else "F"
    return "F1R2" if read.is_read1 != read.is_reverse else "F2R1"


def collect_fragments(bam, chrom, positions, min_bq, min_mq):
    """
    Bases per fragment at 1-based positions: ({pos: {fragment: (base, orientation)}}, {pos: discordant pairs}).
    """
    wanted = {pos - 1 for pos in positions}
    fragments = defaultdict(dict)
    discordant = defaultdict(int)
    dropped = defaultdict(set)
    for read in bam.fetch(chrom, min(wanted), max(wanted) + 1):
        if (
            read.is_unmapped
            or read.is_duplicate
            or read.is_secondary
            or read.is_supplementary
            or read.is_qcfail
            or read.mapping_quality < min_mq
        ):
            continue
        qualities = read.query_qualities
        sequence = read.query_sequence
        read_orientation = orientation(read)
        for qpos, rpos in read.get_aligned_pairs(matches_only=True):
            if rpos not in wanted or qualities[qpos] < min_bq:
                continue
            pos = rpos + 1
            name = read.query_name
            if name in dropped[pos]:
                continue
            seen = fragments[pos].get(name)
            if seen is None:
                fragments[pos][name] = (sequence[qpos], read_orientation)
            elif seen[0] != sequence[qpos]:
                del fragments[pos][name]
                dropped[pos].add(name)
                discordant[pos] += 1
    return fragments, discordant


def pyrimidine_base(base, ref):
    """The base on the strand where the reference base is C."""
    return base if ref == "C" else COMPLEMENT.get(base, "N")


def count_position(row, fragments, discordant):
    counts = {"A": 0, "C": 0, "G": 0, "T": 0}
    alt_orientation = defaultdict(int)
    for base, read_orientation in fragments.values():
        b = pyrimidine_base(base, row["ref"])
        if b not in counts:
            continue
        counts[b] += 1
        if b == "T":
            alt_orientation[read_orientation] += 1
    depth = sum(counts.values())
    return {
        "depth": depth,
        "C": counts["C"],
        "C>T": counts["T"],
        "C>A": counts["A"],
        "C>G": counts["G"],
        "vaf_CT": counts["T"] / depth if depth else 0.0,
        "alt_F1R2": alt_orientation["F1R2"],
        "alt_F2R1": alt_orientation["F2R1"],
        "alt_unpaired": alt_orientation["F"] + alt_orientation["R"],
        "discordant_pairs": discordant.get(row["pos"], 0),
    }


def count_dnv(rows, fragments):
    """CC>TT at CC sites: (fragments covering both Y1 and Y2, fragments with T at both), or (None, None)."""
    hotspots = {row["cls"]: row for row in rows if row["cls"] in HOTSPOT_CLASSES}
    if len(hotspots) != 2:
        return None, None
    y1, y2 = hotspots["Y1"], hotspots["Y2"]
    covering = set(fragments.get(y1["pos"], {})) & set(fragments.get(y2["pos"], {}))
    dnv = sum(
        1
        for name in covering
        if pyrimidine_base(fragments[y1["pos"]][name][0], y1["ref"]) == "T"
        and pyrimidine_base(fragments[y2["pos"]][name][0], y2["ref"]) == "T"
    )
    return len(covering), dnv


def ratio(numerator, denominator):
    return numerator / denominator if denominator else float("nan")


def flag_position(row, params):
    """Add hit and possible_germline to a Y1/Y2/C5 row, or flank_excluded to a flank row."""
    if row["cls"] in SCORED_CLASSES:
        row["hit"] = (
            row["depth"] >= params["min_depth"]
            and row["C>T"] >= params["min_alt"]
            and row["vaf_CT"] >= params["min_vaf"]
        )
        row["possible_germline"] = row["hit"] and row["vaf_CT"] >= params["germline_vaf"]
    elif row["cls"] == "flank":
        non_reference = 1 - row["C"] / row["depth"] if row["depth"] else 0.0
        row["flank_excluded"] = non_reference >= params["flank_max_vaf"]


def summarise_site(rows, dnv):
    hotspots = [row for row in rows if row["cls"] in HOTSPOT_CLASSES]
    c5 = [row for row in rows if row["cls"] == "C5"]
    flank = [row for row in rows if row["cls"] == "flank" and not row["flank_excluded"]]
    flank_plain = [row for row in flank if not row["CpG"]]
    flank_cpg = [row for row in flank if row["CpG"]]
    flank_ct = sum(row["C>T"] for row in flank_plain)
    flank_depth = sum(row["depth"] for row in flank_plain)
    flank_cpg_ct = sum(row["C>T"] for row in flank_cpg)
    flank_cpg_depth = sum(row["depth"] for row in flank_cpg)
    return {
        **{column: rows[0][column] for column in ("site", "gene", "rank", "motif")},
        "hotspot_min_depth": min((row["depth"] for row in hotspots), default=0),
        "hotspot_max_vaf": max((row["vaf_CT"] for row in hotspots), default=0.0),
        "hotspot_CT": sum(row["C>T"] for row in hotspots),
        "hotspot_CA": sum(row["C>A"] for row in hotspots),
        "hotspot_CG": sum(row["C>G"] for row in hotspots),
        "hotspot_depth_sum": sum(row["depth"] for row in hotspots),
        "hit": any(row["hit"] for row in hotspots),
        "possible_germline": any(row["possible_germline"] for row in hotspots),
        "dnv_depth": dnv[0],
        "dnv_CCTT": dnv[1],
        "C5_depth": sum(row["depth"] for row in c5),
        "C5_CT": sum(row["C>T"] for row in c5),
        "C5_vaf": max((row["vaf_CT"] for row in c5), default=None),
        "C5_hit": any(row["hit"] for row in c5),
        "C5_possible_germline": any(row["possible_germline"] for row in c5),
        "flank_CT": flank_ct,
        "flank_depth_sum": flank_depth,
        "flank_CT_rate": ratio(flank_ct, flank_depth),
        "flank_CpG_CT": flank_cpg_ct,
        "flank_CpG_depth_sum": flank_cpg_depth,
        "flank_CpG_CT_rate": ratio(flank_cpg_ct, flank_cpg_depth),
        "flank_excluded": sum(1 for row in rows if row.get("flank_excluded")),
    }


def summarise_sample(sample, position_rows, site_rows, min_depth):
    hotspots = [row for row in position_rows if row["cls"] in HOTSPOT_CLASSES]
    hotspot_ct = sum(site["hotspot_CT"] for site in site_rows)
    hotspot_depth = sum(site["hotspot_depth_sum"] for site in site_rows)
    flank_ct = sum(site["flank_CT"] for site in site_rows)
    flank_depth = sum(site["flank_depth_sum"] for site in site_rows)
    c5_ct = sum(site["C5_CT"] for site in site_rows)
    c5_depth = sum(site["C5_depth"] for site in site_rows)
    f1r2 = sum(row["alt_F1R2"] for row in hotspots)
    f2r1 = sum(row["alt_F2R1"] for row in hotspots)
    hotspot_rate = ratio(hotspot_ct, hotspot_depth)
    flank_rate = ratio(flank_ct, flank_depth)
    c5_rate = ratio(c5_ct, c5_depth)
    return {
        "sample": sample,
        "sites": len(site_rows),
        "sites_covered": sum(1 for site in site_rows if site["hotspot_min_depth"] >= min_depth),
        "median_hotspot_depth": statistics.median(row["depth"] for row in hotspots) if hotspots else 0,
        "sites_hit": sum(site["hit"] for site in site_rows),
        "sites_possible_germline": sum(site["possible_germline"] for site in site_rows),
        "dnv_CCTT_fragments": sum(site["dnv_CCTT"] or 0 for site in site_rows),
        "hotspot_CT": hotspot_ct,
        "hotspot_depth": hotspot_depth,
        "hotspot_CT_rate": hotspot_rate,
        "hotspot_CA_rate": ratio(sum(site["hotspot_CA"] for site in site_rows), hotspot_depth),
        "hotspot_CG_rate": ratio(sum(site["hotspot_CG"] for site in site_rows), hotspot_depth),
        "hotspot_alt_F1R2_fraction": ratio(f1r2, f1r2 + f2r1),
        "flank_CT": flank_ct,
        "flank_depth": flank_depth,
        "flank_CT_rate": flank_rate,
        "flank_CpG_CT_rate": ratio(
            sum(site["flank_CpG_CT"] for site in site_rows), sum(site["flank_CpG_depth_sum"] for site in site_rows)
        ),
        "flank_positions_excluded": sum(site["flank_excluded"] for site in site_rows),
        "hotspot_vs_flank_CT": ratio(hotspot_rate, flank_rate),
        "C5_sites_hit": sum(site["C5_hit"] for site in site_rows),
        "C5_sites_possible_germline": sum(site["C5_possible_germline"] for site in site_rows),
        "C5_CT": c5_ct,
        "C5_depth": c5_depth,
        "C5_CT_rate": c5_rate,
        "C5_vs_flank_CT": ratio(c5_rate, flank_rate),
    }


def uv_hotspots(bam_path, positions_path, sample, params):
    """Count all positions in one BAM: (position rows, site rows, sample summary)."""
    position_rows = []
    site_rows = []
    with pysam.AlignmentFile(bam_path) as bam:
        for rows in read_positions(positions_path).values():
            fragments, discordant = collect_fragments(
                bam, rows[0]["chr"], [row["pos"] for row in rows], params["min_bq"], params["min_mq"]
            )
            site_positions = []
            for row in rows:
                counted = {**row, **count_position(row, fragments.get(row["pos"], {}), discordant)}
                flag_position(counted, params)
                site_positions.append(counted)
            position_rows.extend(site_positions)
            site_rows.append(summarise_site(site_positions, count_dnv(rows, fragments)))
    return position_rows, site_rows, summarise_sample(sample, position_rows, site_rows, params["min_depth"])


def binomial_upper_tail(k, n, p):
    """P(X >= k) for X ~ Binomial(n, p)."""
    if k <= 0:
        return 1.0
    if k > n or p <= 0:
        return 0.0
    if p >= 1:
        return 1.0
    log_p, log_q = math.log(p), math.log1p(-p)

    def pmf(i):
        return math.exp(math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) + i * log_p + (n - i) * log_q)

    if k <= n * p:
        # Tail holds at least about half the mass, so the complement is accurate
        return max(0.0, 1.0 - sum(pmf(i) for i in range(k)))
    # Above the mean the terms only decrease, so sum until they no longer matter
    total = 0.0
    for i in range(k, n + 1):
        term = pmf(i)
        total += term
        if term < total * 1e-16:
            break
    return min(1.0, total)


def low_frequency(position_rows, flank_rate, cutoffs):
    """
    Hotspot C>T below the hit threshold against the flank: (ratio, p-value, result). Hit positions are left out,
    so clonal hits do not drive the ratio.
    """
    rows = [row for row in position_rows if row["cls"] in HOTSPOT_CLASSES and not row["hit"]]
    ct = sum(row["C>T"] for row in rows)
    depth = sum(row["depth"] for row in rows)
    if not depth or flank_rate != flank_rate or flank_rate <= 0:
        return float("nan"), float("nan"), "not evaluable"
    ratio_value = (ct / depth) / flank_rate
    p_value = binomial_upper_tail(ct, depth, flank_rate)
    elevated = ratio_value >= cutoffs["elevated_ratio"] and p_value < cutoffs["elevated_p"]
    return ratio_value, p_value, "elevated" if elevated else "not elevated"


def uv_call(summary, position_rows, cutoffs):
    """Call table rows (metric, value, cutoff, result) with an overall UV call and a manual review flag."""
    somatic_hit_sites = {row["site"] for row in position_rows if row["cls"] in HOTSPOT_CLASSES and row["hit"]}
    germline_sites = {row["site"] for row in position_rows if row["cls"] in HOTSPOT_CLASSES and row["possible_germline"]}
    somatic_hit_sites -= germline_sites
    hits = len(somatic_hit_sites)
    covered = summary["sites_covered"]
    evaluable = covered >= cutoffs["min_sites_covered"]
    ratio_value, p_value, low_frequency_result = low_frequency(position_rows, summary["flank_CT_rate"], cutoffs)

    if not evaluable:
        call = "not evaluable"
    elif hits >= cutoffs["uv_sites_hit"]:
        call = "yes"
    elif hits >= cutoffs["possible_sites_hit"] or low_frequency_result == "elevated":
        call = "possible"
    else:
        call = "no"

    oriented = summary["hotspot_CT"] - sum(row["alt_unpaired"] for row in position_rows if row["cls"] in HOTSPOT_CLASSES)
    f1r2 = summary["hotspot_alt_F1R2_fraction"]
    bias = cutoffs["orientation_bias"]
    if oriented < cutoffs["orientation_min_alt"]:
        orientation = "too few fragments"
    elif bias <= f1r2 <= 1 - bias:
        orientation = "ok"
    else:
        orientation = "possible FFPE artefact"

    reasons = []
    if call == "possible" and hits >= cutoffs["possible_sites_hit"]:
        reasons.append("possible UV call")
    if call == "possible" and low_frequency_result == "elevated":
        reasons.append("elevated hotspot/flank ratio")
    if germline_sites:
        reasons.append("possible germline hit")
    if hits and orientation == "possible FFPE artefact":
        reasons.append("orientation bias")

    hit_cutoff = f">={cutoffs['uv_sites_hit']} yes, >={cutoffs['possible_sites_hit']} possible"
    ratio_cutoff = f">={cutoffs['elevated_ratio']:g} and p<{cutoffs['elevated_p']:g}"
    return [
        {
            "metric": "UV_call",
            "value": call,
            "cutoff": f"provisional, sites_hit {hit_cutoff}; elevated low_frequency_ratio gives possible",
            "result": call,
        },
        {
            "metric": "sites_covered",
            "value": covered,
            "cutoff": f">={cutoffs['min_sites_covered']}",
            "result": "ok" if evaluable else "too few",
        },
        {
            "metric": "sites_hit",
            "value": hits,
            "cutoff": f"{hit_cutoff} (possible germline excluded)",
            "result": call if evaluable else "not evaluable",
        },
        {
            "metric": "sites_possible_germline",
            "value": len(germline_sites),
            "cutoff": "not counted as hits",
            "result": "review" if germline_sites else "none",
        },
        {
            "metric": "dnv_CCTT_fragments",
            "value": summary["dnv_CCTT_fragments"],
            "cutoff": ">=1 supports UV",
            "result": "supports UV" if summary["dnv_CCTT_fragments"] else "none",
        },
        {
            "metric": "hotspot_alt_F1R2_fraction",
            "value": f1r2,
            "cutoff": f"{bias:g}-{1 - bias:g} with >={cutoffs['orientation_min_alt']} C>T fragments",
            "result": orientation,
        },
        {
            "metric": "low_frequency_ratio",
            "value": ratio_value,
            "cutoff": f"{ratio_cutoff}; hotspot C>T at non-hit positions vs non-CpG flank",
            "result": low_frequency_result,
        },
        {"metric": "low_frequency_p", "value": p_value, "cutoff": f"<{cutoffs['elevated_p']:g}", "result": low_frequency_result},
        {"metric": "manual_review", "value": "yes" if reasons else "no", "cutoff": "", "result": ", ".join(reasons)},
    ]


def format_value(value):
    if value is None:
        return ""
    if isinstance(value, float):
        return "NA" if value != value else f"{value:.6g}"
    return value


def write_tsv(path, rows, columns):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t", lineterminator="\n", restval="", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: format_value(value) for key, value in row.items()})


def main(snakemake):
    params = {**DEFAULTS, **{key: getattr(snakemake.params, key) for key in CONFIGURABLE}}
    position_rows, site_rows, summary = uv_hotspots(
        snakemake.input.bam, snakemake.input.positions, snakemake.params.sample, params
    )
    write_tsv(snakemake.output.positions, position_rows, POSITION_COLUMNS)
    write_tsv(snakemake.output.sites, site_rows, SITE_COLUMNS)
    write_tsv(snakemake.output.summary, [summary], SUMMARY_COLUMNS)
    write_tsv(snakemake.output.call, uv_call(summary, position_rows, params), CALL_COLUMNS)


if __name__ == "__main__":
    main(snakemake)

import csv
import math
import os
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest import mock

import pysam

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import uv_hotspots as uv  # noqa: E402
from uv_hotspots import binomial_upper_tail, main, pyrimidine_base, uv_call, uv_hotspots  # noqa: E402

# The synthetic sites have 50-100 fragments, so the depth cut-off is lowered
PARAMS = {**uv.DEFAULTS, "min_depth": 10}

# 120 bp reference: plus-strand CCTTCCG at 41-47, minus-strand CCTTCCG (CGGAAGG on the reference) at 81-87
REFERENCE = list("A" * 120)
for start, motif in ((41, "CCTTCCG"), (81, "CGGAAGG")):
    REFERENCE[start - 1:start - 1 + len(motif)] = motif
REFERENCE[30 - 1] = "C"  # non-CpG flank
REFERENCE[35 - 1:35 + 1] = ["C", "G"]  # CpG flank at 35
REFERENCE[95 - 1] = "G"  # non-CpG flank, minus-strand site
REFERENCE = "".join(REFERENCE)

POSITIONS = [
    # site, motif_strand, cls, pos, ref, CpG
    ("SITE_A", "+", "Y1", 41, "C", False),
    ("SITE_A", "+", "Y2", 42, "C", False),
    ("SITE_A", "+", "C5", 45, "C", False),
    ("SITE_A", "+", "flank", 30, "C", False),
    ("SITE_A", "+", "flank", 35, "C", True),
    ("SITE_B", "-", "Y1", 87, "G", False),
    ("SITE_B", "-", "Y2", 86, "G", False),
    ("SITE_B", "-", "C5", 83, "G", False),
    ("SITE_B", "-", "flank", 95, "G", False),
]


def mutate(changes, start, end):
    """Reference bases start..end (1-based, inclusive) with changes applied."""
    sequence = list(REFERENCE)
    for pos, base in changes:
        sequence[pos - 1] = base
    return "".join(sequence[start - 1:end])


def make_read(name, sequence, start, read1, reverse, mapq=60, duplicate=False):
    read = pysam.AlignedSegment()
    read.query_name = name
    read.query_sequence = sequence
    read.reference_id = 0
    read.reference_start = start - 1
    read.cigarstring = f"{len(sequence)}M"
    read.mapping_quality = mapq
    read.query_qualities = pysam.qualitystring_to_array("I" * len(sequence))
    read.flag = 1 | 2 | (64 if read1 else 128) | (16 if reverse else 32) | (1024 if duplicate else 0)
    read.next_reference_id = 0
    read.next_reference_start = start - 1
    return read


def pair(name, changes, region, mate2_changes=None, **kwargs):
    """Fully overlapping mates covering region (start, end)."""
    start, end = region
    mate2 = changes if mate2_changes is None else mate2_changes
    return [
        make_read(name, mutate(changes, start, end), start, True, False, **kwargs),
        make_read(name, mutate(mate2, start, end), start, False, True, **kwargs),
    ]


SITE_A_REGION = (1, 60)
SITE_B_REGION = (61, 120)


def build_reads():
    reads = []
    # SITE_A: 100 fragments
    for i in range(100):
        changes = []
        if i < 10:
            changes.append((41, "T"))  # C>T at Y1
        elif i < 15:
            changes += [(41, "T"), (42, "T")]  # CC>TT
        if 20 <= i < 26:
            changes.append((45, "T"))  # C>T at C5
        if i < 50:
            changes.append((30, "T"))  # heterozygous germline at a flank position
        if 60 <= i < 62:
            changes.append((35, "T"))  # C>T at the CpG flank
        if i in (90, 91):
            reads += pair(f"a{i}", changes + [(41, "T")], SITE_A_REGION, mate2_changes=changes)  # mates disagree at Y1
        else:
            reads += pair(f"a{i}", changes, SITE_A_REGION)
    # Skipped reads: duplicate and low mapping quality, both with C>T at Y2
    reads += pair("dup", [(42, "T")], SITE_A_REGION, duplicate=True)
    reads += pair("lowmq", [(42, "T")], SITE_A_REGION, mapq=10)
    # SITE_B: 50 fragments, G>A on the reference is C>T on the motif strand
    for i in range(50):
        changes = [(87, "A")] if i < 20 else []
        if i == 0:
            changes.append((95, "A"))
        reads += pair(f"b{i}", changes, SITE_B_REGION)
    return sorted(reads, key=lambda read: read.reference_start)


class TestUvHotspots(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp_dir = tempfile.TemporaryDirectory()
        cls.bam = os.path.join(cls.temp_dir.name, "S_T.bam")
        header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chrT", "LN": len(REFERENCE)}]}
        with pysam.AlignmentFile(cls.bam, "wb", header=header) as f:
            for read in build_reads():
                f.write(read)
        pysam.index(cls.bam)
        cls.positions = os.path.join(cls.temp_dir.name, "positions.tsv")
        with open(cls.positions, "w") as f:
            f.write("site\tgene\trank\tmotif\tmotif_strand\tcls\tchr\tpos\tref\talt\tCpG\n")
            for site, strand, cls_, pos, ref, cpg in POSITIONS:
                alt = "T" if ref == "C" else "A"
                f.write(f"{site}\tGENE_{site}\t1\tCCTTCCG\t{strand}\t{cls_}\tchrT\t{pos}\t{ref}\t{alt}\t{cpg}\n")
        cls.position_rows, cls.site_rows, cls.summary = uv_hotspots(cls.bam, cls.positions, "S_T", PARAMS)
        cls.by_pos = {row["pos"]: row for row in cls.position_rows}
        cls.by_site = {row["site"]: row for row in cls.site_rows}

    @classmethod
    def tearDownClass(cls):
        cls.temp_dir.cleanup()

    def test_pyrimidine_base(self):
        self.assertEqual(pyrimidine_base("T", "C"), "T")
        self.assertEqual(pyrimidine_base("A", "G"), "T")
        self.assertEqual(pyrimidine_base("N", "G"), "N")

    def test_discordant_mates_dropped(self):
        y1 = self.by_pos[41]
        self.assertEqual(y1["depth"], 98)
        self.assertEqual(y1["C>T"], 15)
        self.assertEqual(y1["discordant_pairs"], 2)

    def test_duplicate_and_low_mapq_skipped(self):
        y2 = self.by_pos[42]
        self.assertEqual(y2["depth"], 100)
        self.assertEqual(y2["C>T"], 5)

    def test_minus_strand_reported_as_c_to_t(self):
        y1 = self.by_pos[87]
        self.assertEqual((y1["depth"], y1["C>T"], y1["C"]), (50, 20, 30))
        self.assertAlmostEqual(y1["vaf_CT"], 0.4)
        self.assertTrue(y1["possible_germline"])

    def test_orientation_split(self):
        y1 = self.by_pos[41]
        self.assertEqual((y1["alt_F1R2"], y1["alt_F2R1"]), (15, 0))

    def test_dnv(self):
        site = self.by_site["SITE_A"]
        self.assertEqual((site["dnv_depth"], site["dnv_CCTT"]), (98, 5))

    def test_hits(self):
        site_a, site_b = self.by_site["SITE_A"], self.by_site["SITE_B"]
        self.assertTrue(site_a["hit"])
        self.assertFalse(site_a["possible_germline"])
        self.assertTrue(site_b["hit"])
        self.assertTrue(site_b["possible_germline"])

    def test_c5_reported_separately(self):
        site = self.by_site["SITE_A"]
        self.assertEqual((site["C5_depth"], site["C5_CT"]), (100, 6))
        self.assertTrue(site["C5_hit"])
        self.assertEqual(site["hotspot_CT"], 20)  # Y1 + Y2 only

    def test_flank_exclusion_and_cpg_split(self):
        self.assertTrue(self.by_pos[30]["flank_excluded"])
        self.assertFalse(self.by_pos[35]["flank_excluded"])
        site_a, site_b = self.by_site["SITE_A"], self.by_site["SITE_B"]
        self.assertEqual(site_a["flank_excluded"], 1)
        self.assertEqual(site_a["flank_depth_sum"], 0)
        self.assertTrue(math.isnan(site_a["flank_CT_rate"]))
        self.assertEqual((site_a["flank_CpG_CT"], site_a["flank_CpG_depth_sum"]), (2, 100))
        self.assertEqual((site_b["flank_CT"], site_b["flank_depth_sum"]), (1, 50))

    def test_summary(self):
        s = self.summary
        self.assertEqual((s["sites"], s["sites_covered"], s["sites_hit"], s["sites_possible_germline"]), (2, 2, 2, 1))
        self.assertEqual(s["dnv_CCTT_fragments"], 5)
        self.assertEqual((s["hotspot_CT"], s["hotspot_depth"]), (40, 298))
        self.assertEqual((s["flank_CT"], s["flank_depth"]), (1, 50))
        self.assertAlmostEqual(s["hotspot_vs_flank_CT"], (40 / 298) / (1 / 50))
        self.assertEqual((s["C5_sites_hit"], s["C5_CT"], s["C5_depth"]), (1, 6, 150))
        self.assertEqual(s["flank_positions_excluded"], 1)

    def test_main_writes_outputs(self):
        out = {
            name: os.path.join(self.temp_dir.name, f"S_T.uv_hotspots.{name}.tsv")
            for name in ("positions", "sites", "summary", "call")
        }
        snakemake = SimpleNamespace(
            input=SimpleNamespace(bam=self.bam, positions=self.positions),
            params=SimpleNamespace(sample="S_T", **{key: uv.DEFAULTS[key] for key in uv.CONFIGURABLE}),
            output=SimpleNamespace(**out),
        )
        with mock.patch.dict(uv.DEFAULTS, {"min_depth": 10, "min_sites_covered": 2}):
            main(snakemake)
        with open(out["call"]) as f:
            call = {row["metric"]: row for row in csv.DictReader(f, delimiter="\t")}
        # SITE_B's hit is possible germline, so one somatic hit remains
        self.assertEqual(call["UV_call"]["value"], "possible")
        self.assertEqual(call["sites_hit"]["value"], "1")
        self.assertEqual(call["manual_review"]["value"], "yes")
        with open(out["summary"]) as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["sample"], "S_T")
        self.assertEqual(rows[0]["sites_hit"], "2")
        with open(out["positions"]) as f:
            self.assertEqual(sum(1 for _ in f), len(POSITIONS) + 1)
        with open(out["sites"]) as f:
            sites = list(csv.DictReader(f, delimiter="\t"))
        self.assertEqual(sites[0]["flank_CT_rate"], "NA")
        self.assertEqual(sites[1]["dnv_depth"], "50")


def hotspot_row(site, hit=False, germline=False, unpaired=0, ct=0, depth=1000):
    return {
        "cls": "Y1", "site": site, "hit": hit, "possible_germline": germline, "alt_unpaired": unpaired,
        "C>T": ct, "depth": depth,
    }


# Background: 40 C>T in 200000 fragments, the flank rate of 2e-4
BACKGROUND = hotspot_row("BG", ct=40, depth=200000)


def call_summary(**values):
    summary = {
        "sites_covered": 100,
        "hotspot_CT": 40,
        "hotspot_depth": 200000,
        "hotspot_alt_F1R2_fraction": 0.5,
        "flank_CT_rate": 2e-4,
        "hotspot_vs_flank_CT": 1.0,
        "dnv_CCTT_fragments": 0,
    }
    summary.update(values)
    return summary


class TestUvCall(unittest.TestCase):
    def call(self, summary, rows):
        return {row["metric"]: row for row in uv_call(summary, rows, uv.DEFAULTS)}

    def test_binomial_upper_tail(self):
        from scipy.stats import binom

        for k, n, p in ((0, 10, 0.3), (1, 2, 0.5), (2, 2, 0.5), (40, 200000, 2e-4), (80, 200000, 2e-4), (500, 200000, 2e-4)):
            expected = binom.sf(k - 1, n, p)
            self.assertAlmostEqual(binomial_upper_tail(k, n, p), expected, delta=1e-12 + 1e-9 * expected)
        self.assertEqual(binomial_upper_tail(5, 4, 0.5), 0.0)

    def test_yes_at_three_somatic_hits(self):
        rows = [hotspot_row(f"S{i}", hit=True) for i in range(3)] + [BACKGROUND]
        call = self.call(call_summary(), rows)
        self.assertEqual(call["UV_call"]["value"], "yes")
        self.assertEqual(call["manual_review"]["value"], "no")

    def test_germline_hits_not_counted(self):
        rows = [hotspot_row(f"S{i}", hit=True) for i in range(2)] + [hotspot_row("S2", hit=True, germline=True), BACKGROUND]
        call = self.call(call_summary(), rows)
        self.assertEqual((call["UV_call"]["value"], call["sites_hit"]["value"]), ("possible", 2))
        self.assertEqual(call["sites_possible_germline"]["value"], 1)
        self.assertEqual(call["manual_review"]["result"], "possible UV call, possible germline hit")

    def test_no_and_not_evaluable(self):
        self.assertEqual(self.call(call_summary(), [BACKGROUND])["UV_call"]["value"], "no")
        call = self.call(call_summary(sites_covered=79), [hotspot_row(f"S{i}", hit=True) for i in range(5)])
        self.assertEqual(call["UV_call"]["value"], "not evaluable")
        self.assertEqual(call["sites_covered"]["result"], "too few")

    def test_orientation(self):
        rows = [hotspot_row("S0", hit=True), BACKGROUND]
        call = self.call(call_summary(hotspot_alt_F1R2_fraction=0.95), rows)
        self.assertEqual(call["hotspot_alt_F1R2_fraction"]["result"], "possible FFPE artefact")
        self.assertIn("orientation bias", call["manual_review"]["result"])
        call = self.call(call_summary(hotspot_CT=5, hotspot_alt_F1R2_fraction=1.0), rows)
        self.assertEqual(call["hotspot_alt_F1R2_fraction"]["result"], "too few fragments")

    def test_low_frequency_raises_no_to_possible(self):
        # 160 C>T in 200000 is 4x the flank rate of 2e-4
        call = self.call(call_summary(), [hotspot_row("S0", ct=160, depth=200000)])
        self.assertEqual(call["low_frequency_ratio"]["result"], "elevated")
        self.assertAlmostEqual(call["low_frequency_ratio"]["value"], 4.0)
        self.assertLess(call["low_frequency_p"]["value"], 1e-20)
        self.assertEqual(call["UV_call"]["value"], "possible")
        self.assertEqual(call["manual_review"]["result"], "elevated hotspot/flank ratio")

    def test_low_frequency_does_not_give_yes(self):
        rows = [hotspot_row("S0", hit=True, ct=300), hotspot_row("S1", ct=160, depth=200000)]
        call = self.call(call_summary(), rows)
        self.assertEqual(call["UV_call"]["value"], "possible")
        self.assertEqual(call["manual_review"]["result"], "possible UV call, elevated hotspot/flank ratio")

    def test_low_frequency_leaves_out_hit_positions(self):
        rows = [hotspot_row(f"S{i}", hit=True, ct=300) for i in range(3)] + [BACKGROUND]
        call = self.call(call_summary(), rows)
        self.assertAlmostEqual(call["low_frequency_ratio"]["value"], 1.0)
        self.assertEqual(call["low_frequency_ratio"]["result"], "not elevated")
        self.assertEqual(call["UV_call"]["value"], "yes")

    def test_low_frequency_not_elevated_or_not_evaluable(self):
        call = self.call(call_summary(), [hotspot_row("S0", ct=48, depth=200000)])
        self.assertEqual(call["low_frequency_ratio"]["result"], "not elevated")
        self.assertEqual(call["UV_call"]["value"], "no")
        call = self.call(call_summary(flank_CT_rate=float("nan")), [BACKGROUND])
        self.assertEqual(call["low_frequency_ratio"]["result"], "not evaluable")

    def test_dnv_supports_uv(self):
        call = self.call(call_summary(dnv_CCTT_fragments=2), [BACKGROUND])
        self.assertEqual(call["dnv_CCTT_fragments"]["result"], "supports UV")


if __name__ == "__main__":
    unittest.main()

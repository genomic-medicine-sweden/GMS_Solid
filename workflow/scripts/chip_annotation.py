__author__ = "Jonas A"
__copyright__ = "Copyright 2025, Jonas A"
__email__ = "jonas.almlof@igp.uu.se"
__license__ = "GPL-3"

import re
import yaml
import numpy as np
import pysam

LOF_CONSEQUENCES = frozenset({
    "frameshift_variant", "stop_gained", "stop_lost",
    "splice_donor_variant", "splice_acceptor_variant", "start_lost",
})

NEW_INFO = [
    ("CHIP_GENE", "1", "Integer",
     "CHIP gene tier: 0=not CHIP gene, 1=minor CHIP gene, 2=major CHIP gene (DNMT3A/TET2/ASXL1)"),
    ("CHIP_VAF", "1", "Integer",
     "CHIP VAF tier: 0=outside range; 1=vaf_tier1_min to vaf_tier1_max (default 0-10%); "
     "2=vaf_tier2_min to vaf_tier2_max (default 1-3%)"),
    ("CHIP_FRAG", "1", "Integer",
     "CHIP fragment length flag: 1=alt reads both above absolute length threshold (frag_abs_threshold) "
     "and above length ratio relative to ref reads (frag_ratio_threshold); "
     "-1=alt reads below short fragment threshold (frag_short_threshold), suggesting tumour-derived origin; "
     "0=not flagged; missing if fewer than min_alt_reads supporting reads or no ref reads"),
    ("CHIP_FRAG_ALT", "1", "Float",
     "Mean fragment length of alt-supporting reads"),
    ("CHIP_FRAG_REF", "1", "Float",
     "Mean fragment length of ref-supporting reads"),
    ("CHIP_HOTSPOT", "1", "Integer",
     "1 if variant matches a known CHIP hotspot position"),
    ("CHIP_COSMIC_HEMATO", "1", "Integer",
     "Count of haematological tumour samples with this variant in COSMIC (0 if not found)"),
    ("CHIP_CONSEQUENCE", "1", "Integer",
     "1 if LoF consequence in a LoF-CHIP gene or missense in a missense-CHIP gene"),
    ("CHIP_SAME_CLONE", "1", "Integer",
     "1 if this variant is in a clonal_context_gene (e.g. TP53) and co-occurs in this sample "
     "with a CHIP-gene variant at a similar VAF (within vaf_clonal_tolerance), "
     "indicating evidence of a shared haematopoietic clone. Contributes 3 to CHIP_SCORE."),
    ("CHIP_SCORE", "1", "Integer",
     "CHIP combined score: CHIP_GENE (0-2) + CHIP_VAF (0-2) + CHIP_FRAG (-1/0/1) + "
     "min(CHIP_HOTSPOT*4, 4) + min(CHIP_COSMIC_HEMATO, 4) + CHIP_CONSEQUENCE (0-1) + "
     "CHIP_SAME_CLONE*3. Filter threshold: >= 5"),
]


def get_vep_field_index(header, field_name):
    for record in header.records:
        if record.key == "INFO" and record.get("ID") == "CSQ":
            match = re.search(r"Format: ([^\"]+)", record.get("Description", ""))
            if match:
                fields = match.group(1).strip().split("|")
                if field_name in fields:
                    return fields.index(field_name)
    return None


def lookup_tabix_vcf(tabix, chrom, pos, ref, alt):
    """Return the INFO string if exact variant (chrom/pos/ref/alt) is found, else None."""
    try:
        for row in tabix.fetch(chrom, pos - 1, pos):
            cols = row.split("\t")
            if int(cols[1]) == pos and cols[3] == ref and cols[4] == alt:
                return cols[7] if len(cols) > 7 else "."
    except (ValueError, StopIteration):
        pass
    return None


def pileup_at_variant(bam, chrom, pos, ref, alt):
    """
    Single-pass BAM pileup at a variant position.
    Returns (alt_lengths, ref_lengths, alt_fwd, alt_rev).
    For SNVs, reads are assigned to alt/ref by base call.
    For indels, all reads are collected as alt (strand counts included).
    """
    alt_lengths, ref_lengths = [], []
    alt_fwd = alt_rev = 0
    is_snv = len(ref) == 1 and len(alt) == 1

    for col in bam.pileup(chrom, pos - 1, pos, truncate=True,
                          min_base_quality=0, stepper="nofilter"):
        if col.reference_pos != pos - 1:
            continue
        for pread in col.pileups:
            if pread.is_del or pread.is_refskip:
                continue
            read = pread.alignment
            if not read.is_proper_pair or read.is_unmapped:
                continue
            tlen = abs(read.template_length)
            if tlen == 0 or tlen > 1000:
                continue
            if is_snv:
                base = read.query_sequence[pread.query_position]
                if base == alt:
                    alt_lengths.append(tlen)
                    alt_rev += read.is_reverse
                    alt_fwd += not read.is_reverse
                elif base == ref:
                    ref_lengths.append(tlen)
            else:
                alt_lengths.append(tlen)
                alt_rev += read.is_reverse
                alt_fwd += not read.is_reverse

    return alt_lengths, ref_lengths, alt_fwd, alt_rev


# --- per-flag functions ---

def chip_gene_flag(rec, symbol_idx, major_genes, minor_genes):
    if symbol_idx is None or "CSQ" not in rec.info:
        return 0
    best = 0
    for csq in rec.info["CSQ"]:
        parts = csq.split("|")
        if len(parts) <= symbol_idx:
            continue
        sym = parts[symbol_idx]
        if sym in major_genes:
            best = 2
        elif sym in minor_genes and best < 2:
            best = 1
    return best


def chip_vaf_flag(rec, vaf_tier1_min, vaf_tier1_max, vaf_tier2_min, vaf_tier2_max):
    if "AF" not in rec.info:
        return 0
    af = rec.info["AF"]
    if isinstance(af, (tuple, list)):
        af = af[0] if af else None
    if af is None:
        return 0
    if vaf_tier2_min <= af <= vaf_tier2_max:
        return 2
    if vaf_tier1_min <= af <= vaf_tier1_max:
        return 1
    return 0


def chip_frag_and_sbs_flags(rec, bam, min_alt_reads, frag_abs_threshold, frag_ratio_threshold, frag_short_threshold):
    alts = [a for a in rec.alts if a != "<NON_REF>"] if rec.alts else []
    if not alts:
        return None, None, None

    alt_lens, ref_lens, alt_fwd, alt_rev = pileup_at_variant(
        bam, rec.contig, rec.pos, rec.ref, alts[0]
    )

    if len(alt_lens) < min_alt_reads:
        return None, None, None

    alt_mean = float(np.mean(alt_lens))
    ref_mean = float(np.mean(ref_lens)) if ref_lens else None

    if alt_mean < frag_short_threshold:
        chip_frag = -1
    elif ref_mean is not None and (alt_mean / ref_mean) > frag_ratio_threshold and alt_mean > frag_abs_threshold:
        chip_frag = 1
    else:
        chip_frag = 0
    frag_alt_out = round(alt_mean, 1)
    frag_ref_out = round(ref_mean, 1) if ref_mean is not None else None

    return chip_frag, frag_alt_out, frag_ref_out


def chip_hotspot_flag(rec, hotspot_tabix):
    alts = [a for a in rec.alts if a != "<NON_REF>"] if rec.alts else []
    if not alts or hotspot_tabix is None:
        return 0
    return 1 if lookup_tabix_vcf(hotspot_tabix, rec.contig, rec.pos, rec.ref, alts[0]) is not None else 0


def chip_cosmic_hemato_flag(rec, cosmic_tabix, symbol_idx, cosmic_exclude_genes, cosmic_min_count):
    alts = [a for a in rec.alts if a != "<NON_REF>"] if rec.alts else []
    if not alts or cosmic_tabix is None:
        return 0
    if cosmic_exclude_genes and symbol_idx is not None and "CSQ" in rec.info:
        for csq in rec.info["CSQ"]:
            parts = csq.split("|")
            if len(parts) > symbol_idx and parts[symbol_idx] in cosmic_exclude_genes:
                return 0
    info_str = lookup_tabix_vcf(cosmic_tabix, rec.contig, rec.pos, rec.ref, alts[0])
    if info_str is None:
        return 0
    for field in info_str.split(";"):
        if field.startswith("CNT="):
            try:
                cnt = int(field[4:])
                return cnt if cnt >= cosmic_min_count else 0
            except ValueError:
                pass
    return 0


def chip_consequence_flag(rec, symbol_idx, consequence_idx, lof_genes, missense_genes):
    if None in (symbol_idx, consequence_idx) or "CSQ" not in rec.info:
        return 0
    max_idx = max(symbol_idx, consequence_idx)
    for csq in rec.info["CSQ"]:
        parts = csq.split("|")
        if len(parts) <= max_idx:
            continue
        sym = parts[symbol_idx]
        consequences = set(parts[consequence_idx].split("&"))
        if sym in lof_genes and consequences & LOF_CONSEQUENCES:
            return 1
        if sym in missense_genes and "missense_variant" in consequences:
            return 1
    return 0


def get_af(rec):
    """Return raw AF from INFO field, or None if absent."""
    if "AF" not in rec.info:
        return None
    af = rec.info["AF"]
    if isinstance(af, (tuple, list)):
        af = af[0] if af else None
    return af


def is_in_gene_set(rec, symbol_idx, gene_set):
    """Return True if any CSQ transcript for this record has SYMBOL in gene_set."""
    if not gene_set or symbol_idx is None or "CSQ" not in rec.info:
        return False
    for csq in rec.info["CSQ"]:
        parts = csq.split("|")
        if len(parts) > symbol_idx and parts[symbol_idx] in gene_set:
            return True
    return False


def chip_same_clone_flag(vaf, chip_partners, tolerance):
    """Return 1 if a CHIP-gene variant exists at a similar VAF (shared clone evidence)."""
    if vaf is None or vaf <= 0:
        return 0
    for partner_vaf, _ in chip_partners:
        if partner_vaf is not None and partner_vaf > 0:
            if abs(vaf - partner_vaf) / max(vaf, partner_vaf) <= tolerance:
                return 1
    return 0


def main(snakemake):
    with open(snakemake.input.chip_genes_file) as fh:
        chip_genes = yaml.safe_load(fh)
    major_genes = set(chip_genes.get("major_chip_genes", []))
    minor_genes = set(chip_genes.get("minor_chip_genes", []))
    lof_genes = set(chip_genes.get("lof_chip_genes", []))
    missense_genes = set(chip_genes.get("missense_chip_genes", []))
    cosmic_exclude_genes = set(chip_genes.get("cosmic_exclude_genes", []))
    clonal_context_genes = set(chip_genes.get("clonal_context_genes", []))
    min_alt_reads = snakemake.params.min_alt_reads
    frag_abs_threshold = snakemake.params.frag_abs_threshold
    frag_ratio_threshold = snakemake.params.frag_ratio_threshold
    frag_short_threshold = snakemake.params.frag_short_threshold
    cosmic_min_count = snakemake.params.cosmic_min_count
    vaf_tier1_min = snakemake.params.vaf_tier1_min
    vaf_tier1_max = snakemake.params.vaf_tier1_max
    vaf_tier2_min = snakemake.params.vaf_tier2_min
    vaf_tier2_max = snakemake.params.vaf_tier2_max
    vaf_clonal_tolerance = snakemake.params.vaf_clonal_tolerance
    chip_partner_min_score = snakemake.params.chip_partner_min_score

    vcf_in = pysam.VariantFile(snakemake.input.vcf)
    bam = pysam.AlignmentFile(snakemake.input.bam, "rb")

    hotspot_tabix = (
        pysam.TabixFile(snakemake.input.chip_hotspot_vcf)
        if snakemake.input.chip_hotspot_vcf else None
    )
    cosmic_tabix = (
        pysam.TabixFile(snakemake.input.cosmic_hemato_vcf)
        if snakemake.input.cosmic_hemato_vcf else None
    )

    for id_, number, type_, desc in NEW_INFO:
        vcf_in.header.add_meta(
            "INFO",
            items=[("ID", id_), ("Number", number), ("Type", type_), ("Description", desc)],
        )

    symbol_idx = get_vep_field_index(vcf_in.header, "SYMBOL")
    consequence_idx = get_vep_field_index(vcf_in.header, "Consequence")

    # --- Pass 1: annotate all records, collect CHIP-gene variants for same-clone check ---
    records = []  # (rec, is_non_ref, raw_vaf, chip_gene, in_clonal_context)

    for rec in vcf_in:
        if rec.alts and all(a == "<NON_REF>" for a in rec.alts):
            records.append((rec, True, None, 0, False, 0))
            continue

        raw_vaf = get_af(rec)
        in_clonal_context = is_in_gene_set(rec, symbol_idx, clonal_context_genes)
        chip_gene = chip_gene_flag(rec, symbol_idx, major_genes, minor_genes)
        chip_vaf = chip_vaf_flag(rec, vaf_tier1_min, vaf_tier1_max, vaf_tier2_min, vaf_tier2_max)
        chip_hotspot = chip_hotspot_flag(rec, hotspot_tabix)
        chip_cosmic = chip_cosmic_hemato_flag(rec, cosmic_tabix, symbol_idx, cosmic_exclude_genes, cosmic_min_count)
        chip_consequence = chip_consequence_flag(rec, symbol_idx, consequence_idx, lof_genes, missense_genes)

        prelim_score = (
            chip_gene
            + chip_vaf
            + min(chip_hotspot * 4, 4)
            + min(chip_cosmic, 4)
            + chip_consequence
        )
        # Run pileup if prelim_score is high enough OR variant is in a clonal_context_gene
        # (e.g. TP53) so CHIP_FRAG is always available for same-clone candidates.
        if prelim_score >= 3 or in_clonal_context:
            chip_frag, frag_alt, frag_ref = chip_frag_and_sbs_flags(
                rec, bam, min_alt_reads, frag_abs_threshold, frag_ratio_threshold, frag_short_threshold
            )
        else:
            chip_frag, frag_alt, frag_ref = None, None, None

        score = prelim_score + (chip_frag if chip_frag is not None else 0)

        rec.info["CHIP_GENE"] = chip_gene
        rec.info["CHIP_VAF"] = chip_vaf
        rec.info["CHIP_FRAG"] = chip_frag
        rec.info["CHIP_FRAG_ALT"] = frag_alt
        rec.info["CHIP_FRAG_REF"] = frag_ref
        rec.info["CHIP_HOTSPOT"] = chip_hotspot
        rec.info["CHIP_COSMIC_HEMATO"] = chip_cosmic
        rec.info["CHIP_CONSEQUENCE"] = chip_consequence
        rec.info["CHIP_SAME_CLONE"] = 0  # updated in pass 2
        rec.info["CHIP_SCORE"] = score   # updated in pass 2

        records.append((rec, False, raw_vaf, chip_gene, in_clonal_context, score))

    bam.close()

    # Build partner list: only confirmed CHIP_candidates (score >= chip_partner_min_score)
    # with CHIP_GENE > 0 and a valid VAF. Using only high-confidence CHIP variants as
    # evidence prevents noise-floor variants from triggering the same-clone flag.
    chip_partners = [
        (vaf, gene) for _, is_nr, vaf, gene, _, sc in records
        if not is_nr and gene > 0 and sc >= chip_partner_min_score
        and vaf is not None and vaf > 0
    ]

    # --- Pass 2: compute same-clone flag and write output ---
    with pysam.VariantFile(snakemake.output.vcf, "w", header=vcf_in.header) as vcf_out:
        for rec, is_non_ref, raw_vaf, chip_gene, in_clonal_context, _ in records:
            if is_non_ref:
                vcf_out.write(rec)
                continue

            same_clone = chip_same_clone_flag(raw_vaf, chip_partners, vaf_clonal_tolerance) \
                if in_clonal_context else 0
            rec.info["CHIP_SAME_CLONE"] = same_clone
            rec.info["CHIP_SCORE"] = rec.info["CHIP_SCORE"] + same_clone * 3

            vcf_out.write(rec)


if "snakemake" in dir():
    main(snakemake)

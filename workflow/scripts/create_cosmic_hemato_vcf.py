"""
Convert COSMIC CosmicMutantExport.tsv(.gz) to a bgzip+tabix-ready VCF
containing only haematopoietic_and_lymphoid_tissue SNV entries.

Usage:
    python create_cosmic_hemato_vcf.py CosmicMutantExport.tsv.gz cosmic_hemato.vcf

Then:
    bcftools sort -O z -o cosmic_hemato.vcf.gz cosmic_hemato.vcf
    tabix -p vcf cosmic_hemato.vcf.gz

Only substitution SNVs are included (HGVSG pattern: N:g.POSREF>ALT).
Indels are skipped — they require reference genome lookup for VCF REF bases.
Chromosome names are converted from COSMIC style (9, X, MT) to hg19 style
(chr9, chrX, chrM).
"""

import gzip
import re
import sys
from collections import defaultdict

HEMATO_SITE = "haematopoietic_and_lymphoid_tissue"
SNV_RE = re.compile(r"^(\w+):g\.(\d+)([ACGT])>([ACGT])$")

# Map COSMIC chromosome names to hg19 chr-prefixed names
def normalise_chrom(chrom):
    chrom = chrom.strip()
    if chrom == "MT":
        return "chrM"
    return "chr" + chrom


def parse_hgvsg_snv(hgvsg):
    m = SNV_RE.match(hgvsg.strip())
    if not m:
        return None
    chrom = normalise_chrom(m.group(1))
    pos = int(m.group(2))
    ref = m.group(3)
    alt = m.group(4)
    return chrom, pos, ref, alt


def chrom_sort_key(chrom):
    c = chrom.replace("chr", "")
    if c.isdigit():
        return (0, int(c))
    return (1, c)


def read_fai_contigs(fai_path):
    """Return list of '##contig=<ID=...,length=...>' lines from a .fai file."""
    contigs = []
    with open(fai_path) as fh:
        for line in fh:
            parts = line.split("\t")
            if len(parts) >= 2:
                contigs.append(f"##contig=<ID={parts[0]},length={parts[1].strip()}>\n")
    return contigs


def main(infile, outfile, fai_path=None):
    # variant key -> (sample count, first COSV ID seen)
    counts = defaultdict(int)
    cosv_ids = {}

    opener = gzip.open if infile.endswith(".gz") else open
    with opener(infile, "rt", encoding="utf-8", errors="replace") as fh:
        header = next(fh).rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}

        idx_site = col["Primary site"]
        idx_cosv = col["GENOMIC_MUTATION_ID"]
        idx_hgvsg = col["HGVSG"]

        skipped_indels = 0
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= idx_hgvsg:
                continue
            if parts[idx_site] != HEMATO_SITE:
                continue

            hgvsg = parts[idx_hgvsg]
            parsed = parse_hgvsg_snv(hgvsg)
            if parsed is None:
                skipped_indels += 1
                continue

            chrom, pos, ref, alt = parsed
            key = (chrom, pos, ref, alt)
            counts[key] += 1
            if key not in cosv_ids:
                cosv_ids[key] = parts[idx_cosv]

    print(
        f"Parsed {len(counts)} unique SNV positions "
        f"({skipped_indels} non-SNV entries skipped)",
        file=sys.stderr,
    )

    sorted_keys = sorted(
        counts.keys(),
        key=lambda k: (chrom_sort_key(k[0]), k[1]),
    )

    contig_lines = read_fai_contigs(fai_path) if fai_path else []

    with open(outfile, "w") as out:
        out.write("##fileformat=VCFv4.1\n")
        out.write("##source=COSMIC_v97_haematopoietic_and_lymphoid_tissue_SNVs\n")
        for cl in contig_lines:
            out.write(cl)
        out.write(
            '##INFO=<ID=CNT,Number=1,Type=Integer,'
            'Description="Number of haematological tumour samples with this '
            'mutation in COSMIC">\n'
        )
        out.write(
            '##INFO=<ID=COSV,Number=1,Type=String,'
            'Description="COSMIC genomic mutation ID">\n'
        )
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for key in sorted_keys:
            chrom, pos, ref, alt = key
            cnt = counts[key]
            cosv = cosv_ids[key]
            out.write(f"{chrom}\t{pos}\t{cosv}\t{ref}\t{alt}\t.\t.\tCNT={cnt};COSV={cosv}\n")


if __name__ == "__main__":
    if len(sys.argv) not in (3, 4):
        sys.exit(f"Usage: {sys.argv[0]} CosmicMutantExport.tsv.gz output.vcf [reference.fasta.fai]")
    main(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) == 4 else None)

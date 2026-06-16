__author__ = "Jonas A"
__copyright__ = "Copyright 2025, Jonas A"
__email__ = "jonas.almlof@igp.uu.se"
__license__ = "GPL-3"


rule chip_annotation:
    input:
        vcf="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.artifact_annotated.hotspot_annotated.background_annotated.include.exon.vcf.gz",
        bam=get_deduplication_bam_input,
        bai=lambda wildcards: get_deduplication_bam_input(wildcards) + ".bai",
        chip_genes_file=config.get("chip_annotation", {}).get("chip_genes_file", ""),
        chip_hotspot_vcf=config.get("chip_annotation", {}).get("chip_hotspot_vcf", []),
        chip_hotspot_vcf_tbi=config.get("chip_annotation", {}).get("chip_hotspot_vcf", "") + ".tbi"
        if config.get("chip_annotation", {}).get("chip_hotspot_vcf", "")
        else [],
        cosmic_hemato_vcf=config.get("chip_annotation", {}).get("cosmic_hemato_vcf", []),
        cosmic_hemato_vcf_tbi=config.get("chip_annotation", {}).get("cosmic_hemato_vcf", "") + ".tbi"
        if config.get("chip_annotation", {}).get("cosmic_hemato_vcf", "")
        else [],
    output:
        vcf="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.artifact_annotated.hotspot_annotated.background_annotated.include.exon.chip_annotated.vcf",
    params:
        min_alt_reads=config.get("chip_annotation", {}).get("min_alt_reads", 5),
        frag_abs_threshold=config.get("chip_annotation", {}).get("frag_abs_threshold", 150),
        frag_ratio_threshold=config.get("chip_annotation", {}).get("frag_ratio_threshold", 1.1),
        frag_short_threshold=config.get("chip_annotation", {}).get("frag_short_threshold", 145),
        cosmic_min_count=config.get("chip_annotation", {}).get("cosmic_min_count", 20),
        vaf_tier1_min=config.get("chip_annotation", {}).get("vaf_tier1_min", 0.0),
        vaf_tier1_max=config.get("chip_annotation", {}).get("vaf_tier1_max", 0.10),
        vaf_tier2_min=config.get("chip_annotation", {}).get("vaf_tier2_min", 0.01),
        vaf_tier2_max=config.get("chip_annotation", {}).get("vaf_tier2_max", 0.03),
        vaf_clonal_tolerance=config.get("chip_annotation", {}).get("vaf_clonal_tolerance", 0.5),
        chip_partner_min_score=config.get("chip_annotation", {}).get("chip_partner_min_score", 5),
    log:
        "snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.chip_annotated.vcf.log",
    benchmark:
        repeat(
            "snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.chip_annotated.vcf.benchmark.tsv",
            config.get("chip_annotation", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("chip_annotation", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("chip_annotation", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("chip_annotation", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("chip_annotation", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("chip_annotation", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("chip_annotation", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("chip_annotation", {}).get("container", config["default_container"])
    wildcard_constraints:
        type="T|N",
    script:
        "../scripts/chip_annotation.py"


use rule whatshap_phase from snv_indels as snv_indels_whatshap_phase_chip_annotated with:
    input:
        bam=lambda wildcards: get_deduplication_bam_input(wildcards),
        bai=lambda wildcards: "%s.bai" % get_deduplication_bam_input(wildcards),
        fasta=config.get("reference", {}).get("fasta", ""),
        vcf="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.artifact_annotated.hotspot_annotated.background_annotated.include.exon.chip_annotated.filter.{tag}.vcf",
    output:
        vcf=temp(
            "snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.artifact_annotated.hotspot_annotated.background_annotated.include.exon.chip_annotated.filter.{tag}.phased.vcf"
        ),
    log:
        "snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.artifact_annotated.hotspot_annotated.background_annotated.include.exon.chip_annotated.filter.{tag}.phased.vcf.log",
    benchmark:
        repeat(
            "snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.artifact_annotated.hotspot_annotated.background_annotated.include.exon.chip_annotated.filter.{tag}.phased.vcf.benchmark.tsv",
            config.get("whatshap_phase", {}).get("benchmark_repeats", 1),
        )

__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2024, Jonas Almlöf"
__email__ = "jonas.almlof@sclifelab.uu.se"
__license__ = "GPL-3"


rule estimate_ctdna_fraction:
    input:
        germline_vcf="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.filter.germline.exclude.blacklist.vcf.gz",
        germline_vcf_tabix="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.filter.germline.exclude.blacklist.vcf.gz.tbi",
        segments="cnv_sv/jumble_vcf/{sample}_{type}.pathology_purecn.vcf",
        vcf="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.artifact_annotated.hotspot_annotated.background_annotated.include.exon.filter.snv_hard_filter_umi.codon_snvs.sorted.vep_annotated.qci.vcf",
    output:
        ctDNA_fraction=temp("cnv_sv/estimate_ctdna_fraction/{sample}_{type}.ctDNA_fraction.tsv"),
        ctDNA_fraction_info=temp("cnv_sv/estimate_ctdna_fraction/{sample}_{type}.ctDNA_fraction_info.tsv"),
    params:
        artifact_limit=config.get("estimate_ctdna_fraction", {}).get("artifact_limit", 0),
        callers=config.get("estimate_ctdna_fraction", {}).get("callers", ["vardict"]),
        chip_genes=config.get("estimate_ctdna_fraction", {}).get("chip_genes", ["DNMT3A", "TET2", "ASXL1"]),
        excluded_consequences=config.get("estimate_ctdna_fraction", {}).get(
            "excluded_consequences",
            ["synonymous_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "non_coding_transcript_exon_variant"],
        ),
        max_af=config.get("estimate_ctdna_fraction", {}).get("max_af", 0.4),
        max_gnomad_af=config.get("estimate_ctdna_fraction", {}).get("max_gnomad_af", 0.0002),
        max_msi=config.get("estimate_ctdna_fraction", {}).get("max_msi", 4),
        max_nm=config.get("estimate_ctdna_fraction", {}).get("max_nm", 3.0),
        max_odd_ratio=config.get("estimate_ctdna_fraction", {}).get("max_odd_ratio", 1.5),
        max_panel_median=config.get("estimate_ctdna_fraction", {}).get("max_panel_median", 0.002),
        min_germline_af=config.get("estimate_ctdna_fraction", {}).get("min_germline_af", 0.1),
        min_mq=config.get("estimate_ctdna_fraction", {}).get("min_mq", 40),
        min_nr_SNPs_per_segment=config.get("estimate_ctdna_fraction", {}).get("min_nr_SNPs_per_segment", 35),
        min_position_nr_sd=config.get("estimate_ctdna_fraction", {}).get("min_position_nr_sd", 10),
        min_pmean=config.get("estimate_ctdna_fraction", {}).get("min_pmean", 25),
        min_qual=config.get("estimate_ctdna_fraction", {}).get("min_qual", 40),
        min_sbf=config.get("estimate_ctdna_fraction", {}).get("min_sbf", 0.05),
        min_segment_length=config.get("estimate_ctdna_fraction", {}).get("min_segment_length", 10000000),
        min_sn=config.get("estimate_ctdna_fraction", {}).get("min_sn", 50),
        vaf_baseline=config.get("estimate_ctdna_fraction", {}).get("vaf_baseline", 0.48),
    log:
        "twist_solid/estimate_ctdna_fraction/{sample}_{type}.ctDNA_tc.tsv.log",
    benchmark:
        repeat(
            "twist_solid/estimate_ctdna_fraction/{sample}_{type}.ctDNA_tc.tsv.benchmark.tsv",
            config.get("estimate_ctdna_fraction", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("estimate_ctdna_fraction", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("estimate_ctdna_fraction", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("estimate_ctdna_fraction", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("estimate_ctdna_fraction", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("estimate_ctdna_fraction", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("estimate_ctdna_fraction", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("estimate_ctdna_fraction", {}).get("container", config["default_container"])
    message:
        "{rule}: estimate ctdna fraction based on CNV and SNV data into {output.ctDNA_fraction}"
    script:
        "../scripts/estimate_ctDNA_fraction.py"

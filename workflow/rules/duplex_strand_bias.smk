rule duplex_strand_bias:
    input:
        vcf="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.artifact_annotated.hotspot_annotated.background_annotated.vcf",
        bam=get_deduplication_bam_input,
        bai=get_deduplication_bam_input_bai,
    output:
        vcf=temp(
            "snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.artifact_annotated."
            "hotspot_annotated.background_annotated.duplex_annotated.vcf"
        ),
    params:
        min_per_strand_reads=config.get("duplex_strand_bias", {}).get("min_per_strand_reads", 1),
        min_total_reads=config.get("duplex_strand_bias", {}).get("min_total_reads", 3),
    wildcard_constraints:
        type="T|N",
    log:
        "snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.artifact_annotated."
        "hotspot_annotated.background_annotated.duplex_annotated.vcf.log",
    benchmark:
        repeat(
            "snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.artifact_annotated."
            "hotspot_annotated.background_annotated.duplex_annotated.vcf.benchmark.tsv",
            config.get("duplex_strand_bias", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("duplex_strand_bias", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("duplex_strand_bias", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("duplex_strand_bias", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("duplex_strand_bias", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("duplex_strand_bias", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("duplex_strand_bias", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("duplex_strand_bias", {}).get("container", config["default_container"])
    message:
        "{rule}: annotate strand bias stats excluding fgbio duplex-canonicalized reads in {input.vcf}"
    script:
        "../scripts/duplex_strand_bias.py"

rule jumble_gis_score:
    input:
        gis_csv="cnv_sv/jumble_run/{sample}_{type}/{sample}_{type}.jumble_gis.csv",
        purecn_csv="cnv_sv/purecn/{sample}_{type}.csv",
        ichorcna_params="cnv_sv/ichorcna_offtarget_run/{sample}_{type}/{sample}_{type}.params.txt",
    output:
        gis_score=temp("cnv_sv/jumble_gis_score/{sample}_{type}.{tc_method}.predicted_gis.txt"),
    params:
        pathology_tc=lambda wildcards: get_sample(samples, wildcards).get("tumor_content"),
        hrd_cutoff=config.get("jumble_gis_score", {}).get("hrd_cutoff", 42),
        review_margin=config.get("jumble_gis_score", {}).get("review_margin", 5),
        dna_fraction_review_margin=config.get("jumble_gis_score", {}).get("dna_fraction_review_margin", 0.1),
    log:
        "cnv_sv/jumble_gis_score/{sample}_{type}.{tc_method}.predicted_gis.log",
    benchmark:
        repeat(
            "cnv_sv/jumble_gis_score/{sample}_{type}.{tc_method}.predicted_gis.benchmark.tsv",
            config.get("jumble_gis_score", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("jumble_gis_score", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("jumble_gis_score", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("jumble_gis_score", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("jumble_gis_score", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("jumble_gis_score", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("jumble_gis_score", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("jumble_gis_score", {}).get("container", config["default_container"])
    message:
        "{rule}: Extract predicted GIS at the pathology, PureCN and ichorCNA tumor DNA fractions for {wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/extract_gis_score.py"

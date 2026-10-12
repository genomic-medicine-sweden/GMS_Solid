rule ichorcna_purity:
    input:
        params_txt="cnv_sv/ichorcna_offtarget_run/{sample}_{type}/{sample}_{type}.params.txt",
    output:
        purity=temp("cnv_sv/ichorcna_purity/{sample}_{type}.ichorcna_purity.tsv"),
    log:
        "cnv_sv/ichorcna_purity/{sample}_{type}.ichorcna_purity.tsv.log",
    benchmark:
        repeat(
            "cnv_sv/ichorcna_purity/{sample}_{type}.ichorcna_purity.tsv.benchmark.tsv",
            config.get("ichorcna_purity", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("ichorcna_purity", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("ichorcna_purity", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("ichorcna_purity", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("ichorcna_purity", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("ichorcna_purity", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("ichorcna_purity", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("ichorcna_purity", {}).get("container", config["default_container"])
    message:
        "{rule}: Summarise the selected ichorCNA solution with predicted gender, subclone fraction and noise for {wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/ichorcna_purity.py"

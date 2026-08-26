
rule somalier_best_match_report:
    input:
        pairs="qc/somalier_ungrouped/somalier_relate.pairs.tsv",
    output:
        report=temp("qc/somalier_ungrouped/somalier_best_match.tsv"),
    params:
        extra=config.get("somalier_best_match_report", {}).get("extra", ""),
        match_cutoff=config.get("somalier_best_match_report", {}).get("match_cutoff", 0.7),
    log:
        "qc/somalier_ungrouped/somalier_best_match.tsv.log",
    benchmark:
        repeat(
            "qc/somalier_ungrouped/somalier_best_match.tsv.benchmark.tsv",
            config.get("somalier_best_match_report", {}).get("benchmark_repeats", 1),
        )
    container:
        config.get("somalier_best_match_report", {}).get("container", config["default_container"])
    threads: config.get("somalier_best_match_report", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("somalier_best_match_report", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("somalier_best_match_report", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("somalier_best_match_report", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("somalier_best_match_report", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("somalier_best_match_report", {}).get("time", config["default_resources"]["time"]),
    message:
        "{rule}: generating best match report from somalier pairs"
    script:
        "../scripts/somalier_best_match.py"


rule somalier_dna_dna_report:
    input:
        pairs="qc/somalier_ungrouped/somalier_relate.pairs.tsv",
    output:
        report=temp("qc/somalier_ungrouped/somalier_dna_dna_match.tsv"),
    params:
        extra=config.get("somalier_dna_dna_report", {}).get("extra", ""),
        match_cutoff=config.get("somalier_dna_dna_report", {}).get("match_cutoff", 0.6),
        ibs0_cutoff=config.get("somalier_dna_dna_report", {}).get("ibs0_cutoff", 2),
    log:
        "qc/somalier_ungrouped/somalier_dna_dna_match.tsv.log",
    benchmark:
        repeat(
            "qc/somalier_ungrouped/somalier_dna_dna_match.tsv.benchmark.tsv",
            config.get("somalier_dna_dna_report", {}).get("benchmark_repeats", 1),
        )
    container:
        config.get("somalier_dna_dna_report", {}).get("container", config["default_container"])
    threads: config.get("somalier_dna_dna_report", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("somalier_dna_dna_report", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("somalier_dna_dna_report", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("somalier_dna_dna_report", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("somalier_dna_dna_report", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("somalier_dna_dna_report", {}).get("time", config["default_resources"]["time"]),
    message:
        "{rule}: generating DNA-DNA best match report from somalier pairs"
    script:
        "../scripts/somalier_dna_dna_match.py"


rule somalier_dna_dna_matches_mqc:
    input:
        report="qc/somalier_ungrouped/somalier_dna_dna_match.tsv",
        mqc_config=config.get("somalier_dna_dna_matches_mqc", {}).get("mqc_config", ""),
    output:
        mqc=temp("qc/somalier_ungrouped/somalier_dna_dna_matches_mqc.tsv"),
    log:
        "qc/somalier_ungrouped/somalier_dna_dna_matches_mqc.tsv.log",
    benchmark:
        repeat(
            "qc/somalier_ungrouped/somalier_dna_dna_matches_mqc.tsv.benchmark.tsv",
            config.get("somalier_dna_dna_matches_mqc", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("somalier_dna_dna_matches_mqc", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("somalier_dna_dna_matches_mqc", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("somalier_dna_dna_matches_mqc", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("somalier_dna_dna_matches_mqc", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("somalier_dna_dna_matches_mqc", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("somalier_dna_dna_matches_mqc", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("somalier_dna_dna_matches_mqc", {}).get("container", config["default_container"])
    message:
        "{rule}: creating custom MultiQC content for the somalier DNA-DNA matches table"
    shell:
        """(sed 's/^/# /' {input.mqc_config} > {output.mqc} && cat {input.report} >> {output.mqc}) &> {log}"""

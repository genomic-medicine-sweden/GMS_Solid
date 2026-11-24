__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2026, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule uv_hotspots:
    input:
        bam=get_deduplication_bam_input,
        bai=get_deduplication_bam_input_bai,
        positions=config.get("uv_hotspots", {}).get("positions", ""),
    output:
        positions=temp("biomarker/uv_hotspots/{sample}_{type}.uv_hotspots.positions.tsv"),
        sites=temp("biomarker/uv_hotspots/{sample}_{type}.uv_hotspots.sites.tsv"),
        summary=temp("biomarker/uv_hotspots/{sample}_{type}.uv_hotspots.summary.tsv"),
        call=temp("biomarker/uv_hotspots/{sample}_{type}.uv_hotspots.call.tsv"),
    params:
        sample=lambda wildcards: f"{wildcards.sample}_{wildcards.type}",
        min_alt=config.get("uv_hotspots", {}).get("min_alt", 3),
        min_vaf=config.get("uv_hotspots", {}).get("min_vaf", 0.02),
        uv_sites_hit=config.get("uv_hotspots", {}).get("uv_sites_hit", 3),
        possible_sites_hit=config.get("uv_hotspots", {}).get("possible_sites_hit", 1),
        elevated_ratio=config.get("uv_hotspots", {}).get("elevated_ratio", 1.5),
        elevated_p=config.get("uv_hotspots", {}).get("elevated_p", 0.01),
    log:
        "biomarker/uv_hotspots/{sample}_{type}.uv_hotspots.log",
    benchmark:
        repeat(
            "biomarker/uv_hotspots/{sample}_{type}.uv_hotspots.benchmark.tsv",
            config.get("uv_hotspots", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("uv_hotspots", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("uv_hotspots", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("uv_hotspots", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("uv_hotspots", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("uv_hotspots", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("uv_hotspots", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("uv_hotspots", {}).get("container", config["default_container"])
    message:
        "{rule}: Count C>T at UV hotspot positions (experimental) in {wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/uv_hotspots.py"

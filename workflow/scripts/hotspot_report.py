import os
import tempfile

import yaml
from hydra_genetics.utils.io.hotspot_report import generate_hotspot_report

hotspot_file = snakemake.input.hotspots
vcf_file = snakemake.input.vcf
vcf_file_wo_pick = snakemake.input.vcf_file_wo_pick
gvcf_file = snakemake.input.gvcf

report = snakemake.output.report

levels = snakemake.params.levels
report_config = snakemake.params.report_config
sample_name = snakemake.params.sample_name
chr_translation_file = snakemake.params.chr_translation_file

# Support a single path (str) or a list of paths that are deep-merged in order.
if isinstance(report_config, str):
    config_file = report_config
    tmp_file = None
else:
    merged = {}
    for path in report_config:
        with open(path) as fh:
            data = yaml.safe_load(fh)
        if not merged:
            merged = data
        else:
            merged.setdefault("columns", {}).update(data.get("columns", {}))
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False)
    yaml.dump(merged, tmp)
    tmp.close()
    config_file = tmp.name
    tmp_file = tmp.name

try:
    generate_hotspot_report(sample=sample_name,
                            output=report,
                            levels=levels,
                            hotspot_file=hotspot_file,
                            vcf_file=vcf_file,
                            vcf_file_wo_pick=vcf_file_wo_pick,
                            gvcf_file=gvcf_file,
                            chr_mapping=chr_translation_file,
                            column_yaml_file=config_file)
finally:
    if tmp_file and os.path.exists(tmp_file):
        os.unlink(tmp_file)

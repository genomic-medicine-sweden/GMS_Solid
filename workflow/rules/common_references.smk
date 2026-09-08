__author__ = "Jonas A"
__copyright__ = "Copyright 2021, Jonas A"
__email__ = "jonas.almlof@igp.uu.se"
__license__ = "GPL-3"

import json
import linecache
import pandas
import pandas as pd
import re
import sys
import typing
import yaml
import logging
from snakemake.utils import validate
from snakemake.utils import min_version

from hydra_genetics.utils.misc import get_module_snakefile
from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *

log = logging.getLogger()

min_version("9.0.0")

### Set and validate config file

if not workflow.overwrite_configfiles:
    sys.exit("At least one config file must be passed using --configfile/--configfiles, by command line or a profile!")


validate(config, schema="../schemas/config_references.schema.yaml")
config = load_resources(config, config["resources"])
config = load_resources(config, config["resources_references"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file

samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "type"], drop=False)
validate(units, schema="../schemas/units.schema.yaml")

### Set wildcard constraints


wildcard_constraints:
    sample="|".join(re.escape(s) for s in samples.index),
    unit="N|T|R",


def compile_output_list(wildcards):
    output_files = []
    types = set([unit.type for unit in units.itertuples()])
    for filedef in output_spec["files"]:
        output_files += set(
            [
                filedef["output"].format(
                    sample=sample, type=unit_type, caller=caller, design=config["reference"]["design_bed"].split("/")[-1]
                )
                for sample in get_samples(samples)
                for unit_type in get_unit_types(units, sample)
                if unit_type in set(filedef["types"]).intersection(types)
                for caller in config["bcbio_variation_recall_ensemble"]["callers"]
            ]
        )
    return list(set(output_files))


with open(config["output"]) as output:
    if output.name.endswith(".json"):
        output_spec = json.load(output)
    elif output.name.endswith(".yaml") or output.name.endswith(".yml"):
        output_spec = yaml.safe_load(output)
    else:
        raise ValueError(f"output specification should be JSON or YAML: {output.name}")


class _IdentityLineMap(dict):
    """Line mapping for generated code, whose compiled and source line numbers are the same."""

    def __contains__(self, lineno):
        return True

    def __missing__(self, lineno):
        return lineno


def generate_copy_code(workflow, output_spec):
    code = ""
    for filedef in output_spec["files"]:
        if filedef["input"] is None:
            continue

        input_file = filedef["input"]
        output_file = filedef["output"]
        rule_name = "_copy_{}".format("_".join(re.sub(r"[\"'-.,]", "", filedef["name"].strip().lower()).split()))
        result_file = os.path.basename(filedef["output"])

        if input_file.startswith("{") and input_file.endswith("}"):
            # An input of the form "{some_function}" is an input-function reference. The functions it
            # can name (get_deduplication_bam_input and friends) live in rules/common.smk, which is not
            # in scope here. Skip it: Snakefile_references.smk imports Snakefile as a module, and that
            # generator — which does understand the form, and does have the functions — creates it.
            continue

        if workflow.is_rule(rule_name):
            # The references pipeline includes common_references.smk and then imports Snakefile as a
            # module, so both generators run in one process against the same output spec. Snakemake 9
            # refuses a duplicate rule name (Snakemake 7 silently overwrote it); the second definition
            # would be byte-identical anyway, so keep the first.
            continue

        mem_mb = config.get("_copy", {}).get("mem_mb", config["default_resources"]["mem_mb"])
        mem_per_cpu = config.get("_copy", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"])
        partition = config.get("_copy", {}).get("partition", config["default_resources"]["partition"])
        threads = config.get("_copy", {}).get("threads", config["default_resources"]["threads"])
        time = config.get("_copy", {}).get("time", config["default_resources"]["time"])
        copy_container = config.get("_copy", {}).get("container", config["default_container"])

        code += f'@workflow.rule(name="{rule_name}")\n'
        code += f'@workflow.input("{input_file}")\n'
        code += f'@workflow.output("{output_file}")\n'
        code += f'@workflow.log("logs/{rule_name}_{result_file}.log")\n'
        code += f'@workflow.container("{copy_container}")\n'
        code += f'@workflow.resources(time = "{time}", threads = {threads}, mem_mb = {mem_mb}, mem_per_cpu = {mem_per_cpu}, partition = "{partition}")\n'
        code += '@workflow.shellcmd("cp --preserve=timestamps {input} {output}")\n\n'
        code += "@workflow.run\n"
        code += (
            f"def __rule_{rule_name}(input, output, params, wildcards, threads, resources, log, rule, "
            "conda_env, container_img, singularity_args, use_singularity, env_modules, bench_record, jobid, is_shell, "
            "bench_iteration, cleanup_scripts, shadow_dir, edit_notebook, conda_base_path, basedir, sourcecache_path, "
            "runtime_sourcecache_path, runtime_paths, __is_snakemake_rule_func=True):\n"
            '\tshell ( "(cp --preserve=timestamps {input[0]} {output[0]}) &> {log}" , bench_record=bench_record, bench_iteration=bench_iteration)\n\n'
        )

    source_name = "reference_result_to_copy"

    # See the identical comment in rules/common.smk. The name must differ from the one used there:
    # the references pipeline imports Snakefile as a module, so both generators run in one process
    # and a shared name would clobber the other's linecache entry.
    linecache.cache[source_name] = (len(code), None, code.splitlines(True), source_name)
    workflow.linemaps[source_name] = _IdentityLineMap()

    exec(compile(code, source_name, "exec"), workflow.globals)


generate_copy_code(workflow, output_spec)


def get_reference_bam_input(sample: str, type: str) -> str:
    # Mirrors get_deduplication_bam_input() in rules/common.smk: the reference
    # pipeline must build PoN/background/PureCN/etc. from the same final,
    # deduplicated bam the main pipeline actually reports on, rather than the
    # pre-dedup/pre-consensus alignment/samtools_merge_bam bam (which for FFPE
    # still contains PCR duplicates, and for ctDNA/UMI predates consensus
    # calling entirely).
    if config.get("deduplication") == "umi":
        return "alignment/samtools_merge_bam_umi/%s_%s.bam" % (sample, type)
    return "alignment/samtools_merge_bam_all_final/%s_%s.bam" % (sample, type)


def get_bams(units: pandas.DataFrame, name: str) -> typing.List[str]:
    return get_files(units, name, get_reference_bam_input)


def get_counts(units: pandas.DataFrame, name: str) -> typing.List[str]:
    return get_files(units, name, "references/jumble_count/%s_%s.bam.counts.RDS")


def get_hdf5(units: pandas.DataFrame, name: str) -> typing.List[str]:
    return get_files(units, name, "references/collect_read_counts/%s_%s.counts.hdf5")


def get_wigs(units: pandas.DataFrame, name: str) -> typing.List[str]:
    return get_files(units, name, "references/ichorcna_offtarget_read_counter/%s_%s.wig")


def get_vcfs(units: pandas.DataFrame, name: str) -> typing.List[str]:
    return get_files(units, name, "snv_indels/bcbio_variation_recall_ensemble/%s_%s.ensembled.vep_annotated.vcf.gz")


def get_gvcfs(units: pandas.DataFrame, name: str) -> typing.List[str]:
    return get_files(units, name, "qc/add_mosdepth_coverage_to_gvcf/%s_%s.mosdepth.g.vcf.gz")


def get_cnv_vcfs(units: pandas.DataFrame, name: str) -> typing.List[str]:
    return get_files(units, name, "cnv_sv/svdb_merge/%s_%s.pathology_purecn.merged.vcf")


def get_cnvkit_target(units: pandas.DataFrame, name: str) -> typing.List[str]:
    return get_files(units, name, "%s_%s.targetcoverage.cnn")


def get_cnvkit_antitarget(units: pandas.DataFrame, name: str) -> typing.List[str]:
    return get_files(units, name, "%s_%s.antitargetcoverage.cnn")


def get_files(units: pandas.DataFrame, name: str, string_path):
    types = []
    sample_list = get_samples(samples)
    for i in output_spec["files"]:
        if i["name"] == name:
            types = i["types"]
    to_path = string_path if callable(string_path) else lambda sample, type: string_path % (sample, type)
    data = [to_path(t.sample, t.type) for t in units[units["type"].isin(types)].itertuples() if t.sample in sample_list]
    if not data:
        log.warning(f"No files matching the output files found for rules using name: {name}, {string_path}")
    return set(data)

# References, panel of normals and design files

## Easy setup

### Download data
Use hydra-genetics to setup reference files. Remember to update config/config.data.hg19.yaml and include it when running an analysis.
```bash
# make sure hydra-genetics is available
# make sure that TMPDIR points to a location with a lot of storage, it
# will be required to fetch reference data
export TMPDIR=/PATH_TO_STORAGE
# NextSeq
 hydra-genetics --debug references download -o design_and_ref_files -v config/references/design_files.hg19.yaml -v config/references/nextseq.hg19.pon.yaml -v config/references/references.hg19.yaml

 #NovaSeq, not all files are prepare for novaseq
 hydra-genetics references download -o design_and_ref_files -v config/references/design_files.hg19.yaml -v config/references/novaseq.hg19.pon.yaml -v config/references/references.hg19.yaml

```

### Validate if data requires update

To validate if all design and reference files are up to date the following command can be run, assuming that they are store at the same parent folder.
```bash
# This will make sure that all design and reference files exists and haven't changed
# Warnings for possible file PATH/hydra-genetics and missing tbi files in config can be ignored
hydra-genetics --debug references validate -c config/config.yaml -c config/config.data.hg19.yaml -v config/references/design_files.hg19.yaml -v config/references/nextseq.hg19.pon.yaml -v config/references/references.hg19.yaml  -p ${PATH_TO_design_and_ref_files} 
```

## References overview
The following reference files, panel of normals and design files are needed to run the Twist Solid Pipeline:

| Rule | Config name | File |
|-|-|-|
| reference | <div id="background_db"> background </div> | `background_panel_nextseq_noUmea_27_dp500_af015.tsv` |
| | <div id="artifact_db"> artifacts </div> | `artifact_panel_nextseq_36.tsv` |
| | <div id="reference_fasta">fasta</div> | `hg19.with.mt.fasta` |
| | fasta_rna | `GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa` |
| | dict | `hg19.with.mt.dict` |
| | fai | `hg19.with.mt.fai` |
| | <div id="design_bed">design_bed</div> | `pool1_pool2.sort.merged.padded20.cnv200.hg19.split_fusion_genes.reannotated.230222.bed` |
| | design_intervals | `pool1_pool2.sort.merged.padded20.cnv200.hg19.split_fusion_genes.MUC6_31_rm.exon_only.reannotated.210608.interval_list` |
| | <div id="ref_gatk_intervals">design_intervals_gatk_cnv</div> | `pool1_pool2_nochr_3c.sort.merged.padded20.cnv400.hg19.210311.met.annotated.bed.preprocessed.interval_list` |
| | <div id="design_bed_rna">design_bed_rna</div> | `Twist_RNA_Design5.annotated.bed` |
|_ _| design_intervals_rna | `Twist_RNA_Design5.annotated.20230630.interval_list` |
| <div id="arriba_reference">arriba</div> | assembly | `hg19.with.mt.fasta` |
| | <div id="arriba_blacklist">blacklist</div> | `arriba/arriba_v2.3.0/database/blacklist_hg19_hs37d5_GRCh37_v2.3.0.tsv.gz` |
| | <div id="arriba_gtf">gtf</div> | `hg19.refGene.gtf` |
| | <div id="arriba_extra_gff3_">extra</div> | `arriba/arriba_v2.3.0/database/protein_domains_hg19_hs37d5_GRCh37_v2.3.0.gff3` |
|_ _| <div id="arriba_extra_tsv">extra | `arriba/arriba_v2.3.0/database/known_fusions_hg19_hs37d5_GRCh37_v2.3.0.tsv.gz` |
| <div id="arriba_draw_fusion_cytobands">arriba_draw_fusion</div> | cytobands | `arriba/arriba_v2.3.0/database/cytobands_hg19_hs37d5_GRCh37_v2.3.0.tsv` |
| | <div id="arriba_draw_fusion_gtf">gtf</div> | `hg19.refGene.gtf` |
|_ _| <div id="arriba_draw_fusion_protein_domains">protein_domains</div> | `arriba/arriba_v2.3.0/database/protein_domains_hg19_hs37d5_GRCh37_v2.3.0.gff3` |
| <div id="cnv_amp_genes">annotate_cnv</div> | cnv_amp_genes | `cnv_amp_genes_240307.bed` |
|_ _| <div id="cnv_loh_genes">cnv_loh_genes</div> | `cnv_loh_genes.bed` |
| bcftools_annotate | annotation_db | `small_exac_common_3.hg19.vcf.gz` |
| <div id="bcftools_filter">bcftools_filter_include_region</div> | exon | `pool1_pool2.sort.merged.padded20.cnv200.hg19.split_fusion_genes.MUC6_31_rm.exon_only.reannotated.210608.bed` |
| <div id="bcftools_filter_exclude_region">bcftools_filter_exclude_region</div> | blacklist | `cnvkit_germline_blacklist_20221221.bed` |
| <div id="bwa_me_ref">bwa_mem | amb </div> | `hg19.with.mt.amb` |
| | ann | `hg19.with.mt.ann` |
| | bwt | `hg19.with.mt.bwt` |
| | pac | `hg19.with.mt.pac` |
|_ _| sa | `hg19.with.mt.sa` |
| <div id="call_small_cnv_deletions">call_small_cnv_deletions</div> | regions_file | `cnv_deletion_genes_240618.tsv` |
| <div id="call_small_cnv_amplifications">call_small_cnv_amplifications</div> | regions_file | `cnv_amplification_genes_240307.tsv` |
| <div id="cnvkit_ref">cnvkit_batch</div> | normal_reference | `cnvkit_nextseq_36.cnn` |
| <div id="normal_reference_hrd">cnvkit_batch_hrd</div> | normal_reference_hrd | `cnvkit_nextseq_27_HRD.cnn` |
| <div id="exon_skipping">exon_skipping</div> | design_bed | `Twist_RNA_Design5.annotated.bed` |
| <div id="fusioncatcher">fusioncatcher</div> | genome_path | `human_v102/` |
| <div id="gatk_collect_allelic_counts">gatk_collect_allelic_counts</div> | SNP_interval | `gnomad_SNP_0.001_target.annotated.interval_list` |
| <div id="gatk_denoise_read_counts_pon">gatk_denoise_read_counts</div> | normal_reference | `gatk_cnv_nextseq_36.hdf5` |
| <div id="gatk_get_pileup_summaries_sites">gatk_get_pileup_summaries</div> | sites | `small_exac_common_3.hg19.vcf.gz` |
|_ _| variants | `small_exac_common_3.hg19.vcf.gz` |
| <div id="genefuse_transcripts">gene_fuse</div> | genes | `GMS560_fusion_w_pool2.hg19.221117.csv` |
|_ _| <div id="genefuse_fasta">fasta</div> | `hg19.with.mt.fasta` |
| <div id="fuseq_wes_json">FuSeq_WES</div> | transcript annotation | `UCSC_hg19_wes_contigSize3000_bigLen130000_r100.json` |
| | <div id="fuseq_wes_sqlite">transcript database</div> | `UCSC_hg19_wes_contigSize3000_bigLen130000_r100.sqlite` |
| | <div id="fuseq_wes_fusion_db">fusion database</div> | `Mitelman_fusiondb.RData` |
|_ _| <div id="fuseq_wes_paralog_db">paralog database</div> | `ensmbl_paralogs_grch37.RData` |
| <div id="filter_report_fuseq_wes">filter_report_fuseq_wes</div> | transcript annotation | `hg19.refGene.gtf` |
| | <div id="fuseq_wes_white_list">gene white list</div> | `fuseq_wes_gene_white_list.txt` |
| | <div id="gene_fusion_black_list">fusion gene pair black list</div> | `false_positive_fusion_pairs.txt` |
|_ _| <div id="fuseq_wes_transcript_black_list">transcript black list</div> | `fuseq_wes_transcript_black_list.txt` |
| <div id="hotspot_file">hotspot_annotation</div> | hotspots | `Hotspots_combined_regions_nodups.csv` |
| <div id="hotspot_report">hotspot_report</div> | hotspot_mutations | `Hotspots_combined_regions_nodups.csv` |
| <div id="jumble_run">jumble_run</div> | normal_reference | `jumble.combined.filtered.50.PoN.hg19.RDS` |
| <div id="manta_design_bed">manta_config_t</div> | extra | `pool1_pool2.sort.merged.padded20.cnv200.hg19.split_fusion_genes.210608.bed.gz` |
| <div id="msisensor_pro_pon">msisensor_pro</div> | PoN | `Msisensor_pro_reference_nextseq_36.list_baseline` |
| <div id="purecn_estimation_mapping_pon">purecn</div> | extra | `mapping_bias_nextseq_27_hg19.rds` |
| | <div id="purecn_estimation_normaldb">normaldb</div> | `normalDB_nextseq_27_hg19.rds` |
|_ _| <div id="purecn_estimation_intervals">intervals</div> | `targets_twist-gms-st_hg19_25000_intervals.txt` |
| <div id="purecn_coverage_intervals">purecn_coverage</div> | intervals | `targets_twist-gms-st_hg19_25000_intervals.txt` |
| <div id="report_fusions">report_fusions</div> | annotation_bed | `Twist_RNA_fusionpartners.bed` |
| <div id="genefuse_filter_fusions">report_gene_fuse</div> | filter_fusions | `filter_fusions_20230214.csv` |
| <div id="star_genome_index">star</div> | genome_index | `v2.7.10a_hg19/` |
|_ _| <div id="star_genome_extra">extra</div> | `hg19.refGene.gtf` |
| <div id="star_fusion">star_fusion</div> | genome_path | `GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/` |
| <div id="svdb_query">svdb_query</div> | db_string | `all_TN_292_svdb_0.8_20220505.vcf` |
| <div id="vep_cache">vep</div> | vep_cache | `VEP/` |

## Downloadable reference files

### Fasta reference genome
hg19 with mitochondria but without HLA and without decoys
```bash
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/*'
gunzip *.fa.gz
cat *.fa > hg19.with.mt.fasta
```

### BWA indexes
```bash
bwa index hg19.with.mt.fasta
```

### VEP v105
```bash
wget http://ftp.ensembl.org/pub/grch37/release-105/variation/vep/homo_sapiens_refseq_vep_105_GRCh37.tar.gz
```

### GNOMAD common SNPs
Used by GATK pileup (CNV calling) and GATK contamination
```bash
gsutil cp gs://gatk-best-practices/somatic-b37/small_exac_common_3.vcf* .
```
Adapt file to hg19 (add chr at all lines and in header)

### Mappability file
Used in CNVkit PoN creation
```bash
wget https://raw.githubusercontent.com/etal/cnvkit/master/data/access-5k-mappable.hg19.bed
```

### Arriba v2.3.0
```bash
singularity exec -B /path/to/references:/references docker://hydragenetics/arriba:2.3.0 download_references.sh hs37d5+RefSeq
wget https://github.com/suhrig/arriba/releases/download/v2.3.0/arriba_v2.3.0.tar.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.refGene.gtf.gz
```

### Star genome index
```bash
singularity exec docker://hydragenetics/star:2.7.10a STAR --runThreadN 8 --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles Human_genome.fasta
```

### Fusioncather v102
```bash
wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.aa
wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.ab
wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.ac
wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.ad
cat human_v102.tar.gz.* | tar xz
ln -s human_v102 current
```

### Star-fusion
```bash
wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh37_gencode_v19_CTAT_lib_Mar012021.plug-n-play.tar.gz
```

## Panel of normals

Files needed by Twist Solid that are generated by the reference pipeline are listed below. Many of the panel of normals are available to download from the Uppsala Owncloud solution but should preferably be generated on in-house data.  
<br />
**Result files**

* `result/cnvkit.PoN.cnn`
* `result/gatk_cnv_panel_of_normal.hdf5`
* `result/design.preprocessed.interval_list`
* `result/jumble.PoN.RDS`
* `result/Msisensor_pro_reference.list_baseline`
* `result/background_panel.tsv`
* `result/artifact_panel.tsv`
* `result/svdb_cnv.vcf`
* `result/mapping_bias.rds`
* `result/purecn_normal_db.rds`
* `result/purecn_targets_intervals.txt`

### Create samples and units
Files and samples used in the generation of the panel of normals are specified in samples.tsv and units.tsv. Required files are listed down below.
Adapt out file specification (`workflow/rules/common_references.smk`) and comment out files that should not be generated.

### Run command
Run the pipeline in the same way as the [standard pipeline](running.md), using a reference specific profile and Snakefile:
```bash
snakemake --profile profiles/uppsala_ref/ -s workflow/Snakefile_references.smk
```
<br />
**Note**  
The `units.tsv` file needs to be adapted depending which panel of normals are created and should contain all the samples needed to create the panel of normals.

### CNVkit

**Reference files**

* design bedfile
* fasta genome reference
* mappability file

### GATK CNV

**Reference files**

* design bedfile
* fasta reference genome
* fasta reference genome dictionary

### MSISensor-pro

**Reference files**

* design bedfile
* fasta reference genome

**Software settings**

| Options | Value | Description |
|-|-|-|
| extra | -c 50 | minimal coverage, recommended for WES: 20; WGS: 15 |

### SVDB

Should be made up of both normal and tumor FFPE samples!

**Software settings**

| Options | Value | Description |
|-|-|-|
| extra | --overlap 0.8 | Overlap used to cluster variants (default 0.8) |

### Artifacts

Based on unfiltered and merged vcf files from normal FFPE samples

### Background

Based on genome vcf files from Mutect2 from normal FFPE samples

**Software settings**

| Options | Value | Description |
|-|-|-|
| min_dp | 500 | Min read depth to be included (default: 500) |
| max_af | 0.015 | Max allele frequency to be included (default: 0.015) |

### PureCN

Made up of bam and vcf files from normal samples

## Pipeline specific files
Premade panel of normals, design files and references can be download using the hydra-genetics tools. Design files can also be downloaded for our [github repo](https://github.com/genomic-medicine-sweden/Twist_Solid_pipeline_files). Check the [config yaml files](https://github.com/genomic-medicine-sweden/Twist_Solid/tree/develop/config) in the Twist_Solid repo for the latest files.
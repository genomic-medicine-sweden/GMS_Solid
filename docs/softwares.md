# Rule specific to the Twist Solid pipeline that are not defined in hydra

## bcftools_id_snps
Trim `.fastq` files by removing adapter sequences and other unwanted sequences. Adapter sequences are specified in `units.tsv` under the adapter column. See further [ID-SNPs info](rna_id_snps.md).


### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__bcftools__bcftools_id_snps#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__bcftools__bcftools_id_snps#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__bcftools_id_snps#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__bcftools_id_snps#

---

## call_small_cnv_amplifications
The CNVkit and GATK CNV caller often miss small amplifications of four exons or smaller. This rule analyses genes of interest and reports small amplifications in these genes when predefined criteria are met. The caller is based on two rolling window averaged over the region and finds the largest difference in log2ratio between these two windows. If the difference is large enough it is reported. See further [call small cnv amplifications info](dna_cnvs.md#small-cnv-amplifications).

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__call_small_cnv_amplifications__call_small_cnv_amplifications#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__call_small_cnv_amplifications__call_small_cnv_amplifications#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__call_small_cnv_amplifications#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__call_small_cnv_amplifications#

---

## call_small_cnv_deletions
The CNVkit and GATK CNV caller often miss small deletions of four exons or smaller. This rule analyses genes of interest and reports small deletions in these genes when predefined criteria are met. The caller is based on two rolling window averaged over the region and finds the largest difference in log2ratio between these two windows. If the difference is large enough it is reported. See further [call small cnv deletions info](dna_cnvs.md#small-cnv-deletions).

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__call_small_cnv_deletions__call_small_cnv_deletions#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__call_small_cnv_deletions__call_small_cnv_deletions#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__call_small_cnv_deletions#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__call_small_cnv_deletions#

---

## cnv_tsv_report
Collect all CNV calls into an excel friendly text file. Adds potential 1p19q calls. Also add a FP flag for variants called by CNVkit that in GATK CNV does not have any signal in that region. See further [CNV tsv report info](dna_cnvs.md#cnv-report).

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__cnv_tsv_report__cnv_tsv_report#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__cnv_tsv_report__cnv_tsv_report#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__cnv_tsv_report#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__cnv_tsv_report#

---

## estimate_ctdna_fraction
Estimate ctDNA fraction based on CNV and SNP information using an in-house script. For CNV based estimation the VAF-values of germline SNPs (the BAF-plot) for each segment is used to make an density calculation. If two peaks are found on opposite sides of 50% allele frequency the separation of the peak is used to calculate the ctDNA fraction. The highest ctDNA fraction is reported. If no CNV segments have peak separation 0% is reported. For SNV based estimation the ctDNA is reported as the highest variant allele frequency SNV times two.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__estimate_ctdna_fraction__estimate_ctdna_fraction#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__estimate_ctdna_fraction__estimate_ctdna_fraction#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__estimate_ctdna_fraction#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__estimate_ctdna_fraction#

---

## exon_skipping
Calls exon skipping events in RNA data. Only reports MET exon 14 skipping and the EGFRvIII variant. See further [exon skipping info](rna_exon_skipping.md).

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__exon_skipping__exon_skipping#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__exon_skipping__exon_skipping#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__exon_skipping#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__exon_skipping#

---

## fix_vcf_ad_for_qci
Produces a QCI compatible vcf file with SNV and INDEL calls were the AD field has been corrected to so that the numbers corresponds with the AF field. QCI uses the AD field to calculate allele frequency. See further [QCI fix info](dna_snv_indels.md#qci-af-correction-of-vcf).

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__fix_vcf_ad_for_qci__fix_vcf_ad_for_qci#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__fix_vcf_ad_for_qci__fix_vcf_ad_for_qci#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__fix_vcf_ad_for_qci#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__fix_vcf_ad_for_qci#

---

## hotspot_report
A excel friendly text file containing information of all hotspot postions as well as all called variants. For all hotspot positions the coverage and alternative allele information is reported as well as additional annotations. The columns reported are configured by the config file `hotspot_report.yaml`. See further [Hotspot report info](dna_qc.md#coverage-and-mutations).

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__hotspot_report__hotspot_report#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__hotspot_report__hotspot_report#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__hotspot_report#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__hotspot_report#


## hotspot_report_aa_translate
Take the output from hotspot report (coverage_and mutations) and translate the amino acid three letter code into one letter code using an in-house script.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__hotspot_report__hotspot_report_aa_translate#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__hotspot_report__hotspot_report_aa_translate#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__hotspot_report_aa_translate#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__hotspot_report_aa_translate#

---

## jumble_gis_score
Extract the predicted Genomic Instability Score (GIS) from Jumble's output for the current sample's tumor content (TC). 

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__jumble_gis__jumble_gis_score#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__jumble_gis__jumble_gis_score#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__jumble_gis_score#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__jumble_gis_score#

---

## house_keeping_gene_coverage
List the average coverage of the housekeeping genes GAPDH, GUSB, OAZ1, and POLR2A. See further [house keeping gene coverage info](rna_qc.md#house-keeping-gene-coverage).

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__house_keeping_gene_coverage__house_keeping_gene_coverage#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__house_keeping_gene_coverage__house_keeping_gene_coverage#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__house_keeping_gene_coverage#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__house_keeping_gene_coverage#

---

## report_fusions
Make a combined report of filtered fusion calls from all RNA fusion callers. Applies caller specific thresholds based on read support and fusion type. Annotates breakpoints to genes and exons when available. See further [RNA fusions report info](rna_fusions.md#fusion-filtering-and-report).

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__report_fusions__report_fusions#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__report_fusions__report_fusions#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__report_fusions#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__report_fusions#

---

## report_gene_fuse
Make a report of filtered fusion calls from the DNA fusion caller Gene Fuse. Applies configurable thresholds for supporting reads and flags noisy genes as well as filters artifact genes. See further [DNA fusions report info](dna_fusions.md).

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__report_gene_fuse__report_gene_fuse#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__report_gene_fuse__report_gene_fuse#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__report_gene_fuse#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__report_gene_fuse#

---

## purecn_modify_vcf
Increases the MQB (mean base quality) value by 5 as the qualities are so bad for our samples.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__purecn_modify_vcf__purecn_modify_vcf#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__purecn_modify_vcf__purecn_modify_vcf#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__purecn_modify_vcf#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__purecn_modify_vcf#

---

## sample_mixup_check
Compare ID-SNPs in the RNA samples to the DNA samples in the same analysis and report sample similarities to be able to discern sample mixups

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__sample_mixup_check__sample_mixup_check#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__sample_mixup_check__sample_mixup_check#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__sample_mixup_check#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__sample_mixup_check#

---

## somalier_best_match_report
Takes the somalier relate pairs output and finds the best match for each sample (DNA and RNA). The script handles the bidirectional nature of the pairs file to ensure every sample is reported with its top relatedness match.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__somalier__somalier_best_match_report#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__somalier__somalier_best_match_report#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__somalier_best_match_report#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__somalier_best_match_report#



---

## somalier_dna_dna_report
Takes the somalier relate pairs output and finds the best DNA-DNA match for each DNA sample. This is used to detect potential sample mixups between DNA samples.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__somalier__somalier_dna_dna_report#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__somalier__somalier_dna_dna_report#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__somalier_dna_dna_report#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__somalier_dna_dna_report#

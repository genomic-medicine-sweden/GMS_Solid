# Result files
The pipeline generates many files, but only specific ones are defined as final outputs. These are configured in `config/output_files_FFPE.yaml` and `config/output_files_ctDNA.yaml`.
{sample} is replaced with the sample name and {type} is replaced with the sample type (T/N/R, where T=Tumor, N=Normal, R=RNA).

## FFPE Pipeline Output
Files defined in `output_files_FFPE.yaml`.

### DNA

#### Alignment
| **File** | **Description** |
|---|---|
| `bam_dna/{sample}_{type}.bam` | BAM |
| `bam_dna/{sample}_{type}.bam.bai` | BAM index |
| `bam_dna/mutect2_indel_bam/{sample}_{type}.bam` | Mutect2 INDEL BAM |
| `bam_dna/mutect2_indel_bam/{sample}_{type}.bam.bai` | Mutect2 INDEL BAM index |

#### SNV and INDELs
| **File** | **Description** |
|---|---|
| `gvcf_dna/{sample}_{type}.mosdepth.g.vcf.gz` | Mosdepth GVCF |
| `results/dna/{sample}_{type}/additional_files/vcf/{sample}_{type}.annotated.vcf.gz` | Vep annotated VCF |
| `results/dna/{sample}_{type}/additional_files/vcf/{sample}_{type}.annotated.exon_only.filter.soft_filter.vcf` | Exon-only soft-filtered VCF |
| `results/dna/{sample}_{type}/additional_files/vcf/{sample}_{type}.annotated.exon_only.filter.hard_filter.vcf` | Exon-only hard-filtered VCF |
| `results/dna/{sample}_{type}/vcf/{sample}_{type}.annotated.exon_only.filter.hard_filter.codon_snv.vcf` | Codon SNV VCF |
| `results/dna/{sample}_{type}/vcf/{sample}_{type}.annotated.exon_only.filter.hard_filter.codon_snv.qci.vcf` | QCI VCF |
| `results/dna/{sample}_{type}/additional_files/vcf/{caller}_{sample}_{type}.vcf.gz` | Caller VCF |
| `results/dna/{sample}_{type}/id_snps/{sample}_{type}.id_snps.vcf` | ID-SNP VCF DNA |

#### QC
| **File** | **Description** |
|---|---|
| `results/dna/{sample}_{type}/additional_files/qc/{sample}_{type}.duplication_metrics.txt` | Duplication metrics |
| `results/dna/{sample}_{type}/additional_files/qc/{sample}_{type}.alignment_summary_metrics.txt` | Alignment summary metrics |
| `results/dna/{sample}_{type}/additional_files/qc/{sample}_{type}.HsMetrics.txt` | HsMetrics |
| `results/dna/{sample}_{type}/additional_files/qc/{sample}_{type}.insert_size_metrics.txt` | Insert size metrics |
| `results/dna/{sample}_{type}/additional_files/qc/{sample}_{type}.samtools-stats.txt` | Samtools stats |
| `results/dna/{sample}_{type}/qc/{sample}_{type}.coverage_and_mutations.tsv` | Coverage and mutations all |
| `results/dna/{sample}_{type}/qc/{sample}_{type}.coverage_and_mutations.ENC.tsv` | Coverage and mutations ENC |
| `results/dna/{sample}_{type}/qc/{sample}_{type}.coverage_and_mutations.URO.tsv` | Coverage and mutations URO |
| `results/dna/{sample}_{type}/qc/{sample}_{type}.coverage_and_mutations.BRC.tsv` | Coverage and mutations BRC |
| `results/dna/{sample}_{type}/qc/{sample}_{type}.coverage_and_mutations.CRC.tsv` | Coverage and mutations CRC |
| `results/dna/{sample}_{type}/qc/{sample}_{type}.coverage_and_mutations.GIS.tsv` | Coverage and mutations GIS |
| `results/dna/{sample}_{type}/qc/{sample}_{type}.coverage_and_mutations.LUN.tsv` | Coverage and mutations LUN |
| `results/dna/{sample}_{type}/qc/{sample}_{type}.coverage_and_mutations.MEL.tsv` | Coverage and mutations MEL |
| `results/dna/{sample}_{type}/qc/{sample}_{type}.coverage_and_mutations.MTC.tsv` | Coverage and mutations MTC |
| `results/dna/{sample}_{type}/additional_files/qc/{sample}_{type}.contamination.table` | Contamination table |
| `results/dna/qc/multiqc_DNA.html` | MultiQC DNA HTML |

#### Reports
| **File** | **Description** |
|---|---|
| `results/dna/{sample}_{type}/{sample}_{type}.general_report.html` | Results report DNA HTML |
| `results/dna/{sample}_{type}/{sample}_{type}.prio.general_report.html` | Results report small panel DNA HTML |

#### Biomarker
| **File** | **Description** |
|---|---|
| `results/dna/{sample}_{type}/biomarker/{sample}_{type}.msisensor_pro.filtered.score.tsv` | MSI Sensor Pro filtered TSV |
| `results/dna/{sample}_{type}/biomarker/{sample}_{type}.msisensor_pro.unfiltered.score.tsv` | MSI Sensor Pro unfiltered TSV |
| `results/dna/{sample}_{type}/biomarker/{sample}_{type}.TMB.txt` | TMB |
| `results/dna/{sample}_{type}/additional_files/biomarker/{sample}_{type}.purecn.scarhrd_cnvkit_score.txt` | HRD score - PureCN |
| `results/dna/{sample}_{type}/additional_files/biomarker/{sample}_{type}.pathology.scarhrd_cnvkit_score.txt` | HRD score - pathology |
| `results/dna/{sample}_{type}/biomarker/{sample}_{type}.pathology_purecn.scarhrd_cnvkit_score.txt` | HRD score - pathology_purecn |
| `results/dna/{sample}_{type}/biomarker/{sample}_{type}.pathology_purecn.predicted_gis.txt` | Predicted GIS score |

#### Fusions
| **File** | **Description** |
|---|---|
| `results/dna/{sample}_{type}/fusion/{sample}_{type}.fuseq_wes.report.csv` | fuseq WES |
| `results/dna/{sample}_{type}/additional_files/fusion/{sample}_{type}.fuseq_wes.unfiltered.results.csv` | fuseq WES unfiltered |
| `results/dna/{sample}_{type}/fusion/{sample}_{type}.juli.filtered.fusions.txt` | JuLI fusions |
| `results/dna/{sample}_{type}/additional_files/fusion/{sample}_{type}.juli.fusions.txt` | JuLI unfiltered fusions |

#### CNV
| **File** | **Description** |
|---|---|
| `results/dna/{sample}_{type}/additional_files/cnv/{sample}_{type}.cnvkit.scatter.png` | CNVkit scatter PNG |
| `results/dna/{sample}_{type}/additional_files/cnv/{sample}_{type}.cnvkit.diagram.pdf` | CNVkit diagram PDF |
| `results/dna/{sample}_{type}/additional_files/cnv/{sample}_{type}.purecn.svdb_query.vcf` | SVDB query - PureCN |
| `results/dna/{sample}_{type}/additional_files/cnv/{sample}_{type}.pathology.svdb_query.vcf` | SVDB query - pathology |
| `results/dna/{sample}_{type}/additional_files/cnv/{sample}_{type}.pathology_purecn.svdb_query.vcf` | SVDB query - pathology_purecn |
| `results/dna/{sample}_{type}/additional_files/cnv/{sample}_{type}.purecn.cnv.html` | CNV HTML report - PureCN |
| `results/dna/{sample}_{type}/additional_files/cnv/{sample}_{type}.pathology.cnv.html` | CNV HTML report - pathology |
| `results/dna/{sample}_{type}/cnv/{sample}_{type}.pathology_purecn.cnv.html` | CNV HTML report - pathology_purecn |
| `results/dna/{sample}_{type}/cnv/{sample}_{type}.purecn_purity_ploidity.csv` | PureCN purity |
| `results/dna/{sample}_{type}/additional_files/cnv/{sample}_{type}.purecn.amp_all_del_all.cnv_report.tsv` | CNV TSV report - PureCN |
| `results/dna/{sample}_{type}/additional_files/cnv/{sample}_{type}.pathology.amp_all_del_all.cnv_report.tsv` | CNV TSV report - pathology |
| `results/dna/{sample}_{type}/additional_files/cnv/{sample}_{type}.pathology_purecn.amp_all_del_all.cnv_report.tsv` | CNV TSV report - pathology_purecn |
| `results/dna/{sample}_{type}/additional_files/cnv/{sample}_{type}.purecn.amp_all_del_validated.cnv_report.tsv` | CNV TSV report validated - PureCN |
| `results/dna/{sample}_{type}/additional_files/cnv/{sample}_{type}.pathology.amp_all_del_validated.cnv_report.tsv` | CNV TSV report validated - pathology |
| `results/dna/{sample}_{type}/cnv/{sample}_{type}.pathology_purecn.amp_all_del_validated.cnv_report.tsv` | CNV TSV report validated - pathology_purecn |
| `results/dna/{sample}_{type}/additional_files/cnv/{sample}_{type}.deletions.tsv` | Small CNV deletions |
| `results/dna/{sample}_{type}/additional_files/cnv/{sample}_{type}.amplifications.tsv` | Small CNV amplifications |
| `results/dna/{sample}_{type}/additional_files/cnv/{sample}_{type}.pathology.jumble.vcf` | Jumble vcf pathology |
| `results/dna/{sample}_{type}/additional_files/cnv/{sample}_{type}.purecn.jumble.vcf` | Jumble vcf purecn |
| `results/dna/{sample}_{type}/additional_files/cnv/{sample}_{type}.pathology_purecn.jumble.vcf` | Jumble vcf pathology_purecn |
| `results/dna/{sample}_{type}/additional_files/cnv/{sample}_{type}.germline.vcf.gz` | Germline vcf used in CNV analysis |
| `results/dna/{sample}_{type}/additional_files/cnv/{sample}_{type}.manta_tumorSV.vcf.gz` | Manta VCF |
| `results/dna/{sample}_{type}/cnv/{sample}_{type}.ctDNA_fraction.tsv` | ctDNA fraction |
| `results/dna/{sample}_{type}/additional_files/cnv/{sample}_{type}.ctDNA_fraction_info.tsv` | ctDNA fraction_info |
| `results/dna/{sample}_{type}/cnv/{sample}_{type}.ITD.hard_filter.vcf` | ITD filtered |

### RNA

#### Alignment
| **File** | **Description** |
|---|---|
| `bam_rna/{sample}_{type}.star_fusion.bam` | STAR Fusion BAM |
| `bam_rna/{sample}_{type}.star_fusion.bam.bai` | STAR Fusion BAM index |

#### SNV and INDELs
| **File** | **Description** |
|---|---|
| `results/rna/{sample}_{type}/id_snps/{sample}_{type}.id_snps.vcf` | ID-SNP VCF RNA |

#### QC
| **File** | **Description** |
|---|---|
| `results/rna/{sample}_{type}/qc/{sample}_{type}.house_keeping_gene_coverage.tsv` | Housekeeping gene coverage |
| `results/rna/qc/multiqc_RNA.html` | MultiQC RNA HTML |
| `results/sample_mixup_check.tsv` | RNA DNA sample mixup report |

#### Fusions
| **File** | **Description** |
|---|---|
| `results/rna/{sample}_{type}/fusion/{sample}_{type}.fusion_report.tsv` | Fusion report TSV |
| `results/rna/{sample}_{type}/additional_files/fusion/{sample}_{type}.arriba.fusions.tsv` | Arriba fusion TSV |
| `results/rna/{sample}_{type}/fusion/{sample}_{type}.arriba.fusions.pdf` | Arriba fusion PDF |
| `results/rna/{sample}_{type}/fusion/{sample}_{type}.exon_skipping.tsv` | Exon skipping TSV |
| `results/rna/{sample}_{type}/additional_files/fusion/{sample}_{type}.ctat.exon_skipping.unfiltered.tsv` | CTAT-splicing unfiltered TSV |
| `results/rna/{sample}_{type}/fusion/{sample}_{type}.ctat.exon_skipping.filtered.tsv` | CTAT-splicing filtered TSV |
| `results/rna/{sample}_{type}/fusion/{sample}_{type}.ctat.exon_skipping.igv.html` | CTAT-splicing html |
| `results/rna/{sample}_{type}/additional_files/fusion/{sample}_{type}.star-fusion.fusion_predictions.tsv` | STAR Fusion prediction TSV |
| `results/rna/{sample}_{type}/additional_files/fusion/{sample}_{type}.fusioncatcher.fusion_predictions.txt` | Fusioncatcher prediction TSV |

## ctDNA Pipeline Output
Files defined in `output_files_ctDNA.yaml` that differ or are specific to ctDNA.

### Alignment
| **File** | **Description** |
|---|---|
| `bam_dna/{sample}_{type}.umi.bam` | BAM umi |
| `bam_dna/{sample}_{type}.umi.bam.bai` | BAM umi index |

### SNV and INDELs
| **File** | **Description** |
|---|---|
| `results/dna/{sample}_{type}/additional_files/vcf/{sample}_{type}.annotated.exon_only.filter.soft_filter.umi.vcf` | Exon-only soft-filtered VCF UMI |
| `results/dna/{sample}_{type}/vcf/{sample}_{type}.annotated.exon_only.filter.hard_filter.codon_snv.umi.qci.vcf` | QCI umi VCF |

### Reports
| **File** | **Description** |
|---|---|
| `results/dna/{sample}_{type}/{sample}_{type}.umi.general_report.html` | Results report umi DNA HTML |

### Biomarker
| **File** | **Description** |
|---|---|
| `results/dna/{sample}_{type}/biomarker/{sample}_{type}.umi.TMB.txt` | TMB umi |
| `results/dna/{sample}_{type}/hla/{sample}_{type}_hla_type_result.tsv` | HLA optitype |
| `results/dna/{sample}_{type}/hla/{sample}_{type}_hla_type_coverage_plot.pdf` | HLA coverage optitype |

### Fusions
| **File** | **Description** |
|---|---|
| `results/dna/{sample}_{type}/fusion/{sample}_{type}.gene_fuse_report.umi.tsv` | GeneFuse report |
| `results/dna/{sample}_{type}/additional_files/fusion/{sample}_{type}.gene_fuse.unfiltered.results.txt` | GeneFuse fusions unfiltered |

### CNV
| **File** | **Description** |
|---|---|
| `results/dna/{sample}_{type}/additional_files/cnv/{sample}_{type}.ITD.vcf` | ITD |

### Fragmentomics
| **File** | **Description** |
|---|---|
| `results/dna/{sample}_{type}/fragmentomics/{sample}_{type}.end-motifs.tsv` | Fragmentomics end-motifs |
| `results/dna/{sample}_{type}/fragmentomics/{sample}_{type}.mds.txt` | Fragmentomics mds |
| `results/dna/{sample}_{type}/fragmentomics/{sample}_{type}.interval-end-motifs.tsv` | Fragmentomics interval-end-motifs |
| `results/dna/{sample}_{type}/fragmentomics/{sample}_{type}.interval-mds.txt` | Fragmentomics interval-mds |
| `results/dna/{sample}_{type}/fragmentomics/{sample}_{type}.frag-length-bins.tsv` | Fragmentomics fragment length bins |
| `results/dna/{sample}_{type}/fragmentomics/{sample}_{type}.fragment_length_score.txt` | Fragmentomics fragment length score |

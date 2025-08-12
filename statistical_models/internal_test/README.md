# NGS Variant Calling Pipeline

A comprehensive Next-Generation Sequencing (NGS) pipeline for processing SRA data through quality control, alignment, variant calling, and annotation with gene panel filtering.

## Overview

This pipeline processes SRA (Sequence Read Archive) files through a complete variant calling workflow, from raw sequencing data to filtered, annotated variants. It implements best practices for NGS data processing including GATK's recommended workflows for variant discovery.

## Pipeline Workflow

### 1. Quality Control (QC)
- **fasterq-dump**: Convert SRA files to FASTQ format
- **fastp**: Adapter trimming and quality filtering
- **BBDuk**: Additional adapter removal and quality trimming
- **IlluQC**: Quality control assessment
- **FastQC**: Generate quality reports
- **MultiQC**: Aggregate QC results

### 2. Alignment
- **BWA MEM**: Align reads to reference genome (human_g1k_v37)
- **Samtools**: SAM/BAM file processing

### 3. Post-Alignment Processing
- **Picard SortSam**: Sort BAM files by coordinate
- **Coverage Analysis**: Calculate depth of coverage statistics
- **Picard MarkDuplicates**: Remove PCR duplicates
- **Picard AddOrReplaceReadGroups**: Add read group information
- **Picard BuildBamIndex**: Index BAM files

### 4. Base Quality Score Recalibration (BQSR)
- **GATK BaseRecalibrator**: Generate recalibration table
- **GATK ApplyBQSR**: Apply base quality recalibration

### 5. Variant Calling
- **GATK HaplotypeCaller**: Call SNPs and indels
- **GATK IndexFeatureFile**: Index VCF files

### 6. Variant Quality Score Recalibration (VQSR)
- **GATK VariantRecalibrator**: Build recalibration model for SNPs and indels
- **GATK ApplyVQSR**: Apply quality score recalibration

### 7. Annotation and Filtering
- **ANNOVAR**: Comprehensive variant annotation
- **Custom filtering**: Gene panel-based filtering with quality thresholds
- **PASS filtering**: Retain only variants that pass quality filters

## Prerequisites

### Required Software
- **SRA Toolkit** (fasterq-dump)
- **fastp** (version with adapter trimming support)
- **BBTools** (BBDuk)
- **IlluQC**
- **FastQC**
- **MultiQC**
- **BWA** (Burrows-Wheeler Aligner)
- **Samtools**
- **Picard Tools**
- **GATK** (Genome Analysis Toolkit)
- **ANNOVAR**
- **GNU Parallel**

### Required Reference Files
- Human reference genome (human_g1k_v37.fasta)
- GATK resource bundle files:
  - dbsnp_138.b37.vcf.gz
  - Mills_and_1000G_gold_standard.indels.b37.vcf.gz
  - hapmap_3.3.b37.vcf.gz
  - 1000G_omni2.5.b37.vcf.gz
  - 1000G_phase1.snps.high_confidence.b37.vcf.gz
- BBTools adapter sequences (adapters.fa)
- ANNOVAR database files
- Gene panel file (list of genes of interest)

## Installation and Setup

1. Install all required software packages
2. Download and index reference genome and GATK resource files
3. Set up ANNOVAR database
4. Update directory paths in the script to match your system
5. Prepare gene panel file with one gene symbol per line
6. Create SRR accession list file

## Usage

```bash
# Make the script executable
chmod +x ngs_pipeline.sh

# Run the pipeline
./ngs_pipeline.sh
```

The script reads SRR accession numbers from a file (`/path/to/SRR_Acc_List.txt`) and processes each sample in parallel.

## Configuration

### Directory Structure
Update the following paths in the script:
- `BASE_DIR`: Working directory for processing
- `RESULTS_DIR`: Final output directory
- `SRA_DATA_DIR`: Location of SRA files
- Reference genome and resource file paths
- ANNOVAR database path
- Gene panel file path

### Processing Parameters
- **Threads**: Set to 12 for most tools (adjust based on available cores)
- **Quality thresholds**: 
  - fastp: mean quality ≥20, minimum length 36bp
  - BBDuk: quality trimming at Q10
  - IlluQC: length ≥70bp, quality ≥20
- **Variant filtering**:
  - Heterozygous variants: ref count ≥8, alt count ≥8, alt frequency ≥25%
  - Homozygous variants: alt count ≥8
  - Excludes Y chromosome variants and benign variants
  - PASS filter: Only variants with "PASS" in FILTER column retained

## Output Files

### Quality Control
- FastQC reports (HTML and zip files)
- MultiQC summary report
- Coverage statistics (average_depth.txt)

### Variant Files
- Raw variants (var.vcf.gz)
- Quality score recalibrated variants (recal_snp_recal_indel.vcf.gz)
- Annotated variants (ANNOVAR output files)
- Filtered variants:
  - `matched_pgl.vcf`: Gene panel matched variants
  - `AD_GT_matched_pgl.vcf`: Quality filtered variants
  - `filtered_without_Y_and_benign_AD_GT_matched_pgl.vcf`: Clinical filtered variants
  - `passed_filtered_without_Y_and_benign_AD_GT_matched_pgl.vcf`: Final filtered variants

### Log Files
- Process logs: `{SRR}_process.log`
- Error logs: `{SRR}_process.err`

## Pipeline Features

- **Parallel processing**: Uses GNU Parallel for efficient sample processing
- **Comprehensive QC**: Multiple quality control steps with detailed reporting
- **Best practices**: Implements GATK best practices for variant calling
- **Gene panel filtering**: Focuses analysis on clinically relevant genes
- **Quality-based filtering**: Applies stringent quality thresholds
- **Automated cleanup**: Removes intermediate files to save disk space
- **Detailed logging**: Comprehensive logging for troubleshooting

## Performance Considerations

- **CPU**: Pipeline uses 12 threads by default (adjust based on available cores)
- **Memory**: GATK steps may require substantial RAM (≥8GB recommended)
- **Storage**: Significant disk space required for intermediate files
- **Processing time**: Varies by sample size (typically several hours per sample)

## Troubleshooting

- Check log files for detailed error messages
- Ensure all required software is installed and in PATH
- Verify reference file paths and formats
- Confirm adequate disk space for processing
- Monitor system resources during execution

## Citation

If you use this pipeline in your research, please cite the appropriate tools and resources:
- GATK
- BWA
- Picard Tools
- ANNOVAR
- FastQC/MultiQC
- And other tools used in the pipeline

## License

This pipeline is provided as-is for research purposes. Please ensure compliance with individual tool licenses.

## Support

For questions or issues with this pipeline, please check the documentation of individual tools or create an issue in this repository.

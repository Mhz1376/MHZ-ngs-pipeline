# MHZ-ngs-pipeline

A comprehensive NGS processing and variant analysis pipeline developed for the research published in *"Exploratory analysis of germline variants in immunodeficiency-associated genes identifies TYR and FANCD2 as candidate loci in early-onset breast cancer"*.

## Overview

This repository contains the complete bioinformatics workflow used to analyze whole-exome sequencing data from 62 early-onset breast cancer patients (\<50 years) against a curated panel of 1,672 immunodeficiency-associated genes. The pipeline integrates quality control, alignment, variant calling, annotation, statistical prioritization, and pathway enrichment to identify candidate susceptibility loci contributing to the missing heritability of early-onset breast cancer.

## Key Findings

Our analysis identified two statistically significant candidate variants:
- **TYR** (rs1126809, 11-89017961-G-A): 58.1% carrier frequency, FDR < 0.05
- **FANCD2** (rs750338758, 3-10088407-AG-A): 14.5% carrier frequency, FDR < 0.05

Both variants were classified as VUS leaning toward likely pathogenic and demonstrated robust statistical evidence across multiple models. Pathway analysis revealed complementary mechanisms: FANCD2 enriched in DNA-repair networks (BRCA1, ATM, TP53 pathways) and FANCD2 linked to PKR-mediated innate immune signaling, while TYR mapped to metabolic/pigmentation pathways with connections to IL-12, IL-4, and IL-13 cytokine networks.

## Repository Structure

### `ngs_pipeline/`
Complete NGS data processing workflow from raw reads to annotated variants:

- **Quality Control & Preprocessing**
  - Raw data conversion (fasterq-dump)
  - Adapter trimming and quality filtering (fastp, BBDuk, IlluQC)
  - QC reporting (FastQC, MultiQC)

- **Alignment & Post-Processing**
  - BWA-MEM alignment to GRCh37/hg19
  - SAM/BAM conversion and sorting (Samtools, Picard)
  - PCR duplicate removal
  - Base quality score recalibration (GATK BQSR)

- **Variant Calling & Annotation**
  - Germline variant calling (GATK HaplotypeCaller)
  - Variant quality score recalibration (VQSR)
  - Multi-database annotation (ANNOVAR)
  - Clinical classification (Franklin/Genoox)
  - Pathogenicity scoring (CADD, REVEL, BayesDel, AlphaMissense, SpliceAI)

- **Filtering Strategy**
  - Stringent quality filters (GQ ≥ 30, AD ≥ 8, VAF ≥ 0.25)
  - Panel gene list (PGL) restriction: 1,672 genes
  - Clinical classification: VUS or higher pathogenicity
  - Exclusion: Y chromosome, benign variants

Detailed documentation and parameters are provided in the directory README.

### `statistical_models/`
Comprehensive statistical prioritization framework for variant analysis:

#### `Variants_Models/` (NEW)
**Per-variant statistical testing pipeline** using three complementary approaches:

1. **Fisher's Exact Test**
   - Classical 2×2 contingency table analysis
   - Two-sided p-values with FDR correction
   - Conservative test for small sample sizes

2. **Barnard's Exact Test** (Fixed implementation for large samples)
   - More powerful unconditional exact test
   - Intelligent handling based on sample size:
     - **Exact method**: Total N ≤ 10,000
     - **Downsampling**: Controls downsampled to 5,000 maintaining carrier proportion
     - **Chi-squared approximation**: Fallback for very large datasets
   - Method tracking for transparency

3. **Bayesian Analysis**
   - Beta-Binomial conjugate model
   - 1,000,000 posterior samples
   - Non-informative priors (α=0.5, β=0.5)
   - Outputs:
     - Posterior median OR with 95% credible intervals
     - Posterior probability OR > 1
     - Two-sided Bayesian p-value

**Features:**
- Multi-threaded parallel processing (configurable workers)
- Chunked processing with memory management
- Thread-safe CSV output with progressive writing
- Homozygote-aware control carrier estimation
- Hardy-Weinberg equilibrium-based calculations when homozygote data unavailable
- Multiple testing corrections (Bonferroni, FDR-BH)
- Comprehensive error handling and skipped variant logging

**Input Requirements:**
- CSV with variant data including:
  - Case carrier counts (SRR_Files or File_Count)
  - gnomAD Allele Count (AC)
  - gnomAD Allele Number (AN)
  - gnomAD Homozygotes (optional but recommended)

**Outputs:**
- `variant_fisher.csv`: Fisher's test results with OR and CI
- `variant_barnard.csv`: Barnard's test results with method tracking
- `variant_bayesian.csv`: Bayesian posterior distributions and credible intervals
- `skipped_variants.csv`: Variants excluded with reasons

Complete usage instructions, configuration options, and technical details in `Variants_Models/README.md`.

#### Legacy Models (Not used in final publication)
- Model 1: Internal weighted scoring
- Model 2: Alternative Fisher implementation
- Model 3: Alternative Barnard implementation
- Model 4: Alternative Bayesian implementation

*Note: These models were part of exploratory analyses but were superseded by the `Variants_Models` approach in the final manuscript.*

## Computational Environment

- **Platform**: WSL2 (Windows 11)
- **Languages**: Bash, Python 3.7+, R
- **Key Dependencies**: 
  - BWA-MEM, Samtools, Picard, GATK
  - Python: pandas, numpy, scipy (≥1.9.0), statsmodels, tqdm
  - R: Statistical analysis packages
  
See `version` file in data directory for complete software specifications.

## Data Availability

### Processed Variant Dataset
The complete processed dataset including VCF files and statistical analysis results is publicly available:

**Zenodo DOI:** https://doi.org/10.5281/zenodo.17561209

Dataset includes:
- Filtered VCF files for all 62 cases
- Annotated variant tables
- Project / protein folder contents

### Raw Sequencing Data
Source WES data obtained from NCBI Sequence Read Archive (SRA):
- Reference: Vellichirammal et al. (2023), *Hum Genomics* 17:64
- Study focus: US Midwestern breast cancer cohort
- Access: Via SRA database (NCBI)

**Important:** Please read the README file in the data directory before use.

## Analysis Strategy

### Gene Panel Construction (1,672 genes)
Integrated curation from:
- **PanelApp**: Canonical immunodeficiency genes
- **DisGeNET**: Disease-gene associations
- **GeneCards**: Extended immune-related annotations

*Note: HLA genes excluded due to insufficient WES coverage for reliable variant calling.*

### Statistical Prioritization Workflow
1. Case cohort: 62 early-onset breast cancer patients (age < 50)
2. Control: gnomAD population frequencies
3. 2×2 contingency table construction (Cases carriers/non-carriers vs. Controls carriers/non-carriers)
4. Three complementary statistical tests per variant
5. Multiple testing correction (FDR-BH)
6. Prioritization threshold: FDR < 0.05

### Pathway & Functional Enrichment
- **Databases**: GeneAnalytics, PathCards, MalaCards, GeneCards, Pharos
- **Focus**: DNA repair, immune signaling, metabolic pathways
- **Disease associations**: Breast cancer, inherited cancer syndromes
- **GWAS integration**: NHGRI-EBI GWAS Catalog with LD analysis (1000G Phase 3)

### Structural & Computational Analysis
- **AlphaFold modeling**: Protein structure prediction (TYR p.Arg402Gln)
- **In silico tools**: REVEL, BayesDel, AlphaMissense, SpliceAI, CADD
- **Mutational signatures**: Alignment with COSMIC SBS patterns

## Key Results Summary

### Variant Statistics
| Variant | Gene | rsID | Carriers | Frequency | Fisher FDR | Barnard FDR | Bayesian Pr(OR>1) |
|---------|------|------|----------|-----------|------------|-------------|-------------------|
| 11-89017961-G-A | TYR | rs1126809 | 36/62 | 58.1% | 1.3×10⁻³ | 1.7×10⁻³ | 0.99 |
| 3-10088407-AG-A | FANCD2 | rs750338758 | 9/62 | 14.5% | 2.4×10⁻¹¹ | 1.9×10⁻³ | 1.00 |

### BRCA1/2 Findings
- **BRCA1**: 1 pathogenic frameshift mutation (1/62 cases)
- **BRCA2**: 1 pathogenic missense mutation (1/62 cases)
- Near-absence highlights "missing heritability" gap

### Pathway Enrichment
**FANCD2** (DNA Repair Axis):
- BRCA1 pathway (22 breast cancer diseases)
- ATM signaling (13 breast cancer diseases)
- TP53 regulation (13 breast cancer diseases)
- PKR-mediated innate immune signaling

**TYR** (Metabolic/Immune Axis):
- Pigmentation pathways (melanin/pheomelanin biosynthesis)
- Glucose/energy metabolism
- Links to inherited breast cancer syndromes (Cowden, Peutz-Jeghers)
- Cytokine networks (IL-12, IL-4, IL-13)
- COSMIC SBS38 association (UV/oxidative damage)

## Usage

### Prerequisites
```bash
# Install Python dependencies
pip install pandas numpy scipy statsmodels tqdm

# Verify scipy version (≥1.9.0 required for barnard_exact)
python -c "import scipy; print(scipy.__version__)"
```

### Running the Pipeline

1. **NGS Processing**
```bash
cd ngs_pipeline/
# Follow step-by-step instructions in ngs_pipeline/README.md
```

2. **Statistical Analysis**
```bash
cd statistical_models/Variants_Models/

# Edit CONFIG section in the script to set:
# - Input file path
# - Output file paths
# - Total case count
# - Thread/chunk parameters

python variant_analysis_pipeline.py
```

3. **Configuration Example**
```python
CONFIG = {
    'input_file': 'variant_analysis_full.csv',
    'total_cases': 62,
    'n_workers': 12,
    'chunk_size': 12,
    'significance_level': 0.05,
    # ... see script for full options
}
```

### Expected Runtime
- NGS pipeline: ~2-4 hours per sample (depends on coverage)
- Statistical analysis: ~30-60 minutes for 622 variants (12 threads)

## Reproducibility

All analyses are fully reproducible with:
- Fixed random seeds for Bayesian sampling
- Deterministic sorting and processing order
- Version-controlled software dependencies
- Public data availability
- Complete parameter documentation

## Study Limitations

1. **Sample size**: Limited to 62 cases, reducing statistical power
2. **Control cohort**: gnomAD-derived frequencies introduce ancestry heterogeneity
3. **Targeted panel**: 1,672 genes; non-targeted regions not assessed
4. **VUS classification**: Prioritized variants lack functional validation
5. **Correlative**: Pathway associations are hypothesis-generating, not causal

## Future Directions

- Large-scale replication in ancestry-matched cohorts (n > 1,000)
- Functional validation:
  - Minigene splicing assays (FANCD2)
  - DNA repair functional assays (γH2AX, comet assay)
  - Protein stability and activity assays (TYR)
  - Cytokine profiling in patient-derived cells
- Multi-omics integration (RNA-seq, proteomics)
- Epistasis testing (FANCD2 × TYR interactions)
- Clinical validation for risk stratification

## Citation

If you use this pipeline or data in your research, please cite:

```
Zabihi MH, Kalantar SM, Dehghan HR, Sheikhha MH. (2025). 
Exploratory analysis of germline variants in immunodeficiency-associated 
genes identifies TYR and FANCD2 as candidate loci in early-onset breast cancer. 
[Journal Name]. [DOI]
```

**Primary data source:**
```
Vellichirammal NN, et al. (2023). The mutational landscape of a US Midwestern 
breast cancer cohort reveals subtype-specific cancer drivers and prognostic markers. 
Hum Genomics 17:64. https://doi.org/10.1186/s40246-023-00511-6
```

## Contributing

We welcome contributions, bug reports, and feature requests:
- Open an issue for bugs or questions
- Submit pull requests for improvements
- Contact us for collaboration opportunities

## License

[MIT License]

## Ethics Statement

This study was approved by the Shahid Sadoughi University of Medical Sciences Ethics Committee (IR.SSU.MEDICINE.REC.1403.108). All analyses were conducted on de-identified, publicly available sequencing data.

## Acknowledgments

- Bash pipeline automation and error handling: Research team with ChatGPT-3.5 assistance
- Python statistical models: Research team with ChatGPT-5 for optimization
- R statistical scripts: Research team with ChatGPT-5 for refactoring
- Manuscript preparation: Research team with Google NotebookLM and language polishing
- All AI-assisted outputs were critically reviewed and approved by authors

## Contact

For questions, collaboration requests, or bug reports:

**Mohammad Hosein Zabihi**  
Department of Medical Genetics  
Shahid Sadoughi University of Medical Sciences, Yazd, Iran

📧 **Mohammad.hosein.zabihi@gmail.com**

**Subject Line:** Please include "MHZ-ngs-pipeline" in your email subject for efficient tracking.

**ORCID:** https://orcid.org/0000-0001-9710-3046

---

*This pipeline was developed to investigate the genetic architecture underlying missing heritability in early-onset breast cancer through comprehensive analysis of immunodeficiency-associated germline variants.*
 
**Data DOI:** https://doi.org/10.5281/zenodo.17561209

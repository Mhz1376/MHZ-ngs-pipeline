# Gene-Level Burden Analysis using Fisher's Exact Test

A comprehensive R script for performing gene-level burden analysis on variant data using Fisher's exact test with population frequency-based controls and advanced statistical corrections.

## Overview

This script performs gene-level burden analysis to identify genes with significant enrichment of rare variants in case samples compared to population-based expected frequencies. It implements Fisher's exact test with multiple testing corrections and generates publication-ready visualizations.

## Key Features

- **Parallel processing**: Efficient multi-core data processing
- **Robust data handling**: Safe CSV reading with automatic column detection
- **Population-based controls**: Uses gnomAD and population frequency data for control estimates
- **Statistical rigor**: Fisher's exact test with Bonferroni and FDR corrections
- **Publication-ready plots**: Journal-style volcano plots, bar charts, and distribution plots
- **Comprehensive logging**: Detailed progress tracking and error handling

## Analysis Workflow

### 1. Data Loading and Preprocessing
- Reads variant data from CSV files for multiple samples
- Handles missing columns and standardizes gene names
- Processes allele frequencies (gnomAD_AF, PopFreqMax)
- Filters out invalid entries and standardizes data format

### 2. Sample-Level Burden Calculation
- Identifies genes with at least one qualifying variant per sample
- Calculates binary burden status (presence/absence) per gene per sample
- Aggregates burden counts across all case samples

### 3. Population Burden Estimation
- Uses exact product method for population burden estimation
- Incorporates variant allele frequencies from gnomAD and population databases
- Calculates expected burden frequencies using Hardy-Weinberg principles
- Accounts for multiple variants per gene using probability theory

### 4. Statistical Testing
- Performs Fisher's exact test for each gene
- Compares observed case burden vs. expected population burden
- Calculates odds ratios with confidence intervals
- Handles edge cases and statistical exceptions

### 5. Multiple Testing Correction
- Applies Bonferroni correction for family-wise error rate control
- Implements Benjamini-Hochberg FDR correction
- Identifies significant genes at different significance thresholds

### 6. Visualization and Output
- Generates volcano plots with selective gene labeling
- Creates bar charts of top significant genes
- Produces burden distribution histograms
- Exports results in multiple formats

## Prerequisites

### Required R Packages
```r
# Core packages (automatically installed if missing)
- parallel      # Multi-core processing
- data.table    # Efficient data manipulation
- dplyr         # Data wrangling
- readr         # Fast CSV reading
- ggplot2       # Advanced plotting
- ggrepel       # Smart label positioning
- scales        # Plot scaling utilities
```

### Input Data Requirements

#### File Structure
```
/path/to/results/directory/
├── SRR_ID_1/
│   └── ANNOTATION/
│       └── Franklin_Results_Without_Benign_LikelyBenign_0.001.csv
├── SRR_ID_2/
│   └── ANNOTATION/
│       └── Franklin_Results_Without_Benign_LikelyBenign_0.001.csv
└── ...
```

#### Required CSV Columns
- **Gene Name**: Gene symbol or identifier
- **gnomAD_AF**: gnomAD allele frequency (numeric)
- **PopFreqMax**: Maximum population frequency (numeric)
- **Variant**: Unique variant identifier (auto-generated if missing)

#### Sample List File
- Text file with one SRR identifier per line
- Located at `/path/to/srr_names.txt`

## Installation and Setup

1. **Install R** (version ≥ 4.0 recommended)

2. **Clone repository and set up paths**:
```r
# Update these paths in the script:
base_path <- "/path/to/results/directory"           # NGS results location
srr_names_file <- "/path/to/srr_names.txt"          # Sample list file
results_dir <- "/path/to/output/directory"          # Output location
```

3. **Prepare input data**:
   - Ensure CSV files are properly formatted
   - Create SRR names file with sample identifiers
   - Verify directory structure matches expected format

## Usage

### Basic Execution
```r
# Run the complete analysis
source("gene_burden_analysis.R")
```

### Advanced Usage
```r
# Run individual components
results <- main_burden_analysis()
create_summary_plots(results)
```

### Parallel Processing Control
The script automatically detects available cores and uses up to 12 cores. To modify:
```r
num_cores <- min(8, available_cores - 1)  # Use 8 cores maximum
```

## Output Files

### Statistical Results
- **Gene_Burden_Analysis_Results.csv**: Complete analysis results
- **Significant_Genes_FDR_0.05.csv**: FDR-significant genes
- **Significant_Genes_Raw_P_0.05.csv**: Nominally significant genes

### Visualizations
- **Gene_Burden_Volcano_Plot_annotated.png**: Volcano plot with selective labeling
- **Top_Significant_Genes_Styled.png**: Bar chart of top genes
- **Burden_Distribution_Styled.png**: Distribution histogram

### Result Columns
- **Gene**: Gene symbol
- **Cases_Burden**: Number of cases with burden variants
- **Controls_Burden**: Expected number from population
- **P_Value**: Raw Fisher's exact test p-value
- **P_Value_FDR**: FDR-corrected p-value
- **P_Value_Bonferroni**: Bonferroni-corrected p-value
- **Odds_Ratio**: Effect size estimate
- **CI_Lower/CI_Upper**: 95% confidence intervals

## Statistical Methods

### Fisher's Exact Test
- Tests association between gene burden and case status
- Uses exact probabilities for small sample sizes
- Provides unbiased p-values and confidence intervals

### Population Control Estimation
- Calculates gene-level burden probability using variant frequencies
- Uses product formula: P(no burden) = ∏(1 - 2×AF_i)
- Accounts for multiple variants per gene
- Applies Hardy-Weinberg equilibrium assumptions

### Multiple Testing Correction
- **Bonferroni**: Controls family-wise error rate (FWER)
- **FDR**: Controls false discovery rate using Benjamini-Hochberg procedure
- Recommended threshold: FDR < 0.05

## Visualization Features

### Volcano Plot
- **X-axis**: log2(Odds Ratio)
- **Y-axis**: -log10(FDR p-value)
- **Point size**: Number of cases with burden
- **Colors**: Significance status (FDR < 0.05)
- **Labels**: Selected high-priority genes (HLA-DQA1, HLA-DRB1, HLA-DQB1, KMT2C, HLA-DRB5, PROM1)

### Bar Chart
- Shows top 25 most significant genes
- Ranked by -log10(raw p-value)
- Journal-style formatting with consistent theme

### Distribution Plot
- Histogram of burden counts across genes
- Shows distribution of variant burden in the dataset

## Performance Considerations

### System Requirements
- **RAM**: ≥8GB recommended for large datasets
- **CPU**: Multi-core processor for parallel processing
- **Storage**: Sufficient space for intermediate files and outputs

### Optimization Tips
- Use parallel processing for large sample sets
- Filter input data to reduce memory usage
- Monitor system resources during execution

## Troubleshooting

### Common Issues

1. **Missing files**: Check file paths and directory structure
2. **Column name mismatches**: Verify CSV column headers
3. **Memory issues**: Reduce parallel cores or filter input data
4. **No significant results**: Check input data quality and sample size

### Error Handling
- Automatic package installation for missing dependencies
- Graceful handling of missing or corrupted files
- Comprehensive error logging with detailed messages
- Fallback options for failed statistical tests

## Quality Control

### Input Validation
- Checks for required columns and data types
- Handles missing values and invalid entries
- Validates file existence before processing

### Statistical Validity
- Excludes genes with insufficient data
- Handles edge cases in Fisher's test
- Provides confidence intervals for effect sizes

## Citation

If you use this analysis pipeline in your research, please cite:
- Fisher, R.A. (1922). On the interpretation of χ² from contingency tables, and the calculation of P. *Journal of the Royal Statistical Society*, 85(1), 87-94.
- Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society*, 57(1), 289-300.

## License

This script is provided for research purposes. Please ensure compliance with individual package licenses.

## Support

For questions or issues:
1. Check the console output for detailed error messages
2. Verify input data format and file paths
3. Ensure all required packages are installed
4. Review the troubleshooting section above


# Gene-Level Burden Analysis using Fisher's Exact Test

A comprehensive R pipeline for performing gene-level burden analysis using Fisher's exact test with high-quality publication-ready visualizations.

## Overview

This pipeline analyzes genetic variant burden at the gene level by:
- Loading variant data from multiple samples
- Calculating gene-level burden for each sample
- Estimating population-level burden frequencies
- Performing Fisher's exact tests for case-control comparison
- Applying multiple testing corrections (Bonferroni and FDR)
- Generating publication-quality plots

## Features

- **Parallel Processing**: Automatically detects available cores and parallelizes data loading
- **Robust Error Handling**: Comprehensive error checking and graceful failure recovery
- **High-Quality Graphics**: Publication-ready plots with vector (PDF/EPS) and high-resolution raster (TIFF@1200 DPI) output
- **Font Management**: Automatic system font detection and registration
- **Multiple Testing Correction**: Both Bonferroni and FDR corrections applied
- **Flexible Input**: Handles missing files and inconsistent data formats

## Requirements

### R Version
- R >= 4.0.0

### Required R Packages
```r
# Core packages (automatically installed if missing)
required_packages <- c(
  "parallel", "data.table", "dplyr", "readr",
  "ggplot2", "ggrepel", "scales",
  "showtext", "sysfonts"
)

# Optional but recommended for better graphics
install.packages("ragg")  # Improved TIFF rendering
```

### System Requirements
- Multi-core CPU recommended for parallel processing
- Sufficient RAM for loading all variant files simultaneously
- System fonts (Arial preferred, falls back to system defaults)

## Input Data Format

### Directory Structure
```
/path/to/results/directory/
├── SRR_SAMPLE_1/
│   └── ANNOTATION/
│       └── Franklin_Results_Without_Benign_LikelyBenign_0.001.csv
├── SRR_SAMPLE_2/
│   └── ANNOTATION/
│       └── Franklin_Results_Without_Benign_LikelyBenign_0.001.csv
└── ...
```

### Required Files
1. **SRR Names File**: Text file containing sample identifiers (one per line)
   ```
   SRR_SAMPLE_1
   SRR_SAMPLE_2
   SRR_SAMPLE_3
   ```

2. **Variant Files**: CSV files with the following required columns:
   - `Gene Name` (or similar): Gene symbol/identifier
   - `Variant`: Variant identifier (created automatically if missing)
   - `gnomAD_AF`: gnomAD allele frequency (optional, defaults to 0)
   - `PopFreqMax`: Maximum population frequency (optional, defaults to 0)

## Configuration

Update the following paths in the script before running:

```r
# Main analysis paths
base_path <- "/path/to/results/directory"           # Directory containing sample folders
srr_names_file <- "/path/to/srr_names.txt"          # File listing sample names
results_dir <- "/path/to/output/directory"          # Output directory for results
```

## Usage

### Basic Usage
```r
# Source the script
source("gene_burden_analysis.R")

# The script runs automatically when sourced
# Results will be saved to the specified output directory
```

### Function-Level Usage
```r
# Load the script without auto-execution
# (comment out the main execution block first)

# Run individual components
results <- main_burden_analysis()
create_summary_plots(results)
```

## Output Files

### Results Tables
- `Gene_Burden_Analysis_Results.csv`: Complete results for all genes
- `Significant_Genes_FDR_0.05.csv`: Genes significant after FDR correction
- `Significant_Genes_Raw_P_0.05.csv`: Genes with raw p-value < 0.05

### Visualizations
All plots saved in multiple formats (PDF, EPS, TIFF@1200 DPI):

1. **Fig1_Gene_Burden_Volcano**: Volcano plot showing log2(OR) vs -log10(FDR p-value)
2. **Fig2_Top_Significant_Genes**: Bar plot of top significant genes
3. **Fig3_Burden_Distribution**: Histogram of burden distribution across genes

### Results Columns
- `Gene`: Gene symbol
- `Cases_Burden`: Number of cases with burden variants
- `Cases_No_Burden`: Number of cases without burden variants  
- `Controls_Burden`: Estimated number of controls with burden variants
- `Controls_No_Burden`: Estimated number of controls without burden variants
- `P_Value`: Raw Fisher's exact test p-value
- `P_Value_Bonferroni`: Bonferroni-corrected p-value
- `P_Value_FDR`: FDR-corrected p-value
- `Odds_Ratio`: Estimated odds ratio
- `CI_Lower`, `CI_Upper`: 95% confidence interval bounds
- `Significant_Bonferroni`: Boolean indicating Bonferroni significance
- `Significant_FDR`: Boolean indicating FDR significance

## Methodology

### Burden Calculation
1. **Sample Burden**: Genes with ≥1 qualifying variant per sample
2. **Population Burden**: Estimated using exact product method based on variant allele frequencies
3. **Statistical Test**: Fisher's exact test comparing case vs. estimated control burden

### Population Burden Estimation
```
For each gene with variants [v1, v2, ..., vn]:
- Carrier probability per variant: min(2 × AF, 1)
- Gene burden probability: 1 - ∏(1 - carrier_prob_i)
- Expected burden count(the min amount is one): sample_size × gene_burden_probability
```

### Multiple Testing Correction
- **Bonferroni**: Controls family-wise error rate (FWER)
- **FDR**: Controls false discovery rate using Benjamini-Hochberg method

## Customization

### Parallel Processing
```r
# Modify number of cores (default: auto-detect - 1, max 12)
num_cores <- min(8, available_cores - 1)
```

### Plot Styling
```r
# Modify plot parameters
base_size <- 10          # Base font size
width_cm <- 17.6         # Plot width in cm
height_cm <- 12.0        # Plot height in cm
tiff_res <- 1200         # TIFF resolution (DPI)
```

### Gene Labeling
```r
# Modify genes to highlight in volcano plot
genes_to_label <- c("HLA-DQA1", "HLA-DRB1", "HLA-DQB1", "KMT2C", "HLA-DRB5", "PROM1")
```

## Troubleshooting

### Common Issues

1. **Memory Issues**: Reduce number of parallel cores or process samples in batches
2. **Font Problems**: Script will fall back to system defaults if Arial not found
3. **Missing Files**: Script continues with available files and reports missing ones
4. **Graphics Errors**: Install `ragg` package for better TIFF rendering

### Error Messages
- `"SRR names file not found"`: Check path to sample names file
- `"No sample files were successfully loaded"`: Verify file paths and formats
- `"No valid p-values were calculated"`: Check input data quality

## Performance Notes

- **Runtime**: Typically 5-15 minutes for 50-100 samples with 10,000-20,000 variants
- **Memory**: ~2-8 GB RAM depending on dataset size
- **Cores**: Scales well with additional CPU cores for data loading phase

## Citation

If you use this pipeline in your research, please cite appropriately and consider sharing the methodology details in your methods section.

## License

This pipeline is provided as-is for research purposes. Please ensure compliance with individual tool licenses.

## Support

For questions or issues with this pipeline, please check the documentation of individual tools or create an issue in this repository.

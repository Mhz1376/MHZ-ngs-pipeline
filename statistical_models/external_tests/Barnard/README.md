# Gene-Level Burden Analysis using Barnard's Test

A comprehensive R pipeline for performing gene-level burden analysis using Barnard's exact test with Fisher's test fallback, featuring high-quality publication-ready visualizations.

## Overview

This pipeline analyzes genetic variant burden at the gene level by:
- Loading variant data from multiple samples
- Calculating gene-level burden for each sample  
- Estimating population-level burden frequencies using exact product method
- Performing Barnard's exact tests for case-control comparison (with Fisher's test fallback)
- Applying multiple testing corrections (Bonferroni and FDR)
- Generating publication-quality plots with vector and high-resolution raster formats

## Key Features

- **Advanced Statistics**: Uses Barnard's exact test (more powerful than Fisher's test) with automatic Fisher's test fallback
- **Parallel Processing**: Automatically detects available cores and parallelizes computations
- **Robust Error Handling**: Comprehensive error checking and graceful failure recovery
- **High-Quality Graphics**: Publication-ready plots in PDF/EPS vector formats plus TIFF@1200 DPI raster
- **Font Management**: Automatic system font detection with showtext integration
- **Haldane Continuity Correction**: Configurable continuity correction for odds ratio estimation
- **Flexible Input**: Handles missing files and inconsistent data formats

## Requirements

### R Version
- R >= 4.0.0

### Required R Packages
```r
# Core packages (automatically installed if missing)
required_packages <- c(
  "parallel", "data.table", "dplyr", "readr",
  "ggplot2", "ggrepel", "scales", "DescTools",
  "showtext", "sysfonts"
)

# Optional but recommended for better graphics
install.packages("ragg")  # Improved TIFF rendering
```

### System Requirements
- Multi-core CPU recommended for parallel processing
- Sufficient RAM for loading all variant files simultaneously (typically 2-8 GB)
- System fonts (Arial preferred, falls back to system defaults)

## Configuration

### Key Parameters
Update these configuration variables at the top of the script:

```r
# Primary configuration
RESULTS_DIR <- "/path/to/output/directory"           # Output directory
BASE_PATH   <- "/path/to/results/directory"          # Sample data directory  
SRR_FILE    <- "/path/to/srr_names.txt"             # Sample names file
VERBOSE     <- TRUE                                  # Enable detailed logging

# Statistical parameters
HALDANE_PC  <- 0.6    # Haldane continuity correction for OR
LOG2OR_MIN  <- 0      # Minimum log2(OR) for volcano plot
LOG2OR_MAX  <- 15     # Maximum log2(OR) for volcano plot (clipped axis)
TOP_LABEL_N <- 17     # API compatibility parameter
```

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
1. **SRR Names File**: Plain text file with sample identifiers (one per line)
   ```
   SRR_SAMPLE_1
   SRR_SAMPLE_2
   SRR_SAMPLE_3
   ```

2. **Variant Files**: CSV files with columns:
   - `Gene Name` (or similar): Gene symbol/identifier (required)
   - `Variant`: Variant identifier (auto-generated if missing)
   - `gnomAD_AF`: gnomAD allele frequency (optional, defaults to 0)
   - `PopFreqMax`: Maximum population frequency (optional, defaults to 0)

## Usage

### Standard Usage
```r
# Simply source the script - it runs automatically
source("gene_burden_analysis_barnard.R")

# Results and plots will be saved to RESULTS_DIR
```

### Function-Level Usage
```r
# For custom workflows, comment out the main execution block and use:
results <- main_burden_analysis()
create_summary_plots(results, top_n = 17)
```

## Output Files

### Results Tables
- `Gene_Burden_Analysis_Results.csv`: Complete results for all analyzed genes
- `Significant_Genes_FDR_0.05.csv`: Genes significant after FDR correction (q < 0.05)
- `Significant_Genes_Raw_P_0.05.csv`: Genes with nominal significance (p < 0.05)

### Visualizations
All plots saved in multiple high-quality formats (PDF, EPS, TIFF@1200 DPI):

1. **Fig1_Gene_Burden_Volcano**: 
   - Volcano plot: log2(OR) vs -log10(FDR p-value)
   - X-axis clipped to [0,15] for better visualization
   - Highlights specific genes: HLA-DQA1, HLA-DRB1, HLA-DQB1, KMT2C, HLA-DRB5, PROM1

2. **Fig2_Top_Significant_Genes**: 
   - Horizontal bar plot of top significant genes
   - Ordered by statistical significance

3. **Fig3_Burden_Distribution**: 
   - Histogram showing distribution of burden counts across genes

## Results Columns

| Column | Description |
|--------|-------------|
| `Gene` | Gene symbol/identifier |
| `Cases_Burden` | Number of cases with burden variants |
| `Cases_No_Burden` | Number of cases without burden variants |
| `Controls_Burden` | Estimated number of controls with burden (integer) |
| `Controls_No_Burden` | Estimated number of controls without burden |
| `Controls_Burden_Fraction` | Exact fractional control burden estimate |
| `P_Value` | Raw Barnard's test p-value (Fisher's fallback if needed) |
| `P_Value_Bonferroni` | Bonferroni-corrected p-value |
| `P_Value_FDR` | FDR-corrected p-value (Benjamini-Hochberg) |
| `Odds_Ratio` | Estimated odds ratio (with Haldane correction) |
| `CI_Lower`, `CI_Upper` | 95% confidence interval bounds |
| `Significant_Bonferroni` | Boolean: significant after Bonferroni correction |
| `Significant_FDR` | Boolean: significant after FDR correction |

## Statistical Methodology

### Burden Calculation
1. **Sample Burden**: Genes with ≥1 qualifying variant per sample
2. **Population Burden**: Estimated using exact product method:
   ```
   For gene with variants [v1, v2, ..., vn]:
   - Carrier probability per variant: min(2 × AF, 1)
   - Gene burden probability: 1 - ∏(1 - carrier_prob_i)
   - Expected burden count: sample_size × gene_burden_probability
   ```

### Statistical Testing
- **Primary**: Barnard's exact test (z-pooled method via DescTools::BarnardTest)
- **Fallback**: Fisher's exact test with Haldane continuity correction
- **Odds Ratio**: Calculated with configurable Haldane correction (default: 0.6)

### Multiple Testing Correction
- **Bonferroni**: Controls family-wise error rate (FWER)
- **FDR**: Controls false discovery rate using Benjamini-Hochberg procedure

## Customization

### Statistical Parameters
```r
HALDANE_PC  <- 0.6    # Continuity correction (0.5 = standard, 0.6 = default)
LOG2OR_MAX  <- 15     # Maximum OR for volcano plot (higher values clipped)
```

### Parallel Processing
```r
# Modify core usage (default: auto-detect - 1, max 12)
num_cores <- min(8, available_cores - 1)
```

### Gene Labeling
```r
# Modify genes highlighted in volcano plot
genes_to_label <- c("HLA-DQA1", "HLA-DRB1", "HLA-DQB1", "KMT2C", "HLA-DRB5", "PROM1")
```

### Plot Styling
```r
base_size <- 10          # Base font size
width_cm <- 17.6         # Plot width (cm)
height_cm <- 12.0        # Plot height (cm)  
tiff_res <- 1200         # TIFF resolution (DPI)
```

## Advantages of Barnard's Test

- **More Powerful**: Often more sensitive than Fisher's exact test
- **Exact**: Provides exact p-values for 2×2 contingency tables
- **Unconditional**: Doesn't condition on marginal totals like Fisher's test
- **Robust Fallback**: Automatically uses Fisher's test if Barnard's fails

## Performance Notes

- **Runtime**: 10-30 minutes for 50-100 samples with 10,000-20,000 variants
- **Memory**: 2-8 GB RAM depending on dataset size
- **CPU Scaling**: Benefits significantly from multi-core processing
- **Barnard's Test**: Computationally intensive but more powerful than Fisher's

## Troubleshooting

### Common Issues
1. **Memory Problems**: Reduce `num_cores` or process samples in batches
2. **DescTools Installation**: Ensure DescTools package installed for Barnard's test
3. **Long Runtime**: Barnard's test is slower than Fisher's - this is expected
4. **Graphics Issues**: Install `ragg` package for optimal TIFF rendering

### Error Messages
- `"SRR names file not found"`: Check SRR_FILE path
- `"No files to read"`: Verify BASE_PATH and file structure
- `"No valid p-values"`: Check input data quality and variant counts

### Performance Optimization
- Use solid-state storage for input files
- Ensure adequate RAM (8GB+ recommended for large datasets)
- Consider reducing dataset size for initial testing

## Quality Control

The pipeline includes extensive quality checks:
- File existence verification
- Data format validation
- Missing value handling
- Statistical test convergence monitoring
- Automatic fallback mechanisms

## Citation

When using this pipeline, please cite:
- Barnard's exact test implementation (DescTools package)
- Population genetics methodology used
- Any relevant statistical methods papers

Include methodology details in your manuscript's methods section, particularly the use of Barnard's test and population burden estimation approach.

## License

This pipeline is provided as-is for research purposes. Please ensure compliance with individual tool licenses.

## Support

For questions or issues with this pipeline, please check the documentation of individual tools or create an issue in this repository.


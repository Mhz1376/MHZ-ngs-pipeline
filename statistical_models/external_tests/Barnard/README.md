# Gene-Level Burden Analysis using Barnard's Test

A high-performance R pipeline for conducting gene-level variant burden analysis in case-control genomic studies using Barnard's exact test. This tool provides more powerful statistical testing than Fisher's exact test, especially for unbalanced designs, while estimating control burden from population frequency data.

## Features

- **Barnard's Exact Test**: Uses the unconditional exact test (`DescTools::BarnardTest`) for superior statistical power compared to Fisher's exact test
- **Population-Based Controls**: Estimates control burden from gnomAD/PopFreqMax frequencies using Hardy-Weinberg equilibrium
- **Fractional Control Counts**: Maintains fractional expected control counts for accurate population modeling
- **Robust Statistical Framework**: Automatic fallback to Fisher's test with Haldane continuity correction when Barnard's test fails
- **High-Performance Computing**: Parallel processing with automatic core detection and load balancing
- **Comprehensive Visualization**: Generates volcano plots, significance rankings, and burden distributions
- **Flexible Gene Annotation**: Targeted labeling of specific genes of interest in plots

## Requirements

### R Version
- R ≥ 4.0.0

### Required R Packages
```r
# Core packages (auto-installed if missing)
- parallel      # Multi-core processing
- data.table    # High-performance data manipulation
- dplyr         # Data wrangling
- readr         # Fast CSV I/O
- ggplot2       # Advanced plotting
- ggrepel       # Smart plot label positioning
- scales        # Plot scaling utilities
- DescTools     # Barnard's exact test implementation
```

## Statistical Methods

### Barnard's Test vs Fisher's Test
- **Barnard's Test**: Unconditional exact test that considers all possible tables with the same marginal totals, providing higher statistical power
- **Fisher's Test**: Conditional exact test used as fallback with Haldane continuity correction (α = 0.5)
- **Odds Ratio Estimation**: Uses Haldane-corrected odds ratios for stable estimates with small counts

### Population Control Modeling
For genes with variants having allele frequencies AF₁, AF₂, ..., AFₙ:
```
Carrier probability per variant = min(2 × AFᵢ, 1)
Gene burden probability = 1 - ∏(1 - carrier_probᵢ)
Expected controls with burden = sample_size × burden_probability
```

## Input Data Format

### Directory Structure
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

### Required CSV Columns
Each variant annotation file must contain:
- **Gene Name** (or gene-like column that will be auto-detected)
- **Variant** (auto-generated if missing: "var_1", "var_2", ...)
- **gnomAD_AF** (gnomAD allele frequency, defaults to 0 if missing)
- **PopFreqMax** (maximum population frequency, defaults to 0 if missing)

### Sample List File
Create a plain text file with one sample ID per line:
```
SRR_ID_1
SRR_ID_2
SRR_ID_3
...
```

## Configuration

### Key Parameters
Edit these variables in the script header:

```r
# File paths (REQUIRED - update before running)
RESULTS_DIR <- "/path/to/output/directory"     # Output directory
BASE_PATH   <- "/path/to/results/directory"    # Input data directory  
SRR_FILE    <- "/path/to/srr_names.txt"       # Sample list file

# Statistical parameters
HALDANE_PC  <- 0.5    # Haldane continuity correction constant
LOG2OR_MIN  <- 0      # Minimum log2(OR) for volcano plot x-axis
LOG2OR_MAX  <- 15     # Maximum log2(OR) for volcano plot x-axis

# Visualization parameters
TOP_LABEL_N <- 17     # Number of top FDR-significant genes to consider for labeling
VERBOSE     <- TRUE   # Enable detailed logging
```

### Gene Labeling in Plots
The volcano plot specifically labels these genes of interest (if present in results):
```r
genes_to_label <- c("HLA-DQA1", "HLA-DRB1", "HLA-DQB1", "KMT2C", "HLA-DRB5", "PROM1")
```
Modify this list to highlight your genes of interest.

## Usage

### Basic Execution
```bash
Rscript burden_analysis_barnard.R
```

### Parallel Processing
- Automatically detects CPU cores and uses up to 12 cores
- Falls back gracefully to sequential processing if parallel setup fails
- Exports all necessary functions and data to worker processes

## Statistical Workflow

1. **Data Loading**: Parallel reading and validation of variant annotation files
2. **Burden Calculation**: Per-sample gene-level burden (≥1 qualifying variant per gene)
3. **Population Modeling**: Hardy-Weinberg-based estimation of expected control burden
4. **Statistical Testing**: 
   - Primary: Barnard's exact test (z-pooled method)
   - Fallback: Fisher's exact test with Haldane correction
5. **Multiple Testing**: Bonferroni and FDR (Benjamini-Hochberg) corrections
6. **Visualization**: Volcano plots, significance rankings, and distributions

## Output Files

### CSV Results
- **`Gene_Burden_Analysis_Results.csv`**: Complete results for all analyzed genes
- **`Significant_Genes_FDR_0.05.csv`**: Genes passing FDR < 0.05 threshold
- **`Significant_Genes_Raw_P_0.05.csv`**: Genes with nominal p < 0.05

### Key Output Columns
- `Gene`: Gene symbol/identifier
- `Cases_Burden`: Number of cases with ≥1 variant in this gene
- `Cases_No_Burden`: Number of cases without variants in this gene
- `Controls_Burden`: Expected number of controls with burden (rounded for testing)
- `Controls_No_Burden`: Expected number of controls without burden
- `Controls_Burden_Fraction`: Exact fractional expected control burden
- `P_Value`: Barnard's test p-value (or Fisher's fallback)
- `Odds_Ratio`: Haldane-corrected odds ratio estimate
- `CI_Lower`, `CI_Upper`: 95% confidence interval for odds ratio
- `P_Value_FDR`: FDR-corrected p-value
- `Significant_FDR`: Boolean flag for FDR significance

### Visualization Outputs
- **`Gene_Burden_Volcano_Plot_annotated.png`**: Volcano plot with targeted gene annotations
- **`Top_Significant_Genes_Styled.png`**: Bar chart of most significant genes
- **`Burden_Distribution_Styled.png`**: Histogram of case burden counts per gene

## Performance Characteristics

### Runtime Estimates
- ~30-60 seconds per 1,000 genes (12-core system)
- Barnard's test is computationally intensive but more powerful than alternatives
- Automatic load balancing across available cores

### Memory Requirements
- Moderate memory usage: 4-8 GB RAM recommended
- Efficient data.table operations for large datasets
- Automatic garbage collection between test batches

### Optimization Features
- Intelligent chunking of statistical tests across cores
- Early detection and handling of edge cases
- Robust error handling with informative messages

## Advantages over Fisher's Test

1. **Higher Statistical Power**: Barnard's test provides better power, especially for unbalanced designs
2. **No Conditioning**: Unconditional test doesn't assume fixed marginal totals
3. **Better Type I Error Control**: More accurate p-values in small sample scenarios
4. **Automatic Fallback**: Robust implementation with Fisher's test backup

## Troubleshooting

### Common Issues

**Missing DescTools Package**
```
Error: there is no package called 'DescTools'
```
Solution: The script auto-installs missing packages, but you can manually install:
```r
install.packages("DescTools", dependencies = TRUE)
```

**File Path Errors**
```
Error: SRR names file not found
```
Solution: Verify and update the three file path variables at the script top

**Statistical Test Failures**
```
Warning: Falling back to Fisher's test
```
This is normal behavior when Barnard's test encounters numerical issues

**Memory Limitations**
```
Error: cannot allocate vector of size X.X Gb
```
Solution: Reduce the number of parallel cores or process data in smaller batches

### Debugging Mode
The `VERBOSE <- TRUE` setting provides detailed progress information:
- File reading status
- Sample and gene counts
- Parallel processing status  
- Statistical test progress
- Output file locations

## Comparison with Other Methods

| Method | Power | Speed | Exact | Unconditional |
|--------|-------|--------|-------|---------------|
| **Barnard's** | High | Slow | Yes | Yes |
| Fisher's | Medium | Fast | Yes | No |
| Chi-squared | Medium | Very Fast | No | Yes |
| Bayesian | High | Very Slow | No | Yes |

## Best Practices

1. **Quality Control**: Ensure variant annotation quality before analysis
2. **Population Matching**: Use population frequencies matching your study ancestry
3. **Multiple Testing**: Always report FDR-corrected results for gene-level analyses
4. **Effect Size**: Consider both statistical significance and biological relevance (OR magnitude)
5. **Replication**: Validate significant findings in independent cohorts

## Citation

If you use this pipeline in your research, please cite:
- Barnard, G.A. (1947). Significance tests for 2×2 tables. *Biometrika*, 34(1/2), 123-138
- The DescTools R package for Barnard's test implementation
- Population frequency databases used (gnomAD, etc.)

## License

[Specify your license - e.g., MIT, GPL-3, etc.]

## Contact

[Your contact information]


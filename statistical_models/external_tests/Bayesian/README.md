# Bayesian Gene-Level Burden Analysis

A comprehensive R pipeline for performing Bayesian statistical analysis of gene-level variant burden in case-control genomic studies. This tool compares variant burden between cases and population-based controls using Bayesian inference with Jeffreys priors.

## Features

- **Bayesian Statistical Framework**: Uses Beta-Binomial conjugate priors (Jeffreys prior: α=0.5, β=0.5) for robust statistical inference
- **Population Control Estimation**: Estimates control burden from population frequency data (gnomAD, PopFreqMax) instead of requiring matched controls
- **Parallel Processing**: Optimized for multi-core execution with automatic core detection
- **Multiple Testing Correction**: Implements both Bonferroni and FDR (Benjamini-Hochberg) corrections
- **Comprehensive Visualization**: Generates volcano plots, significance bar charts, and burden distribution plots
- **Robust Data Handling**: Safe CSV reading with automatic column detection and data cleaning

## Requirements

### R Version
- R ≥ 4.0.0

### Required R Packages
```r
# Core packages (auto-installed if missing)
- parallel      # Multi-core processing
- data.table    # Fast data manipulation  
- dplyr         # Data wrangling
- readr         # CSV I/O
- ggplot2       # Plotting
- ggrepel       # Plot label positioning
- scales        # Plot scaling utilities
```

## Input Data Format

### File Structure
The pipeline expects the following directory structure:
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
Each variant file must contain:
- **Gene Name** (or similar gene identifier column)
- **Variant** (variant identifier - will be auto-generated if missing)
- **gnomAD_AF** (gnomAD allele frequency - optional, defaults to 0)
- **PopFreqMax** (maximum population frequency - optional, defaults to 0)

### Sample List File
Create a text file (`srr_names.txt`) containing one sample ID per line:
```
SRR_ID_1
SRR_ID_2
SRR_ID_3
...
```

## Configuration

### Key Parameters
Edit these variables at the top of the script:

```r
# Statistical parameters
NSIM_GLOBAL <- 1e6        # Bayesian simulation iterations
PRIOR_ALPHA <- 0.5        # Jeffreys prior alpha (cases)
PRIOR_BETA <- 0.5         # Jeffreys prior beta (cases)
PRIOR_ALPHA_CTRL <- 0.5   # Jeffreys prior alpha (controls)
PRIOR_BETA_CTRL <- 0.5    # Jeffreys prior beta (controls)

# Visualization parameters  
TOP_LABEL_N <- 22         # Number of top genes to label in plots
LOG10OR_CAP <- 6          # Cap for log10(OR) plotting range
NEGLOG10FDR_CAP <- 300    # Cap for -log10(FDR) plotting

# File paths
RESULTS_DIR <- "/path/to/output/directory"
```

### File Paths to Update
Before running, update these paths in the script:
```r
base_path <- "/path/to/results/directory"          # Your data directory
srr_names_file <- "/path/to/srr_names.txt"         # Your sample list file
RESULTS_DIR <- "/path/to/output/directory"         # Output directory
```

## Usage

### Basic Execution
```bash
Rscript burden_analysis.R
```

### Parallel Processing
The script automatically detects available CPU cores and uses up to 12 cores by default. To modify:
```r
num_cores <- min(YOUR_DESIRED_CORES, available_cores - 1)
```

## Methodology

### Statistical Approach
1. **Burden Calculation**: For each gene, counts samples with ≥1 qualifying variant
2. **Population Controls**: Estimates expected burden frequency using Hardy-Weinberg equilibrium from population allele frequencies
3. **Bayesian Inference**: 
   - Uses Beta-Binomial model with Jeffreys priors
   - Samples from posterior distributions of case vs. control rates
   - Calculates odds ratios and credible intervals
4. **Statistical Testing**: Computes two-sided Bayesian p-values and applies multiple testing corrections

### Population Frequency Calculation
For each gene with variants having allele frequencies AF₁, AF₂, ..., AFₙ:
```
Individual burden probability = 1 - ∏(1 - min(2×AFᵢ, 1))
Expected control burden = sample_size × burden_probability
```

## Output Files

### CSV Results
- **`Gene_Burden_Analysis_Results.csv`**: Complete results for all genes
- **`Significant_Genes_FDR_0.05.csv`**: FDR-significant genes only
- **`Significant_Genes_Raw_P_0.05.csv`**: Nominally significant genes (p < 0.05)

### Key Output Columns
- `Gene`: Gene symbol
- `Cases_Burden`: Number of cases with variants in this gene  
- `Controls_Burden`: Expected number of controls with burden
- `Posterior_Median_OR`: Median odds ratio from Bayesian inference
- `OR_CI_Lower`, `OR_CI_Upper`: 95% credible interval for OR
- `Posterior_Prob_OR_gt_1`: Posterior probability that OR > 1
- `Bayes_p`: Two-sided Bayesian p-value
- `P_Value_FDR`: FDR-corrected p-value
- `Significant_FDR`: Boolean flag for FDR significance

### Visualization Outputs
- **`Gene_Burden_Volcano_Plot_annotated.png`**: Volcano plot with gene annotations
- **`Top_Significant_Genes_Styled.png`**: Bar chart of most significant genes
- **`Burden_Distribution_Styled.png`**: Histogram of gene burden counts

## Performance Notes

### Memory Requirements
- Recommended: 8-16 GB RAM for typical datasets
- Memory usage scales with: (number of genes) × (simulation iterations)
- Large datasets may require reducing `NSIM_GLOBAL`

### Runtime Estimates
- ~1-5 minutes per 1000 genes (12-core system)
- Scales approximately linearly with gene count
- Bayesian sampling is the computational bottleneck

### Optimization Tips
- Use SSD storage for faster file I/O
- Increase `CHUNK_SIZE` for more memory, fewer allocations
- Reduce `NSIM_GLOBAL` for faster (less precise) results
- Ensure adequate RAM to avoid swapping

## Troubleshooting

### Common Issues

**Missing Files Error**
```
Error: SRR names file not found
```
Solution: Verify `srr_names_file` path and file existence

**Memory Errors**
```
Error: cannot allocate vector of size X.X Gb  
```
Solution: Reduce `NSIM_GLOBAL` or increase available RAM

**No Significant Results**
```
Warning: No genes pass FDR threshold
```
Solution: Check input data quality, consider less stringent thresholds

**Parallel Processing Failures**
```
Warning: Failed to create parallel cluster
```
Solution: Script automatically falls back to sequential processing

### Debugging Mode
Enable verbose output for troubleshooting:
```r
VERBOSE <- TRUE
```

## Citation

If you use this pipeline in your research, please cite:
- The Bayesian statistical framework implemented
- Population frequency databases used (gnomAD, etc.)
- This analysis pipeline (provide GitHub link or DOI when available)

## License

[Specify your license here - e.g., MIT, GPL-3, etc.]

## Contact

[Your contact information]

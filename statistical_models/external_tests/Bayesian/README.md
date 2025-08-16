# Bayesian Gene-Level Burden Analysis

A comprehensive R pipeline for performing gene-level burden analysis using Bayesian methodology with Jeffreys priors, featuring high-quality publication-ready visualizations.

## Overview

This pipeline analyzes genetic variant burden at the gene level using a fully Bayesian approach by:
- Loading variant data from multiple samples
- Calculating gene-level burden for each sample  
- Estimating population-level burden frequencies using exact product method
- Performing Bayesian hypothesis testing with Beta-Binomial models
- Using Jeffreys priors (α = β = 0.5) for objective analysis
- Applying multiple testing corrections (Bonferroni and FDR)
- Generating publication-quality plots with vector and high-resolution raster formats

## Key Features

- **Bayesian Statistics**: Full Bayesian inference with Beta-Binomial conjugate models
- **Jeffreys Priors**: Objective, non-informative priors (α = β = 0.5)
- **MCMC Sampling**: High-precision Monte Carlo simulation (default: 1M samples)
- **Fractional Controls**: Handles fractional control burden estimates naturally
- **Parallel Processing**: Multi-core MCMC sampling with chunked memory management
- **Robust Graphics**: Publication-ready plots with capped axes and selective labeling
- **Memory Management**: Efficient chunked sampling to handle large simulations

## Statistical Advantages

- **Handles Uncertainty**: Naturally incorporates uncertainty in control estimates
- **No Continuity Corrections**: Fractional control counts handled directly
- **Full Posterior**: Provides complete posterior distribution of odds ratios
- **Credible Intervals**: True Bayesian credible intervals (not confidence intervals)
- **Objective Priors**: Jeffreys priors minimize subjective bias

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
- Multi-core CPU strongly recommended (Bayesian sampling is compute-intensive)
- Sufficient RAM for MCMC sampling (8-16 GB recommended)
- System fonts (Arial preferred, falls back to system defaults)

## Configuration

### Key Parameters
Update these configuration variables at the top of the script:

```r
# Random seed for reproducibility
SEED <- 12345

# MCMC sampling parameters
NSIM_GLOBAL <- 1e6      # Total MCMC samples (1M default)
CHUNK_SIZE <- 5e4       # Memory management chunk size

# Bayesian priors (Jeffreys)
PRIOR_ALPHA <- 0.5      # Cases prior α
PRIOR_BETA  <- 0.5      # Cases prior β  
PRIOR_ALPHA_CTRL <- 0.5 # Controls prior α
PRIOR_BETA_CTRL  <- 0.5 # Controls prior β

# Plotting parameters
LOG10OR_CAP <- 6        # Cap log10(OR) at ±6 for visualization
NEGLOG10FDR_CAP <- 300  # Cap -log10(FDR) at 300
MIN_NONZERO_P <- 1e-300 # Floor for p-values to prevent log(0)

# Paths
RESULTS_DIR <- "/path/to/output/directory"
```

### File Paths
Update these paths in the main function:
```r
srr_file <- "/path/to/srr_names.txt"           # Sample names file
base_path <- "/path/to/results/directory"      # Sample data directory
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
# Source the script - runs automatically
source("gene_burden_analysis_bayesian.R")

# Results and plots saved to RESULTS_DIR
```

### Custom Configuration
```r
# Modify parameters before sourcing
NSIM_GLOBAL <- 5e5    # Reduce for faster testing
LOG10OR_CAP <- 4      # Different visualization range
source("gene_burden_analysis_bayesian.R")
```

## Output Files

### Results Tables
- `Gene_Burden_Analysis_Results.csv`: Complete Bayesian results for all genes
- `Significant_Genes_FDR_0.05.csv`: Genes significant after FDR correction (q < 0.05)
- `Significant_Genes_Raw_P_0.05.csv`: Genes with nominal significance (Bayes p < 0.05)

### Visualizations
All plots saved in multiple high-quality formats (PDF, EPS, TIFF@1200 DPI):

1. **Fig1_Gene_Burden_Volcano**: 
   - Volcano plot: log10(posterior median OR) vs -log10(FDR p-value)
   - Axes capped at ±6 for log10(OR) and 300 for -log10(FDR)
   - Highlights specific genes: HLA-DQA1, HLA-DRB1, HLA-DQB1, KMT2C, HLA-DRB5, PROM1

2. **Fig2_Top_Significant_Genes**: 
   - Horizontal bar plot of top significant genes
   - Ordered by Bayesian p-value significance

3. **Fig3_Burden_Distribution**: 
   - Histogram showing distribution of burden counts across genes

## Results Columns

| Column | Description |
|--------|-------------|
| `Gene` | Gene symbol/identifier |
| `Cases_Burden` | Number of cases with burden variants |
| `Cases_No_Burden` | Number of cases without burden variants |
| `Controls_Burden` | Estimated fractional controls with burden |
| `Controls_No_Burden` | Estimated fractional controls without burden |
| `Posterior_Median_OR` | Posterior median odds ratio |
| `OR_CI_Lower`, `OR_CI_Upper` | 95% Bayesian credible interval |
| `Posterior_Prob_OR_gt_1` | Posterior probability that OR > 1 |
| `Bayes_p` | Two-sided Bayesian p-value |
| `P_Value_Bonferroni` | Bonferroni-corrected p-value |
| `P_Value_FDR` | FDR-corrected p-value (Benjamini-Hochberg) |
| `Significant_Bonferroni` | Boolean: significant after Bonferroni |
| `Significant_FDR` | Boolean: significant after FDR correction |

## Statistical Methodology

### Bayesian Model
For each gene:
- **Cases**: X ~ Binomial(n_cases, p_cases)
- **Controls**: Y ~ Binomial(n_controls, p_controls) [fractional Y allowed]
- **Priors**: p_cases ~ Beta(0.5, 0.5), p_controls ~ Beta(0.5, 0.5)

### Posterior Inference
- **Posteriors**: p_cases|data ~ Beta(0.5 + X, 0.5 + n_cases - X)
- **OR Sampling**: OR = [p_cases/(1-p_cases)] / [p_controls/(1-p_controls)]
- **Bayesian p-value**: 2 × min(P(OR > 1), P(OR ≤ 1))

### Population Burden Estimation
```
For gene with variants [v1, v2, ..., vn]:
- Carrier probability per variant: min(2 × AF, 1)
- Gene burden probability: 1 - ∏(1 - carrier_prob_i)
- Expected burden count: sample_size × gene_burden_probability
```

## Performance Considerations

### Computational Intensity
- **MCMC Sampling**: 1M samples per gene can be time-intensive
- **Runtime**: 30 minutes to several hours depending on dataset size
- **Memory**: Uses chunked sampling to manage memory efficiently
- **Parallelization**: Scales well with multiple CPU cores

### Optimization Tips
```r
# For testing/debugging - reduce simulation size
NSIM_GLOBAL <- 1e5    # 100K samples (faster)

# For production - increase for higher precision
NSIM_GLOBAL <- 5e6    # 5M samples (slower but more precise)

# Chunk size affects memory usage
CHUNK_SIZE <- 1e4     # Smaller chunks = less memory
```

## Advantages over Frequentist Methods

1. **Natural Uncertainty Handling**: Incorporates all sources of uncertainty
2. **Fractional Controls**: No need for artificial rounding
3. **Intuitive Interpretation**: Posterior probabilities are directly interpretable
4. **No Asymptotic Assumptions**: Exact inference regardless of sample size
5. **Rich Output**: Full posterior distribution, not just point estimates

## Customization

### Prior Selection
```r
# Informative priors (if justified)
PRIOR_ALPHA <- 1.0; PRIOR_BETA <- 9.0  # Assumes low baseline burden

# Jeffreys priors (default - objective)
PRIOR_ALPHA <- 0.5; PRIOR_BETA <- 0.5  # Non-informative

# Uniform priors
PRIOR_ALPHA <- 1.0; PRIOR_BETA <- 1.0  # Flat prior
```

### Visualization
```r
# Adjust plotting ranges
LOG10OR_CAP <- 4        # Narrower OR range
NEGLOG10FDR_CAP <- 50   # Lower significance cap

# Modify gene labels
genes_to_label <- c("GENE1", "GENE2", "GENE3")  # Custom gene set
```

## Quality Control

### Convergence Diagnostics
- Large sample sizes (1M) ensure convergence
- Chunked sampling prevents memory issues
- Seed setting ensures reproducibility

### Error Handling
- Robust file reading with error recovery
- Memory management for large simulations
- Graceful handling of edge cases (zero burden genes)

## Troubleshooting

### Common Issues
1. **Long Runtime**: Reduce NSIM_GLOBAL for testing
2. **Memory Problems**: Reduce CHUNK_SIZE or number of cores
3. **Zero Burden Genes**: Handled automatically with prior regularization

### Performance Tuning
```r
# Fast testing configuration
NSIM_GLOBAL <- 1e4; CHUNK_SIZE <- 1e3

# High precision configuration  
NSIM_GLOBAL <- 5e6; CHUNK_SIZE <- 1e5
```

## Citation

When using this Bayesian pipeline, please cite:
- The Bayesian methodology and Beta-Binomial model approach
- Jeffreys prior justification for objective analysis
- Population genetics burden estimation methods

## License

This pipeline is provided as-is for research purposes. Please ensure compliance with individual tool licenses.

## Support

For questions or issues with this pipeline, please check the documentation of individual tools or create an issue in this repository.

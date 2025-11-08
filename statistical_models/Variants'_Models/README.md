# Variant Statistical Analysis Pipeline

A high-performance Python pipeline for conducting Fisher's exact test, Barnard's exact test, and Bayesian analysis on genetic variants using gnomAD population data.

## Features

- **Three Statistical Methods**:
  - Fisher's exact test (fast, conservative)
  - Barnard's exact test (more powerful, with intelligent downsampling for large datasets)
  - Bayesian analysis (Beta-Binomial model with credible intervals)

- **Optimized for Large Datasets**:
  - Multi-threaded parallel processing
  - Chunked processing with memory management
  - Intelligent Barnard test handling (exact, downsampled, or chi-squared approximation)
  - Progress tracking with `tqdm`

- **Homozygote-Aware**: Uses gnomAD homozygote counts when available for more accurate control carrier estimation

- **Multiple Testing Correction**: Automatically applies Bonferroni and FDR (Benjamini-Hochberg) corrections

## Requirements

```bash
pip install pandas numpy scipy statsmodels tqdm
```

**Minimum Requirements**:
- Python 3.7+
- scipy >= 1.9.0 (for `barnard_exact`)

## Input Format

Your input CSV file should contain the following columns:

- `Variant`: Variant identifier (e.g., chr1:12345:A:G)
- `SRR_Files`: Sample IDs carrying the variant (semicolon/comma/pipe separated)
- `File_Count`: Number of case samples carrying the variant (alternative to parsing SRR_Files)
- `gnomAD_Allele_Count`: Allele count from gnomAD
- `gnomAD_Allele_Number`: Allele number from gnomAD
- `gnomAD_Homozygotes`: Number of homozygous individuals (optional but recommended)

**Example**:
```csv
Variant,SRR_Files,File_Count,gnomAD_Allele_Count,gnomAD_Allele_Number,gnomAD_Homozygotes
chr1:12345:A:G,SRR001;SRR002,2,150,250000,0
chr2:67890:C:T,SRR003,1,50,250000,0
```

## Usage

### Basic Usage

1. Place your input CSV file in the same directory as the script or update the path in `CONFIG`:

```python
CONFIG = {
    'input_file': 'variant_analysis_full.csv',
    'total_cases': 62,  # Update this to your case sample size
    ...
}
```

2. Run the script:

```bash
python variant_analysis_pipeline.py
```

### Configuration

Edit the `CONFIG` dictionary at the top of the script to customize:

```python
CONFIG = {
    # File paths
    'input_file': 'variant_analysis_full.csv',
    'output_variant_fisher': 'variant_fisher.csv',
    'output_variant_barnard': 'variant_barnard.csv',
    'output_variant_bayesian': 'variant_bayesian.csv',
    'output_skipped': 'skipped_variants.csv',
    
    # Study parameters
    'total_cases': 62,  # Total number of case samples
    'significance_level': 0.05,
    
    # Column names (adjust if your CSV has different names)
    'gnomad_ac_col': 'gnomAD_Allele_Count',
    'gnomad_an_col': 'gnomAD_Allele_Number',
    'gnomad_hom_col': 'gnomAD_Homozygotes',
    
    # Performance tuning
    'chunk_size': 12,        # Variants per chunk
    'n_workers': 12,         # Parallel threads
    'flush_every_chunk': True,
    
    # Bayesian settings
    'bayesian_nsam': 1000000,      # Number of posterior samples
    'bayesian_prior_alpha': 0.5,   # Beta prior alpha
    'bayesian_prior_beta': 0.5,    # Beta prior beta
    
    # Barnard test limits
    'barnard_max_total': 10000,     # Max sample size for exact test
    'barnard_downsample_to': 5000,  # Downsample size for large controls
    'barnard_timeout': 180          # Timeout in seconds (Unix only)
}
```

## Output Files

### 1. Fisher's Exact Test (`variant_fisher.csv`)
- P-values from Fisher's exact test
- Odds ratios with 95% confidence intervals
- Bonferroni and FDR-corrected p-values

### 2. Barnard's Exact Test (`variant_barnard.csv`)
- P-values from Barnard's exact test
- Method used (exact, downsampled, chi2_approx, or failed)
- Same OR and corrections as Fisher's

### 3. Bayesian Analysis (`variant_bayesian.csv`)
- Posterior median odds ratio
- 95% credible intervals
- Posterior probability OR > 1
- Two-sided Bayesian p-value

### 4. Skipped Variants (`skipped_variants.csv`)
- Variants that couldn't be processed
- Reasons for skipping (missing data, invalid values, etc.)

## Output Columns

All output files include:
- `Variant`: Variant identifier
- `Cases_Carriers`: Number of cases with variant
- `Cases_NonCarriers`: Number of cases without variant
- `Controls_Est_Carriers`: Estimated control carriers
- `Controls_NonCarriers`: Control non-carriers
- `gnomAD_AC`, `gnomAD_AN`, `gnomAD_Homozygotes`: Input data
- `Control_N_indiv`: Total control individuals
- `Odds_Ratio`, `CI_Lower`, `CI_Upper`: OR and confidence intervals
- `P_Value`: Test-specific p-value
- `P_Value_Bonferroni`: Bonferroni-corrected p-value
- `P_Value_FDR`: FDR-corrected p-value
- `Significant_FDR`: Boolean flag for FDR significance

## Performance Tips

1. **Adjust thread count**: Set `n_workers` to match your CPU cores
2. **Tune chunk size**: Larger chunks = less overhead but more memory
3. **Barnard test optimization**: 
   - For very large datasets (>10,000 controls), exact Barnard test may be slow
   - The pipeline automatically uses downsampling or chi-squared approximation
4. **Memory management**: The script uses chunked processing and explicit garbage collection

## Statistical Methods

### Fisher's Exact Test
- Classic 2x2 contingency table test
- Conservative (conditional on margins)
- Fast and reliable

### Barnard's Exact Test
- More powerful than Fisher's test (unconditional)
- Computationally intensive for large samples
- Automatic handling:
  - **Exact**: Total N ≤ 10,000
  - **Downsampled**: Controls downsampled to 5,000
  - **Chi-squared**: Fallback for very large datasets

### Bayesian Analysis
- Beta-Binomial conjugate model
- Non-informative prior (α=0.5, β=0.5)
- Provides posterior distributions and credible intervals

## Control Carrier Estimation

The pipeline estimates control carriers using:

1. **If homozygote data available**:
   ```
   Heterozygotes = AC - 2×Homozygotes
   Carriers = Homozygotes + Heterozygotes
   ```

2. **Fallback (Hardy-Weinberg)**:
   ```
   p = AC/AN
   N = AN/2
   Carriers = N × (2p - p²)
   ```

## Platform Notes

- **Linux/Mac**: Full support including Barnard test timeout
- **Windows**: Timeout feature disabled for Barnard test (will still compute but may be slower)

## Troubleshooting

### Issue: Barnard test is too slow
**Solution**: Reduce `barnard_max_total` or increase `barnard_downsample_to` in CONFIG

### Issue: Out of memory errors
**Solution**: 
- Reduce `chunk_size`
- Reduce `bayesian_nsam`
- Reduce `n_workers`

### Issue: Missing scipy.stats.barnard_exact
**Solution**: Upgrade scipy to >= 1.9.0
```bash
pip install --upgrade scipy
```

### Issue: Many variants skipped
**Solution**: Check `skipped_variants.csv` for reasons. Common causes:
- Missing gnomAD data
- Invalid AC/AN values
- Zero allele numbers

## License

-This pipeline is provided as-is for research purposes. Please ensure compliance with individual tool licenses.

## Support
-For questions or issues with this pipeline, please check the documentation of individual tools or create an issue in this repository.

## Acknowledgments

- gnomAD for population allele frequency data
- scipy, pandas, and numpy development teams

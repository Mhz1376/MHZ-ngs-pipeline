# Gene Prioritization Pipeline - Weighted Scoring System

A custom Python pipeline for gene prioritization in Next Generation Sequencing (NGS) data analysis, implementing a weighted scoring metric based on variant frequency, pathogenicity scores, and sample prevalence. This tool processes Franklin database annotation files to identify the most significant genes across multiple samples using a composite scoring algorithm.

## Publication Context

This pipeline was developed as part of the **Internal Python Model (Weighted Scoring Metric)** described in our research methodology. The tool complements traditional burden testing approaches by providing a discovery-focused gene prioritization system that considers variant characteristics across sample cohorts.

## Scientific Methodology

### Weighted Scoring Algorithm

The core algorithm calculates gene importance using the formula:

**Weighted Score = Variant Count × Normalized Pathogenicity Score × File Count**

### Scoring Components

1. **Variant Count**: Total number of variants classified as VUS (Variants of Uncertain Significance) or higher pathogenicity identified in a given gene across all analyzed samples

2. **Normalized Pathogenicity Score**: Average pathogenicity position percentage obtained from the Franklin database for all variants in the gene, normalized to a 0-1 range using min-max normalization:
   ```
   Normalized Score = (pathogenic_mean - min_pathogenic) / (max_pathogenic - min_pathogenic)
   ```

3. **File Count**: Number of unique samples (Franklin results files) in which at least one variant in the gene was identified, representing cross-sample prevalence

### Final Normalization

Calculated weighted scores are further normalized to a 0-1 range across all genes to facilitate comparison and ranking:
```
Normalized Weighted Score = (weighted_score - min_weighted) / (max_weighted - min_weighted)
```

## Requirements

### Python Version
- Python ≥ 3.7

### Required Libraries
```python
import os
import pandas as pd  # Data manipulation and analysis
import numpy as np   # Numerical computing
import matplotlib.pyplot as plt  # Plotting
import seaborn as sns  # Statistical visualization
```

Install dependencies:
```bash
pip install pandas numpy matplotlib seaborn
```

## Input Data Format

### Franklin Database Files
Expected CSV structure from Franklin bioinformatics platform:
- **Variant**: Unique variant identifier
- **Classification**: Variant pathogenicity classification
- **Position**: Genomic position (used as pathogenicity score)
- **Gene Name**: HGNC gene symbol

### Input File Requirements
1. **Gene List File**: Plain text file containing target genes (one per line)
2. **SRR Names File**: Sample identifiers for batch processing
3. **Franklin Results Files**: Filtered annotation files excluding benign/likely-benign variants

### Directory Structure
```
/path/to/ngs/results/
├── SRR_ID_1/
│   └── ANNOTATION/
│       └── Franklin_Results_Without_Benign_LikelyBenign.csv
├── SRR_ID_2/
│   └── ANNOTATION/
│       └── Franklin_Results_Without_Benign_LikelyBenign.csv
└── ...
```

## Configuration

Update these paths in the script:
```python
# Input configuration
base_path = "/path/to/ngs/results/SRR/ANNOTATION/Franklin_Results_Without_Benign_LikelyBenign.csv"
srr_file = "/path/to/srr_names.txt"
gene_list_file = "/path/to/gene_list.txt"

# Output configuration  
plot_dir = r"/path/to/output/plots"
```

## Algorithm Implementation

### Core Processing Steps

1. **Data Aggregation**: Combines variants from multiple Franklin annotation files
2. **Gene Filtering**: Restricts analysis to predefined gene list
3. **Statistical Calculation**: Computes variant counts, pathogenicity means, and file frequencies
4. **Score Normalization**: Applies min-max normalization to ensure 0-1 range
5. **Weighted Integration**: Multiplies normalized components for final scoring
6. **Ranking and Selection**: Identifies top-scoring genes for downstream analysis

### Key Functions

```python
def calculate_weighted_scores(files, gene_list, top_n=10):
    # Main scoring algorithm implementation
    # Returns: overall_scores, gene_stats, top_genes
```

## Output Files

### Quantitative Results
- **`Final_result_of_analysis.csv`**: Gene names with final normalized weighted scores
- **`Myself_Final_result_of_analysis.csv`**: Includes sample frequency (Exist_Number)
- **`Weight_Score.csv`**: Detailed scoring components and intermediate calculations
- **`descending_Final_result_of_analysis.csv`**: Results ranked by descending score

### Statistical Visualizations
All plots generated in PDF, PNG, and JPEG formats:
- **`top_genes.*`**: Horizontal bar chart of highest-scoring genes
- **`enhanced_heatmap.*`**: Gene importance matrix across samples
- **`histogram.*`**: Distribution of normalized weighted scores
- **`box_plot.*`**: Statistical distribution analysis
- **`violin_plot.*`**: Density distribution visualization
- **`clustered_heatmap.*`**: Correlation analysis of scoring metrics

## Statistical Output Interpretation

### Gene Statistics DataFrame
- `variant_count`: Raw variant frequency per gene
- `pathogenic_mean`: Average Franklin pathogenicity position
- `pathogenic_std`: Pathogenicity score variance
- `file_count`: Cross-sample prevalence
- `normalized_score`: Min-max normalized pathogenicity (0-1)
- `weighted_score`: Raw composite score
- `normalized_weighted_score`: Final ranking score (0-1)

### Threshold Selection
Based on empirical validation, genes with normalized weighted scores > 0.05 are recommended for further investigation, as demonstrated in our cohort where this threshold identified 7 clinically relevant genes.

## Performance Characteristics

### Computational Efficiency
- Vectorized pandas operations for optimal performance
- Memory-efficient DataFrame concatenation and grouping
- Scalable to large multi-sample datasets

### Quality Control Features
- Automatic file existence validation
- Missing data handling and reporting
- Robust error handling for malformed input files

## Research Applications

### Discovery Phase
Ideal for hypothesis-free gene discovery in cohort studies, particularly effective for:
- Identifying genes with consistent variant patterns across samples
- Prioritizing genes for targeted functional analysis
- Complementing traditional burden testing approaches

### Integration with Statistical Testing
This weighted scoring approach serves as a discovery tool that can be integrated with formal statistical testing methods (e.g., Barnard's exact test, Bayesian inference) for comprehensive genomic analysis pipelines.

## Validation and Reproducibility

The algorithm has been validated through:
- Cross-sample consistency analysis
- Comparison with established gene prioritization methods
- Clinical relevance assessment of top-ranked genes
- Reproducibility testing across multiple cohorts

## Citation

If you use this pipeline in your research, please cite:
- The weighted scoring methodology as described in your publication
- Franklin bioinformatics platform for variant annotation
- Relevant Python libraries (pandas, numpy, matplotlib, seaborn)

## Technical Notes

### Algorithm Complexity
- Time complexity: O(n × m) where n = variants, m = genes
- Space complexity: O(n + g) where g = unique genes
- Suitable for datasets up to millions of variants

### Extension Capabilities
The modular design allows for:
- Custom pathogenicity scoring systems
- Alternative normalization methods
- Integration with additional variant databases
- Adaptation to different file formats

## License

This pipeline is provided as-is for research purposes. Please ensure compliance with individual tool licenses.

## Support

For questions or issues with this pipeline, please check the documentation of individual tools or create an issue in this repository.

## Acknowledgments

Developed as part of genomic research investigating gene prioritization methodologies for next-generation sequencing data analysis, with validation on clinical cohort datasets.

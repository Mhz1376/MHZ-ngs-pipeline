# MHZ-ngs-pipeline

A comprehensive NGS processing and variant analysis pipeline developed for the research published in *"Inherited Antigen-Presentation Failure Through Germline HLA Variants: A Fundamental Mechanism in Breast Cancer Predisposition"*.

## Overview

This repository contains the complete bioinformatics workflow used to analyze next-generation sequencing data and identify germline  variants associated with breast cancer predisposition. The pipeline integrates quality control, alignment, variant calling, and statistical modeling to investigate the relationship between  variants and cancer susceptibility.

## Repository Structure

### `ngs_pipeline/`
Contains the complete NGS data processing workflow including:
- Raw data quality control and preprocessing
- Reference genome alignment
- Variant calling and filtering
- Annotation
- Quality metrics and reporting

Detailed documentation of each processing step is provided within the directory.

### `statistical_models/`
Includes four comprehensive statistical models used in the analysis:
- Model 1: [Internal]
- Model 2: [External: Fisher] 
- Model 3: [External: Barnard]
- Model 4: [External: Bayesian]

Each model includes source code, documentation, and validation scripts.

## Data Availability

### Primary Dataset
The complete dataset is publicly available on Zenodo:
**DOI:** https://doi.org/10.5281/zenodo.16811086

### Local Access
Data can also be downloaded directly from this repository. **Important:** Please read the README file located in the data directory before use, and check the `version` file for software version requirements.

## Requirements

See the `version` file in the data directory for complete software version specifications.

## Citation

If you use this pipeline or data in your research, please cite me.

## Contributing

We welcome contributions, bug reports, and feature requests. Please open an issue or submit a pull request.

## License

[Include your license information]

## Contact

For questions, collaboration requests, or bug reports regarding this repository, please contact:

**Mohammad Hosein Zabihi**  
📧 Mohammad.hosein.zabihi@gmail.com

**Subject Line:** Please include "MHZ-ngs-pipeline" in your email subject to help us track and respond to messages efficiently.

---

*This pipeline was developed as part of research into inherited antigen-presentation mechanisms in breast cancer predisposition.*

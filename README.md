# Nasal Microbiome Analysis Pipeline for Dupilumab Treatment Study

## Overview

This repository contains the bioinformatics pipeline and analysis scripts used in the paper "Dupilumab Treatment is Associated with Clinical Improvement and a Shift towards a health-associated Nasal Passage Microbiota in Diffuse Type 2 Chronic Rhinosinusitis" by Ryser et al.

## Study Description

This study investigates the relationship between Dupilumab treatment and changes in the nasal microbiome in patients with Diffuse Type 2 Chronic Rhinosinusitis. The analysis pipeline processes 16S rRNA gene sequencing data to characterize microbial communities and their shifts during treatment.

## Pipeline Features

This analysis utilizes a modified version of the [16S_dada2_snakemake_pipeline](https://github.com/TomasDemeter/16S_dada2_snakemake_pipeline.git) with the following key features:

- DADA2-based processing of raw sequencing reads
- Quality control and filtering of sequence data
- ASV (Amplicon Sequence Variant) generation
- Taxonomic classification using SILVA database
- Comprehensive statistical analysis and visualization

## Requirements

- Snakemake 8.16
- Conda
- R (version 4.0 or higher)
- DADA2 R package
- Additional R dependencies (specified in environment files)

## Installation

1. Clone this repository:

```bash
git clone [repository-url]
```

2. Install required environments:

```bash
cd workflow/envs
conda env create -f snakemake.yml
```

## Usage

### Data Preparation

1. Place raw sequencing data in `raw_reads/<experiment_name>/`
2. Update configuration in `config/config.yml`
3. Prepare metadata file according to template

### Pipeline Execution

```bash
cd workflow
conda activate snakemake
snakemake -s Snakefile.py --profile profiles/default --cores <number_of_cores>
```

## Analysis Notebooks

The repository includes several Jupyter notebooks for different aspects of the analysis:

- `data_cleanup.ipynb`: Initial data processing and quality control
- `differential_abundance.ipynb`: Analysis of microbial abundance changes
- `alpha_beta.ipynb`: Diversity metrics analysis
- `ancombc2.ipynb`: Compositional data analysis
- `clinical_correlations.ipynb`: Correlation with clinical parameters

## Citation

If you use this pipeline, please cite:


## Contributors

- Fabio S. Ryser
- Tomas Demeter
- Judith Bergada Pijuan
- Srikanth Mairpady Shambat
- Catrin Brühlmann
- Tina Mauthe
- Markus Hilty
- Michael B. Soyka
- Urs C. Steiner
- Silvio D. Brugger

## License

GNU General Public License v3.0

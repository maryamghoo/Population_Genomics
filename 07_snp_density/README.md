# SNP Density Analysis Pipeline

## Overview
This workflow calculates genome-wide Single-nucleotide polymorphism SNP density using vcftools and visualizes SNP density distributions and genomic patterns using R.

The pipeline includes:

- SNP density calculation across fixed genomic windows
- Summary statistics calculation
- SNP density distribution visualization
- Genome-wide SNP density plotting

See:
`snp_density.sh`
`snp_density.R`

### Directory Structure
```
project/
│
├── r_data/
│   └── snp_density.snpden
│
├── r_results/
│
├── snp_density.sh
├── snp_density.R
└── README.md
```

## 1. SNP Density Calculation Workflow

### Purpose
The bash workflow calculates SNP density across the genome using fixed-size genomic windows.

The script performs:

- Read filtered SNP VCF file
- Divide genome into fixed windows
- Count SNPs within each window
- Generate SNP density table

### Expected Input
Filtered SNP dataset
`snps_beforeLD.vcf`

### Output
`snp_density.snpden`
`snp_density.log`

### Usage
- Adjust SNP density window size based on genome size, marker density, and downstream visualization goals.
- Update the following:
```
/path/to/snps_beforeLD.vcf
/path/to/output/directory/<file_name>
```

## 2. R SNP Density Visualization Workflow

### Purpose
The R workflow visualizes SNP density distributions and genome-wide SNP density patterns.

The script performs:

- Import SNP density data
- Calculate mean, median, and standard deviation
- Generate SNP density distribution plots
- Calculate cumulative genome positions
- Generate genome-wide SNP density plots

### R Packages
The script requires:

- `ggplot2`
- `dplyr`

### Expected Input
SNP density file generated from the previous step.
`r_data/snp_density.snpden`

### Output
`snp_density_distribution.png`
`genome_snp_density.png`
saved in: `r_results/`

### Usage
- This script was run using `vcftools v0.1.17`. Ensure `vcftools` is installed and available in your environment.
- Adjust SNP density window sizes, smoothing parameters, and figure dimensions as needed for your dataset and visualization goals.
- Update the following in `snp_density.R` if needed:
`data_dir <- "r_data"`
`output_dir <- "r_results"`

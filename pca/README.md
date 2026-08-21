# PCA Analysis Pipeline

## Overview
This workflow performs Principal component analysis (PCA) on filtered SNP datasets using PLINK and visualizes population structure using R.

The pipeline includes:

- PCA generation from genome-wide SNP datasets
- Site-specific PCA analyses using sample subsets
- Conversion of eigenvalue/eigenvector outputs into CSV format
- PCA visualization and figure generation in R

See:
`pca.sh`
`pca.R`

## Directory Structure
```
project/
│
├── bash_data/
│   ├── mainland_bc_samples.txt
│   └── cbc_sbc_samples.txt
│
├── bash_results/
│   ├── all_regions.eigenval
│   ├── all_regions.eigenvec
│   ├── mainland_bc.eigenval
│   ├── mainland_bc.eigenvec
│   ├── cbc_sbc.eigenval
│   └── cbc_sbc.eigenvec
│
├── r_data/
│   ├── all_regions.eigenvalues.csv
│   ├── all_regions.eigenvectors.csv
│   ├── mainland_bc.eigenvalues.csv
│   ├── mainland_bc.eigenvectors.csv
│   ├── cbc_sbc.eigenvalues.csv
│   └── cbc_sbc.eigenvectors.csv
│
├── r_results/
│
├── pca.sh
├── pca.R
└── README.md
```

## 1. Bash PCA Workflow

### Purpose
The bash workflow performs PCA analyses from filtered SNP datasets and prepares outputs for downstream visualization in R.

The script performs:

- PCA on all individuals
- PCA on selected sampling subsets
- Generate `.eigenval` and `.eigenvec` files
- Convert PCA outputs into CSV format for R

### Expected Inputs
- Input SNP dataset
`snps_afterLD.vcf`
- Sample subset files in case of performing PCA on subpopulations. These files contain sample IDs.

### Output
- PLINK PCA outputs
`all_regions.eigenval`
`all_regions.eigenvec`
`mainland_bc.eigenval`
`mainland_bc.eigenvec`
`cbc_sbc.eigenval`
`cbc_sbc.eigenvec`
- CSV files for R
`all_regions.eigenvalues.csv`
`all_regions.eigenvectors.csv`
`mainland_bc.eigenvalues.csv`
`mainland_bc.eigenvectors.csv`
`ccbc_sbc.eigenvalues.csv`
`cbc_sbc.eigenvectors.csv`

### Usage
- This script was run using `plink v1.9`. Ensure `plink` is installed and available in your environment.
- Adjust file paths, output names, and sample subset files based on your project structure and sampling design.
- Update the following in `pca.sh`:
```
/path/to/plink
/path/to/snps_afterLD.vcf
/path/to/output/
/path/to/<site>_samples.txt
/path/to/output/<site>.eigenvalues.csv
/path/to/output/<site>.eigenvectors.csv
```

## 2. R PCA Visualization Workflow

### Purpose
The R workflow generates PCA visualizations from PLINK PCA outputs.

The script performs:

- PCA axis rotation
- Scree plot generation
- Population-specific PCA visualization
- Ellipse visualization around clusters
- Combined figure generation
- Export to PDF, PNG, and SVG formats

### R Packages
The script requires:

- `ggplot2`
- `cowplot`
- `RColorBrewer`

### Expected Input
- CSV PCA outputs generated from the bash workflow:
`r_data/all_regions.eigenvalues.csv`
`r_data/all_regions.eigenvectors.csv`
`r_data/mainland_bc.eigenvalues.csv`
`r_data/mainland_bc.eigenvectors.csv`
`r_data/cbc_sbc.eigenvalues.csv`
`r_data/cbc_sbc.eigenvectors.csv`
	
### Output
`pca_all.pdf`
`pca_bc.pdf`
`pca_cbc_sbc.pdf`
`pca_combined.pdf`
as well as PNG and SVG versions in `r_results/`. SVG format is useful for opening figures in vector graphic editors such as Inkscape, allowing manual adjustment and editing of individual image elements without loss of quality.

### Usage
- Adjust PCA rotation angles, color palettes, factor ordering, and variance labels as needed for your dataset and visualization goals.
- Update the following in `pca.R` if needed:
`data_dir <- "r_data"`
`output_dir <- "r_results"`

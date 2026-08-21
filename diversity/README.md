# Genetic Diversity Analysis Pipeline

## Overview
This workflow estimates and compares genetic diversity among genetic clusters using nucleotide diversity (π), observed heterozygosity (Ho), and expected heterozygosity (He).

The pipeline includes:

- Estimation of nucleotide diversity (π) in non-overlapping 50-kb windows using VCFtools
- Estimation of observed and expected heterozygosity using VCFtools
- Filtering nucleotide diversity windows containing fewer than 10 SNPs
- Calculation of mean nucleotide diversity for each genetic cluster
- Calculation of observed heterozygosity (Ho)
- Calculation of expected heterozygosity (He)
- Correction of expected heterozygosity for sampling bias
- Statistical comparison among genetic clusters using Kruskal-Wallis tests
- Pairwise comparisons using Wilcoxon rank-sum tests with Bonferroni correction
- Automatic generation of significance-group letters
- Visualization of nucleotide diversity and observed heterozygosity using boxplots

See:

`diversity.sh`  
`diversity.R`

### Directory Structure
```
project/
│
├── r_data/
│   ├── `<cluster_name>`.windowed.pi
│   └── `<cluster_name>`.het
│
├── r_results/
│
├── genetic_diversity.sh
├── genetic_diversity.R
└── README.md
```

## 1. VCFtools Genetic Diversity Workflow

### Purpose
The VCFtools workflow calculates nucleotide diversity and individual heterozygosity for each sampling region/genetic cluster.

Two SNP datasets are used:

- **Nucleotide diversity (π):** autosomal biallelic SNPs before LD pruning
- **Heterozygosity:** autosomal biallelic SNPs after LD pruning

Unrelated individuals are used for both analyses.

### Expected Input

- for nucleotide diversity, autosomal biallelic SNP dataset before LD pruning:
```
snps_beforeLD_unrelated.vcf
```
- for heterozygosity, autosomal biallelic SNP dataset after LD pruning:
```
snps_afterLD_unrelated.vcf
```
- sample Lists
A separate sample list is required for each genetic cluster. For example:
```
`<cluster_name>.txt`
```
Each file should contain the sample IDs belonging to that genetic cluster, with one sample ID per line. The sample IDs must exactly match the sample names in the VCF file.

### Output
VCFtools commands generate a `.windowed.pi` and a `.het` file for each genetic cluster. For example:
```
pi_`<cluster_name>`.windowed.pi
het_`<cluster_name>`.het
```

### Usage
- Update the following paths for your dataset:
```
/path/to/snps_beforeLD_unrelated.vcf
/path/to/snps_afterLD_unrelated.vcf
/path/to/<cluster_name>.samples.txt
/path/to/output/
```
- Run the VCFtools commands separately for each genetic cluster.
- Move or copy the resulting `.windowed.pi` and `.het` files into `r_data/` using the filenames expected by the R script.

## 2. R Genetic Diversity Analysis

### Purpose
The R workflow summarizes genetic diversity estimates, statistically compares genetic clusters, and generates a combined figure showing nucleotide diversity and observed heterozygosity.

The script performs:

- Import and combine nucleotide diversity files
- Remove π windows containing fewer than 10 SNPs
- Calculate mean, median, and standard deviation of π
- Import and combine heterozygosity files
- Calculate observed heterozygosity (Ho)
- Calculate expected heterozygosity (He)
- Correct expected heterozygosity for sampling bias
- Calculate summary statistics for each genetic cluster
- Perform Kruskal-Wallis tests
- Perform pairwise Wilcoxon rank-sum tests
- Apply Bonferroni correction for multiple comparisons
- Generate significance-group letters automatically
- Generate boxplots of π and Ho
- Combine the plots into a single figure
- Export statistical results, summary tables, and figures

### R Packages
The script requires:

- `dplyr` — data manipulation, grouping, filtering, joining, and summary statistics
- `ggplot2` — plotting and figure export
- `patchwork` — combining the π and heterozygosity plots
- `multcompView` — generating significance letters from pairwise comparisons

### Expected Input
can be found in `r_data/` folder:
```
`<cluster_name>`.windowed.pi
`<cluster_name>`.het
```

### Output
```
pi_summary.txt
pi_kruskal_wallis.txt
pi_pairwise_wilcoxon.txt
pi_significance_letters.txt
heterozygosity_summary.txt
heterozygosity_individuals.txt
heterozygosity_kruskal_wallis.txt
heterozygosity_pairwise_wilcoxon.txt
heterozygosity_significance_letters.txt
pi_heterozygosity.pdf
pi_heterozygosity.png
pi_heterozygosity.pdf
pi_heterozygosity.png
pi_heterozygosity.svg
```
saved in: `r_results/`. SVG format is useful for opening figures in vector graphic editors such as Inkscape, allowing manual adjustment and editing of individual image elements without loss of quality.

### Usage
- This script was run using `vcftools v0.1.17`. Ensure `vcftools` is installed and available in your environment.
- Adjust the genetic clusters, file names, plotting order, SNP-per-window threshold, or colors change, as needed for your dataset and visualization goals in `diversity.R`
```
input_dir <- "r_data"
output_dir <- "r_results"
min_variants <- 10
clusters <- c(
  "Merritt Island",
  "SL",
  "CBC",
  "SBC",
  "Archbold",
  "Ocala"
)
file_names <- c(
  "Merritt Island" = "merritt",
  "SL"             = "sl",
  "CBC"            = "cbc",
  "SBC"            = "sbc",
  "Archbold"       = "archbold",
  "Ocala"          = "ocala"
)
```

# BayeScan Selection Analysis Pipeline

## Overview
This workflow identifies candidate SNPs potentially under selection using BayeScan. It consists of three bash scripts followed by an R analysis:

1. `01_vcf_to_bayescan.sh`
2. `02_bayescan.sh`
3. `03_prepare_bayescan_results.sh`
4. `bayescan_analysis.R`

The workflow starts from an LD-pruned autosomal biallelic SNP VCF containing unrelated individuals, converts the VCF into BayeScan input format, runs BayeScan, combines BayeScan results with SNP information from the original VCF, and identifies and visualizes candidate outlier SNPs in R.

## Purpose
This workflow is used to detect SNPs showing unusually high genetic differentiation among populations, which may represent candidate loci under selection.

The workflow performs:

- Convert a VCF file into BayeScan input format
- Assign individuals to populations using a manually prepared population file
- Run BayeScan on the converted SNP dataset
- Combine BayeScan results with chromosome, position, SNP ID, reference allele, and alternate allele information from the original VCF
- Map reference sequence names to chromosomes
- Examine the distribution of SNP FST values
- Identify candidate outlier SNPs using FDR
- Generate a Manhattan plot of SNP FST values
- Highlight candidate outlier SNPs in the Manhattan plot


## Step 1: Convert VCF to BayeScan Format

See: `01_vcf_to_bayescan.sh`

This script uses `vcf2bayescan.pl` to convert the SNP VCF into the input format required by BayeScan. Population assignments are provided using a manually prepared sample-to-population file.

`vcf2bayescan.pl` was obtained from: https://github.com/BrianFerolito/BayeScan

### Expected Inputs
- LD-pruned autosomal biallelic SNP VCF containing unrelated individuals:
`snps_afterLD_unrelated.vcf`
- Manually prepared sample population file:
`sample_populations.txt`
The population file associates each sample with its corresponding population. For example:
```
sample_01    Population1
sample_02    Population1
sample_03    Population2
sample_04    Population2
```
The sample IDs should correspond to sample IDs in the input VCF. Adjust the population-file format and `-colNum` option as required by `vcf2bayescan.pl`.
- `vcf2bayescan.pl`

### Output
`bayescan_input.txt`

### Usage

- Download or obtain `vcf2bayescan.pl` from the source repository listed above.
- Prepare `sample_populations.txt` manually according to the population assignments used in the study.
- Update the following in the script:
```
VCF="/path/to/snps_afterLD_unrelated.vcf"
POP_FILE="/path/to/sample_populations.txt"
OUTFILE="/path/to/bayescan_input.txt"
/path/to/vcf2bayescan.pl
```
- Adjust the SLURM settings based on the cluster configuration and dataset size.
- Submit the job:
`sbatch 01_vcf_to_bayescan.sh`


## Step 2: Run BayeScan

See: `02_bayescan.sh`

This script runs BayeScan on the input file generated in Step 1 to estimate locus-specific differentiation and identify candidate loci potentially affected by selection.

### Expected Input
`bayescan_input.txt`

### Output
BayeScan generates multiple output files in the specified output directory: 
`bayescan_input_fst.txt`
`bayescan_input_AccRte.txt`
`bayescan_input_Verif.txt`
`bayescan_input.sel`
The main file subsequently used in this workflow is: `bayescan_input_fst.txt`

### Usage
- This script was run using `bayescan v2.1`. Ensure `bayescan` is installed and available in your environment.
- Update the following in the script:
```
INDIR="/path/to/bayescan"
OUTDIR="/path/to/bayescan/bayescan_output"
/path/to/bayescan_2.1
```
- Adjust BayeScan parameters as appropriate for the dataset and study design.
The current script uses:
```
-pr_odds 10000
-n 5000
-thin 20
-nbp 10
-pilot 2000
-burn 20000
```
- Adjust the SLURM settings, including CPUs, memory, runtime, and log paths, based on the cluster configuration and dataset size.
- Submit the job:
`sbatch 02_bayescan.sh`


## Step 3: Prepare BayeScan Results

See: `03_prepare_bayescan_results.sh`

This script combines SNP information from the original VCF with the corresponding BayeScan results.

SNP chromosome, position, ID, reference allele, and alternate allele information are extracted from the VCF using `bcftools query` and combined with the BayeScan FST output.

### Expected Inputs
- Original SNP VCF used for BayeScan:
`snps_afterLD_unrelated.vcf`
- BayeScan FST output:
`bayescan_input_fst.txt`

### Output
`final_bayescan_fst.txt`
The output contains:
```
CHR
POS
ID
REF
ALT
prob
log10(PO)
qval
alpha
fst
```

### Usage
- This script was run using `bcftools v1.20`. Ensure `bcftools` is installed and available in your environment.
- Ensure that the SNP order in the BayeScan results corresponds to the SNP order in the VCF used to generate the BayeScan input.
- Update the following in the script:
```
VCF="/path/to/snps_afterLD_unrelated.vcf"
BAYESCAN_DIR="/path/to/bayescan/bayescan_output"
OUTPUT="/path/to/final_bayescan_fst.txt"
```
- Adjust the SLURM settings based on the cluster configuration.
- Submit the job:
`sbatch 03_prepare_bayescan_results.sh`


## Step 4: BayeScan Analysis and Visualization in R

See: `bayescan_analysis.R`

This script processes the prepared BayeScan results, maps reference sequence names to chromosomes, identifies candidate outlier SNPs, and generates FST distribution and Manhattan plots.

The script performs:

- Read the prepared BayeScan results
- Read chromosome mapping information
- Map reference sequence names to chromosome numbers and labels
- Remove scaffolds or sequences without chromosome assignments for chromosome-specific analyses
- Calculate the number of SNPs per chromosome
- Plot the distribution of FST values
- Identify candidate outlier SNPs using FDR < 0.05
- Export candidate outlier SNP tables
- Generate a Manhattan plot of SNP FST values
- Highlight FDR outliers in the Manhattan plot
- Export plots in PNG and PDF formats

### R Packages

The script requires:

- `ggplot2`
- `dplyr`
- `qqman`

### Expected Inputs
- Prepared BayeScan results in `bash_results/`:
`final_bayescan_fst.txt`
- Manually prepared chromosome mapping file in `r_data/`:
`chromosomes.txt` which was created manually using information from the reference genome assembly and genome annotation. It maps reference sequence identifiers to chromosome labels and chromosome numbers used for plotting.
These columns are assigned in R as:
```
Chr_Label
Chr_Num
CHR
```
The values in the `CHR` column must match the chromosome/reference sequence names present in `final_bayescan_fst.txt`.

### Output
The R script generates the files below in `r_results/`:
```
fst_distribution.png
fst_chr_distribution.png
fdr05_chr_outliers.txt
fdr05_outliers.txt
bayescan_manhattan.png
bayescan_manhattan.pdf
```

`fdr05_outliers.txt` contains candidate outlier SNPs with FDR < 0.05 across the complete dataset.
`fdr05_chr_outliers.txt` contains FDR < 0.05 candidate SNPs restricted to sequences successfully mapped to chromosomes.

### Usage
- Create `chromosomes.txt` manually based on the reference genome and annotation information.
- Ensure that reference sequence identifiers in `chromosomes.txt` exactly match those in the BayeScan result table.
- Adjust plotting parameters, FDR threshold, colors, figure dimensions, and other visualization settings as needed.
- Update the following in `bayescan_analysis.R` if needed:
```
data_dir <- "r_data"
bash_dir <- "bash_results"
output_dir <- "r_results"
```


## General Usage Notes
- Run the scripts **in order**:
```
01_vcf_to_bayescan.sh
02_bayescan.sh
03_prepare_bayescan_results.sh
bayescan_analysis.R
```
- Adjust all SLURM settings (partition, CPUs, memory, runtime, log paths, etc.) based on the cluster configuration and dataset size.
- Update all `/path/to/...` paths before submitting the bash jobs.
- Ensure the required software is available, including Perl, `bcftools`, BayeScan, and the required R packages.
- `sample_populations.txt` must be prepared manually according to the population assignments used for the analysis.
- `chromosomes.txt` must be prepared manually using the reference genome assembly and annotation information.
- Check output and error logs after each bash step before continuing to the next step.
- BayeScan parameters and the FDR threshold should be adjusted when appropriate for the dataset and study design.

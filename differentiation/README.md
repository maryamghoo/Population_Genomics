# FST and Isolation-by-Distance Analysis Pipeline

## Overview
This workflow estimates pairwise genetic differentiation among genetic clusters using FST, evaluates the significance of observed FST values using permutation tests, and tests for isolation by distance (IBD).

The pipeline includes:

- Pairwise windowed FST estimation
- Permutation testing of pairwise FST
- Calculation of genome-wide average FST
- Comparison of observed and permuted FST values
- Calculation of empirical P-values
- Calculation of genetic and geographic distances
- Mantel test for isolation by distance

See:
`fst_permutation.sh`
`fst.R`
`isolation_by_distance.R`

### Directory Structure
```
project/
│
├── bash_results/
│   ├── fst_results/
│   │   └── *.windowed.weir.fst
│   │
│   └── permutations/
│       ├── *.average_fst_permutation.csv
│       ├── *.shuffled.txt
│       ├── perm_<cluster_1>_<cluster_2>.txt
│       └── perm_<cluster_2>_<cluster_1>.txt
│
├── r_data/
│   └── coordinates.csv
│
├── r_results/
│   ├── pairwise_fst_summary.csv
│   ├── fst_permutation_pvalues.csv
│   ├── fst_permutations.png
│   ├── fst_weighted_permutations.png
│   └── mantel_test_results.csv
│
├── fst_permutation.sh
├── fst.R
├── isolation_by_distance.R
└── README.md
```

## 1. Bash FST and Permutation Workflow

### Purpose
The bash workflow calculates pairwise windowed FST between genetic clusters and performs permutation tests by randomly reassigning individuals between each pair of clusters while preserving the original sample sizes.

The script performs:

- Calculate observed pairwise windowed FST
- Determine the original sample size of each cluster
- Randomly shuffle individuals between two clusters
- Preserve the original cluster sample sizes after shuffling
- Calculate windowed FST for each permutation
- Calculate average weighted and mean FST for each permutation
- Save permutation results in CSV format

### Expected Input
- Autosomal biallelic SNP dataset before LD pruning after removal of related individuals:
`snps_beforeLD_unrelated.vcf`
- Sample ID files for each genetic cluster. Each file should contain the sample IDs belonging to that genetic cluster, with one sample ID per line.
`<cluster_name>.txt`

### Output
- Observed pairwise FST in `bash_results/fst_results`:
`<cluster_1>.<cluster_2>.windowed.weir.fst`
- Permutation summary `bash_results/permutations`:
`<cluster_1>.<cluster_2>.average_fst_permutation.csv`
- Temporary shuffled sample list in `bash_results/permutations/`:
`<cluster_1>.<cluster_2>.shuffled.txt`
- Temporary permuted cluster sample lists in `bash_results/permutations/`:
`perm_<cluster_1>_<cluster_2>.txt`  
`perm_<cluster_2>_<cluster_1>.txt`

The temporary `.txt` files are overwritten during each permutation, so only the files from the final permutation remain after the job is completed. 
The `.average_fst_permutation.csv` file contains the summarized FST results from all permutations and is used for downstream analysis in R.
The permutation CSV contains:
```
permutation
average_weighted_fst
average_mean_fst
```

### Important Notes
- The script should be run for each pairwise comparison by changing `POP1_NAME` and `POP2_NAME`.
- For six genetic clusters, there are 15 unique pairwise comparisons.
- The current workflow uses 50-kb FST windows with a 25-kb step size.
- The current workflow performs 100 permutations for each pairwise comparison.
- Sample lists used for the analysis should correspond to individuals present in the input VCF file.

### Usage
- This script was run using `vcftools v0.1.17`. Ensure `vcftools` is installed and available in your environment.
- Adjust SLURM job settings (partition, CPUs, memory, runtime, output paths, etc.) based on your cluster configuration and dataset size.
- Adjust the number of permutations and FST window size/step as needed for your study.
- Update the following in `fst_permutation.sh`:
```
/path/to/logs/fst_permutation.out
/path/to/logs/fst_permutation.err
/path/to/snps_beforeLD_unrelated.vcf
/path/to/sample_lists
/path/to/fst_results
<site_1>
<site_2>
```
- Submit the job:
  `sbatch fst_permutation.sh`

## 2. R FST Analysis Workflow

### Purpose
The R workflow summarizes observed pairwise FST results and compares them with the permutation distributions.

The script performs:

- Import all pairwise windowed FST files
- Extract cluster names from filenames
- Combine pairwise FST results
- Calculate genome-wide average weighted and mean FST
- Import permutation results
- Plot permutation distributions
- Add observed FST values to permutation plots
- Calculate empirical P-values
- Export summary tables and figures

### R Packages
The script requires:

- `dplyr`
- `ggplot2`
- `patchwork`

### Expected Input
Can be found in `bash_results/`:
- Observed FST files in `bash_results/fst_results/`:
`<cluster_1>.<cluster_2>.windowed.weir.fst`
- Permutation files in `bash_results/permutations/`:
`<cluster_1>.<cluster_2>.average_fst_permutation.csv`

### Output
`pairwise_fst_summary.csv`
`fst_permutation_pvalues.csv`
`fst_permutations.png`
`fst_weighted_permutations.png`
saved in: `r_results/`.

`pairwise_fst_summary.csv` is also used as input for the subsequent isolation-by-distance analysis.

### Usage
- Ensure that observed FST and permutation filenames follow the expected naming format so that cluster names can be extracted automatically.
- Ensure that the same cluster-pair naming is used in observed and permutation results.
- Update the following in `fst.R` if needed:
```
fst_dir <- file.path("bash_results", "fst_results")
permutation_dir <- file.path("bash_results", "permutations")
output_dir <- "r_results"
```
- Plot dimensions, histogram bin width, and other plotting parameters can be adjusted as needed.

## 3. R Isolation-by-Distance Workflow

### Purpose
The isolation-by-distance workflow tests whether genetic differentiation among genetic clusters is associated with geographic distance.

The script performs:

- Import pairwise weighted FST values generated by `fst.R`
- Construct a symmetric pairwise FST matrix
- Transform FST to genetic distance using `FST / (1 - FST)`
- Import sampling coordinates for genetic clusters
- Calculate mean geographic coordinates for each genetic cluster
- Calculate pairwise geographic distances
- Convert geographic distances from meters to kilometers
- Match the order of genetic clusters between genetic and geographic distance matrices
- Perform a Mantel test using Pearson correlation and 9,999 permutations
- Export Mantel test results

### R Packages
The script requires:

- `dplyr`
- `vegan`
- `geosphere`

### Expected Input
- Pairwise FST summary generated by `fst.R`:
`r_results/pairwise_fst_summary.csv`
- Sampling coordinate data:
`r_data/coordinates.csv`
The coordinate file should contain:
```
Cluster
LatY
LatX
```

### Output
`mantel_test_results.csv`
saved in: `r_results/`.
The output contains the Mantel correlation coefficient, P-value, and number of permutations.

### Usage
- Run `fst.R` before `isolation_by_distance.R`, because the isolation-by-distance analysis uses `pairwise_fst_summary.csv` generated by the FST analysis.
- Ensure that genetic cluster names in `coordinates.csv` exactly match the cluster names used in the pairwise FST results.
- Adjust the number of Mantel permutations if needed.
- Update the following in `isolation_by_distance.R` if needed:
```
results_dir <- "r_results"
data_dir <- "r_data"
```

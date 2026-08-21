#!/bin/bash -l
#SBATCH -J prepare_GO
#SBATCH -o /prepare_GO/logs/prepare_GO.out
#SBATCH -e /prepare_GO/logs/prepare_GO.err
#SBATCH -p normal
#SBATCH --cpus-per-task=1
#SBATCH -t 02:00:00

set -euo pipefail

# Define directories
DATA_DIR="bash_data"
OUTDIR="bash_results"

mkdir -p "${OUTDIR}"

# Define input files
ANNOTATION="${DATA_DIR}/fdr_outliers_annotated.fsj_multianno.txt"
GFF3="/path/to/genome_annotation.gff3"


# 1. Extract candidate genes from ANNOVAR results

awk -F'\t' 'NR > 1 {
    n = split($7, genes, ";")
    for (i = 1; i <= n; i++) {
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", genes[i])
        if (genes[i] != "" && genes[i] != "." && genes[i] != "NONE")
            print genes[i]
    }
}' "${ANNOTATION}" | \
sort -u \
> "${OUTDIR}/candidate_genes.txt"


# 2. Extract background genes from the genome annotation

awk -F'\t' '$0 !~ /^#/ && $3 == "gene" {
    if (match($9, /ID=g[0-9]+/)) {
        gene = substr($9, RSTART + 3, RLENGTH - 3)
        print gene
    }
}' "${GFF3}" | \
sort -u \
> "${OUTDIR}/all_genes.txt"


# 3. Create gene-to-GO mapping

perl -ne '
    next unless /GO:\d+/;

    my ($gene) = /(g\d+)/;
    next unless defined $gene;

    my @go_terms = /GO:\d+/g;

    for my $go (@go_terms) {
        print "$gene\t$go\n";
    }
' "${GFF3}" | \
sort -u \
> "${OUTDIR}/gene2go_map.txt"


# 4. Identify candidate genes missing from the background

comm -23 \
  "${OUTDIR}/candidate_genes.txt" \
  "${OUTDIR}/all_genes.txt" \
  > "${OUTDIR}/candidate_genes_missing_from_background.txt"


# 5. Report basic counts

echo "Candidate genes:"
wc -l "${OUTDIR}/candidate_genes.txt"

echo "Background genes:"
wc -l "${OUTDIR}/all_genes.txt"

echo "Gene-GO associations:"
wc -l "${OUTDIR}/gene2go_map.txt"

echo "Candidate genes missing from background:"
wc -l "${OUTDIR}/candidate_genes_missing_from_background.txt"

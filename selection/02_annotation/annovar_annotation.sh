#!/bin/bash -l
#SBATCH -J annovar_annotation
#SBATCH -o /annovar_annotation/logs/annovar_annotation.out
#SBATCH -e /annovar_annotation/logs/annovar_annotation.err
#SBATCH -p normal
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00

set -euo pipefail


# Define directories
DATA_DIR="bash_data"
OUTDIR="bash_results"

mkdir -p "${OUTDIR}"

# Define input files
OUTLIERS="${DATA_DIR}/fdr05_outliers.txt"
CHR_MAP="${DATA_DIR}/gff3_ref_map.txt"

GFF3="/path/to/genome_annotation.gff3"
REFERENCE_FASTA="/path/to/reference_genome.fna"

# Define software and ANNOVAR database
ANNOVAR="/path/to/annovar"
GFF3_TO_GENEPRED="/path/to/gff3ToGenePred"
HUMANDB="${ANNOVAR}/humandb"

BUILD="fsj"

# 1. Convert FDR outlier SNPs to ANNOVAR input format
awk 'BEGIN{OFS="\t"} {print $2,$4,$4,$5,$6}' \
  "${OUTLIERS}" \
  > "${OUTDIR}/fdr_outliers.avinput"

# 2. Prepare GFF3 annotation for ANNOVAR
#
# If chromosome/scaffold names in the GFF3 differ from those in the
# reference genome, rename them using CHR_MAP.
# Also retain only the first Name= value when multiple names occur.

awk 'BEGIN{FS=OFS="\t"}
     NR==FNR {map[$1]=$2; next}
     /^#/ {print; next}
     {if ($1 in map) $1=map[$1]; print}' \
  "${CHR_MAP}" \
  "${GFF3}" | \
awk 'BEGIN{FS=OFS="\t"}
     /^#/ {print; next}
     {sub(/Name=([^;,]+),[^;]+/, "Name=\\1"); print}' \
  > "${OUTDIR}/${BUILD}_annovar_ready.gff3"

# 3. Convert GFF3 to refGene format
"${GFF3_TO_GENEPRED}" \
  "${OUTDIR}/${BUILD}_annovar_ready.gff3" \
  /dev/stdout | \
tr -d '\r' \
  > "${HUMANDB}/${BUILD}_refGene.txt"

# 4. Prepare reference FASTA headers
sed '/^>/ s/ .*$//' \
  "${REFERENCE_FASTA}" \
  > "${OUTDIR}/${BUILD}_reference_short_headers.fna"

# 5. Generate transcript sequences
"${ANNOVAR}/retrieve_seq_from_fasta.pl" \
  -format refGene \
  -seqfile "${OUTDIR}/${BUILD}_reference_short_headers.fna" \
  -out "${HUMANDB}/${BUILD}_refGeneMrna.fa" \
  "${HUMANDB}/${BUILD}_refGene.txt"

# 6. Annotate FDR outlier SNPs
"${ANNOVAR}/table_annovar.pl" \
  "${OUTDIR}/fdr_outliers.avinput" \
  "${HUMANDB}" \
  -buildver "${BUILD}" \
  -out "${OUTDIR}/fdr_outliers_annotated" \
  -remove \
  -protocol refGene \
  -operation g \
  -nastring .
  
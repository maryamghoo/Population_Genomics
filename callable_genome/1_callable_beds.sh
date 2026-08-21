#!/bin/bash
#SBATCH -J callable_beds
#SBATCH -o /callable_beds/logs/callable_beds.out
#SBATCH -e /callable_beds/logs/callable_beds.err
#SBATCH -p normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -t 3-00:00:00

# Stop on errors, undefined variables, and failed pipes
set -euo pipefail

##### Settings #####

MAXJOBS=4
MIN_DEPTH=2
MAX_DEPTH=21

BAMLIST="/path/to/bamlist.txt"
SEX_LIST="/path/to/sex_chromosomes.txt"

OUTTSV="/path/to/callable_sizes.tsv"
BEDDIR="/path/to/callable_beds"

##########

mkdir -p "$BEDDIR"

# output header
echo -e "ID\tcallable_mb" > "$OUTTSV"

# loop over all BAMs
while read -r bam; do

    [[ -z "$bam" ]] && continue

    # wait if MAXJOBS jobs already running
    while [ "$(jobs -r | wc -l)" -ge "$MAXJOBS" ]; do
        sleep 5
    done

    (
        sample=$(basename "$bam" .bam)
        echo "Processing $sample"

        bedtools genomecov -ibam "$bam" -bga \
            | awk -v MIN="$MIN_DEPTH" -v MAX="$MAX_DEPTH" '$4 >= MIN && $4 <= MAX {print $1"\t"$2"\t"$3}' \
            | sort -k1,1 -k2,2n \
            | bedtools merge -i - \
            | awk 'NR==FNR{sex[$1]; next} !($1 in sex)' "$SEX_LIST" - \
            > "${BEDDIR}/${sample}.callable.bed"

        callable_mb=$(awk '{s+=($3-$2)} END{printf "%.6f", s/1e6}' \
            "${BEDDIR}/${sample}.callable.bed")

        echo -e "${sample}\t${callable_mb}" >> "$OUTTSV"

    ) &

done < "$BAMLIST"

wait

echo "Callable region calculation complete."

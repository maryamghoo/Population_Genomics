#!/bin/bash
#SBATCH -J common_callable
#SBATCH -o /common_callable/logs/common_callable.out
#SBATCH -e /common_callable/logs/common_callable.err
#SBATCH -p normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH -t 2-00:00:00

# Stop on errors, undefined variables, and failed pipes
set -euo pipefail

##### settings #####

INPUT="/path/to/callable_multiinter.bed"

COMMON="/path/to/callable_common.bed"
MERGED="/path/to/callable_common.merged.bed"

THRESHOLD=218

##########

echo "Step 1: extracting regions callable in >= ${THRESHOLD} individuals"

awk -v T="$THRESHOLD" '$4>=T {print $1"\t"$2"\t"$3}' "$INPUT" > "$COMMON"

echo "Step 2: sorting and merging intervals"

sort -k1,1 -k2,2n "$COMMON" | bedtools merge -i - > "$MERGED"

echo "Step 3: sanity checks"

echo "Number of regions BEFORE merging:"
wc -l "$COMMON"

echo "Number of regions AFTER merging:"
wc -l "$MERGED"

echo "File size of merged BED:"
du -sh "$MERGED"

echo "Average region size (bp):"
awk '{sum+=($3-$2); n++} END{if(n>0) printf "%.0f\n", sum/n; else print 0}' "$MERGED"

echo "Step 4: calculating callable genome size (denominator)"

awk '{s+=($3-$2)} END{printf "Callable genome size: %.2f Mb\n", s/1e6}' "$MERGED"

echo "Finished."

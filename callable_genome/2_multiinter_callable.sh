#!/bin/bash
#SBATCH -J multiinter_callable
#SBATCH -o /multiinter_callable/logs/multiinter_callable.out
#SBATCH -e /multiinter_callable/logs/multiinter_callable.err
#SBATCH -p normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 2-00:00:00

# Stop on errors, undefined variables, and failed pipes
set -euo pipefail

##### Settings #####

BEDLIST="/path/to/bedlist.txt"

OUT="/path/to/callable_multiinter.bed"

##########

echo "Starting bedtools multiinter..."

mapfile -t BEDS < "$BEDLIST"

bedtools multiinter \
    -i "${BEDS[@]}" \
    > "$OUT"

echo "Finished."

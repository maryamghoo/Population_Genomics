#!/bin/bash
#SBATCH -J prepare_bayescan
#SBATCH -o /prepare_bayescan/logs/prepare_bayescan.out
#SBATCH -e /prepare_bayescan/logs/prepare_bayescan.err
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=8G

# Paths
VCF="/path/to/snps_afterLD_unrelated.vcf"
BAYESCAN_DIR="/path/to/bayescan/bayescan_output"

OUTPUT="/path/to/final_bayescan_fst.txt"

# Combine VCF SNP information with BayeScan results
{
    echo -e "CHR\tPOS\tID\tREF\tALT\tprob\tlog10(PO)\tqval\talpha\tfst"

    paste \
        <(bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' "$VCF") \
        <(tail -n +2 "$BAYESCAN_DIR/bayescan_input_fst.txt" | \
          awk '{$1=""; sub(/^[[:space:]]+/, ""); print}')

} > "$OUTPUT"

echo "Finished."
echo "Output: $OUTPUT"

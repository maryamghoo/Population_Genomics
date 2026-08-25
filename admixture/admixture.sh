#!/bin/bash -l
#SBATCH -J admixture
#SBATCH -o /admixture_results/logs/admixture.out
#SBATCH -e /admixture_results/logs/admixture.err
#SBATCH -p normal
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4000
#SBATCH -t 48:00:00
#SBATCH --array=1-10

cd /path/to/output_directory

K=$SLURM_ARRAY_TASK_ID

# Only task 1 converts VCF to PLINK format
if [ "$K" -eq 1 ]; then

	/path/to/plink \
  		--vcf /path/to/snps_afterLD.vcf \
  		--allow-extra-chr \
  		--make-bed \
  		--out adinput
  		
  	# Save original BIM
    cp adinput.bim adinput_original.bim

    # Replace chromosome column with integers
    awk '{print NR, $2, $3, $4, $5, $6}' \
      adinput.bim > adinput_rename.bim

    # Use renamed BIM for ADMIXTURE
    mv adinput_rename.bim adinput.bim

    # Signal that preparation is completely finished
    touch adinput.ready
fi

# Other array tasks wait until preparation is finished
while [ ! -f adinput.ready ]; do
    sleep 10
done

# Run ADMIXTURE 
/path/to/admixture \
-B \
--cv=10 \
-j10 \
adinput.bed $K | tee log${K}.out 

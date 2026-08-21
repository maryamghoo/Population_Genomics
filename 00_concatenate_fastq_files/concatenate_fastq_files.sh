#!/bin/bash
#SBATCH -J concatenate_fastq_files
#SBATCH -o /concat_results/logs/SAMPLE_ID.out
#SBATCH -e /concat_results/logs/SAMPLE_ID.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH -t 48:00:00


input_dir=/path/to/fastq_files
output_dir=/path/to/output_directory
SAMPLE_NAME="SAMPLE_ID"

cat ${input_dir}/${SAMPLE_NAME}_L00[0-4]_R1_001.fastq.gz > ${output_dir}/${SAMPLE_NAME}_R1.fastq.gz
cat ${input_dir}/${SAMPLE_NAME}_L00[0-4]_R2_001.fastq.gz > ${output_dir}/${SAMPLE_NAME}_R2.fastq.gz

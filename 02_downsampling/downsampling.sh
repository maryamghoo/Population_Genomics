#!/bin/bash
#SBATCH -J downsampling
#SBATCH -o /downsampling_results/downsampling.out
#SBATCH -e /downsampling_results/downsampling.err
#SBATCH -p normal
#SBATCH --cpus-per-task=32
#SBATCH -t 48:00:00
#SBATCH --mem-per-cpu=8000


# Define the maximum number of reads
max_reads=250000000

# Read the list of BAM files
input_list="/path/to/downsampling_target_bamlist.txt"

# Define the output directory for subsampled BAM files
output_dir="/path/to/downsampled/directory"

# Loop through each BAM file and downsample it
while IFS= read -r bam_file; do
  # Calculate the downsample fraction
  frac=$(samtools idxstats "$bam_file" | cut -f3 | awk -v max_reads="$max_reads" 'BEGIN {total=0} {total += $1} END {frac=max_reads/total; if (frac > 1) {print 1} else {print frac}}')
  
  # Get the base name of the BAM file
  bam_basename=$(basename "$bam_file")
  
  # Downsample the BAM file
  output_bam="$output_dir/${bam_basename%.bam}_subsampled.bam"
  samtools view -bs "$frac" "$bam_file" > "$output_bam"
  
  # Index the downsampled BAM file
  samtools index "$output_bam"
  
  # Generate .bamstats
  bamstats_output="$output_dir/${bam_basename%.bam}_subsampled.bamstats"
  samtools stats "$output_bam" > "$bamstats_output"
  
  # Sort the downsampled BAM file
  sorted_bam="$output_dir/${bam_basename%.bam}_subsampled.sorted.bam"
  samtools sort "$output_bam" -o "$sorted_bam"
  
  # Index the sorted BAM file
  samtools index "$sorted_bam"
  
  echo "Processed $bam_file: Downsampled to $output_bam, sorted to $sorted_bam, generated $bamstats_output"
done < "$input_list"

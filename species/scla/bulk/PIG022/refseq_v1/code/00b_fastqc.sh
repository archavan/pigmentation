#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=fastqc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=50G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load FastQC


cd /home/arc78/palmer_scratch/data/bulk-rna/PIG022
mkdir -p fastqc

# Loop through all fastq.gz files in subdirectories
find ./01_fastq -type f -name "*.fastq.gz" | while read fq
do
    # Extract parent folder name and filename
    parent=$(basename $(dirname "$fq"))
    file=$(basename "$fq")
    
    # Create a unique output name
    outname="${parent}_${file}"
    
    # Run FastQC and rename output
    fastqc --outdir=fastqc "$fq"
    
    # Move and rename the output files
    mv "fastqc/${file%.fastq.gz}_fastqc.html" "fastqc/${outname%.fastq.gz}_fastqc.html"
    mv "fastqc/${file%.fastq.gz}_fastqc.zip" "fastqc/${outname%.fastq.gz}_fastqc.zip"
done

#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_S_1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/S_1

cd /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/01_fastq

# run kallisto
kallisto quant \
-i /home/arc78/genomes/species/scla/refseq_v1/results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx \
-t 25 -b 100 -o /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/S_1 \
AS_ds.1deff6413a154cbbbde579e99b69a75c/AS_S25_L001_R1_001.fastq.gz \
AS_ds.1deff6413a154cbbbde579e99b69a75c/AS_S25_L001_R2_001.fastq.gz \
AS_ds.930e88549f914c0dae41ec09f28916ae/AS_S25_L001_R1_001.fastq.gz \
AS_ds.930e88549f914c0dae41ec09f28916ae/AS_S25_L001_R2_001.fastq.gz \
AS_ds.bc08778cd73e473d95a9b982d0541f3e/AS_S25_L001_R1_001.fastq.gz \
AS_ds.bc08778cd73e473d95a9b982d0541f3e/AS_S25_L001_R2_001.fastq.gz \
AS_ds.f577ed3d0ac04f90a89ad9ee8a48573a/AS_S25_L001_R1_001.fastq.gz \
AS_ds.f577ed3d0ac04f90a89ad9ee8a48573a/AS_S25_L001_R2_001.fastq.gz


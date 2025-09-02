#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_I2_2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I2_2

cd /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/01_fastq

# run kallisto
kallisto quant \
-i /home/arc78/genomes/species/scla/refseq_v1/results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx \
-t 25 -b 100 -o /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I2_2 \
B2_ds.3d1ae8c6551443a89c84cf9eb4558eb1/B2_S43_L001_R1_001.fastq.gz \
B2_ds.3d1ae8c6551443a89c84cf9eb4558eb1/B2_S43_L001_R2_001.fastq.gz \
B2_ds.a24c79f8229844abac3ba22b62b4a6be/B2_S43_L001_R1_001.fastq.gz \
B2_ds.a24c79f8229844abac3ba22b62b4a6be/B2_S43_L001_R2_001.fastq.gz \
B2_ds.c42600e3148a47078563b22ca5978265/B2_S43_L001_R1_001.fastq.gz \
B2_ds.c42600e3148a47078563b22ca5978265/B2_S43_L001_R2_001.fastq.gz \
B2_ds.e47bb3cc7eb245d586a1c48a18867a7a/B2_S43_L001_R1_001.fastq.gz \
B2_ds.e47bb3cc7eb245d586a1c48a18867a7a/B2_S43_L001_R2_001.fastq.gz


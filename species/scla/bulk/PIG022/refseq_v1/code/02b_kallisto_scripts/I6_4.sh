#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_I6_4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I6_4

cd /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/01_fastq

# run kallisto
kallisto quant \
-i /home/arc78/genomes/species/scla/refseq_v1/results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx \
-t 25 -b 100 -o /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I6_4 \
E6_ds.910a429249154316b78a6ad5da291b67/E6_S10_L001_R1_001.fastq.gz \
E6_ds.910a429249154316b78a6ad5da291b67/E6_S10_L001_R2_001.fastq.gz \
E6_ds.a33d76c2d6cd4812978159249e1fff1b/E6_S10_L001_R1_001.fastq.gz \
E6_ds.a33d76c2d6cd4812978159249e1fff1b/E6_S10_L001_R2_001.fastq.gz \
E6_ds.b43bf3f405c34af09e9df7cb0f528351/E6_S10_L001_R1_001.fastq.gz \
E6_ds.b43bf3f405c34af09e9df7cb0f528351/E6_S10_L001_R2_001.fastq.gz \
E6_ds.e03f4c4a69554f21b0bd7a4b2479ee15/E6_S10_L001_R1_001.fastq.gz \
E6_ds.e03f4c4a69554f21b0bd7a4b2479ee15/E6_S10_L001_R2_001.fastq.gz


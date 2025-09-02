#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_I4_2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I4_2

cd /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/01_fastq

# run kallisto
kallisto quant \
-i /home/arc78/genomes/species/scla/refseq_v1/results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx \
-t 25 -b 100 -o /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I4_2 \
B4_ds.765d0b675b8c4573b1953caadafa5d18/B4_S47_L001_R1_001.fastq.gz \
B4_ds.765d0b675b8c4573b1953caadafa5d18/B4_S47_L001_R2_001.fastq.gz \
B4_ds.8aab8fc79c734a1e89b7143fc242c2b1/B4_S47_L001_R1_001.fastq.gz \
B4_ds.8aab8fc79c734a1e89b7143fc242c2b1/B4_S47_L001_R2_001.fastq.gz \
B4_ds.a3218369a89444788114ce5f68f89caa/B4_S47_L001_R1_001.fastq.gz \
B4_ds.a3218369a89444788114ce5f68f89caa/B4_S47_L001_R2_001.fastq.gz \
B4_ds.f3f5aaf6f69940e48507981788cc1f18/B4_S47_L001_R1_001.fastq.gz \
B4_ds.f3f5aaf6f69940e48507981788cc1f18/B4_S47_L001_R2_001.fastq.gz


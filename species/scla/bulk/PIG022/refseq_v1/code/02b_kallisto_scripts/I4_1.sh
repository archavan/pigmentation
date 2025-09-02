#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_I4_1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I4_1

cd /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/01_fastq

# run kallisto
kallisto quant \
-i /home/arc78/genomes/species/scla/refseq_v1/results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx \
-t 25 -b 100 -o /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I4_1 \
A4_ds.3c994c321d02446b9c30fd5b90116443/A4_S33_L001_R1_001.fastq.gz \
A4_ds.3c994c321d02446b9c30fd5b90116443/A4_S33_L001_R2_001.fastq.gz \
A4_ds.3d159c8e7323458093665a86f077209b/A4_S33_L001_R1_001.fastq.gz \
A4_ds.3d159c8e7323458093665a86f077209b/A4_S33_L001_R2_001.fastq.gz \
A4_ds.a7fb31cb89a34722821d206a7ee9b607/A4_S33_L001_R1_001.fastq.gz \
A4_ds.a7fb31cb89a34722821d206a7ee9b607/A4_S33_L001_R2_001.fastq.gz \
A4_ds.fdd41365a9a14e13a32968d1adb0fad8/A4_S33_L001_R1_001.fastq.gz \
A4_ds.fdd41365a9a14e13a32968d1adb0fad8/A4_S33_L001_R2_001.fastq.gz


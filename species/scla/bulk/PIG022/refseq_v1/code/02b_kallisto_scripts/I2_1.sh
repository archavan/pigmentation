#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_I2_1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I2_1

cd /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/01_fastq

# run kallisto
kallisto quant \
-i /home/arc78/genomes/species/scla/refseq_v1/results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx \
-t 25 -b 100 -o /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I2_1 \
A2_ds.0d701a927f294ed8b393107437d4c096/A2_S29_L001_R1_001.fastq.gz \
A2_ds.0d701a927f294ed8b393107437d4c096/A2_S29_L001_R2_001.fastq.gz \
A2_ds.3c6e39e14ff349d2be075a72f38e6e72/A2_S29_L001_R1_001.fastq.gz \
A2_ds.3c6e39e14ff349d2be075a72f38e6e72/A2_S29_L001_R2_001.fastq.gz \
A2_ds.8c45bf84ae3243c4aa02dee816ebe787/A2_S29_L001_R1_001.fastq.gz \
A2_ds.8c45bf84ae3243c4aa02dee816ebe787/A2_S29_L001_R2_001.fastq.gz \
A2_ds.d014246a0efa4b0d899a8a88893ad89d/A2_S29_L001_R1_001.fastq.gz \
A2_ds.d014246a0efa4b0d899a8a88893ad89d/A2_S29_L001_R2_001.fastq.gz


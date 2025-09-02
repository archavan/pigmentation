#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_I1_4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I1_4

cd /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/01_fastq

# run kallisto
kallisto quant \
-i /home/arc78/genomes/species/scla/refseq_v1/results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx \
-t 25 -b 100 -o /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I1_4 \
E1_ds.3163f3baba8840018e77cdea7f7f10de/E1_S14_L001_R1_001.fastq.gz \
E1_ds.3163f3baba8840018e77cdea7f7f10de/E1_S14_L001_R2_001.fastq.gz \
E1_ds.620e1b0732b748fba491e95e21c678a8/E1_S14_L001_R1_001.fastq.gz \
E1_ds.620e1b0732b748fba491e95e21c678a8/E1_S14_L001_R2_001.fastq.gz \
E1_ds.6263abb211dc427597ca407d4a6a4bcb/E1_S14_L001_R1_001.fastq.gz \
E1_ds.6263abb211dc427597ca407d4a6a4bcb/E1_S14_L001_R2_001.fastq.gz \
E1_ds.e730e04445484a3b9ee195d9e9ae16aa/E1_S14_L001_R1_001.fastq.gz \
E1_ds.e730e04445484a3b9ee195d9e9ae16aa/E1_S14_L001_R2_001.fastq.gz


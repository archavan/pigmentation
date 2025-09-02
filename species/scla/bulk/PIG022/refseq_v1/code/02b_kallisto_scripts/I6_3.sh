#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_I6_3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I6_3

cd /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/01_fastq

# run kallisto
kallisto quant \
-i /home/arc78/genomes/species/scla/refseq_v1/results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx \
-t 25 -b 100 -o /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I6_3 \
D6_ds.1fecd9f978fd49c08a9a1c96a84e7efd/D6_S22_L001_R1_001.fastq.gz \
D6_ds.1fecd9f978fd49c08a9a1c96a84e7efd/D6_S22_L001_R2_001.fastq.gz \
D6_ds.25fdbdb2ab5f4d0fa60f9ac7d8e0eb5d/D6_S22_L001_R1_001.fastq.gz \
D6_ds.25fdbdb2ab5f4d0fa60f9ac7d8e0eb5d/D6_S22_L001_R2_001.fastq.gz \
D6_ds.7f3701bdf8644dab84b77670431f8543/D6_S22_L001_R1_001.fastq.gz \
D6_ds.7f3701bdf8644dab84b77670431f8543/D6_S22_L001_R2_001.fastq.gz \
D6_ds.c4c8c5d80a0543088d4ae0db26b0282b/D6_S22_L001_R1_001.fastq.gz \
D6_ds.c4c8c5d80a0543088d4ae0db26b0282b/D6_S22_L001_R2_001.fastq.gz


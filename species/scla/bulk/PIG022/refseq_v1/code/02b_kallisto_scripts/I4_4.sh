#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_I4_4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I4_4

cd /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/01_fastq

# run kallisto
kallisto quant \
-i /home/arc78/genomes/species/scla/refseq_v1/results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx \
-t 25 -b 100 -o /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I4_4 \
E4_ds.2a75e0d79f9348b7885474632b5445d3/E4_S23_L001_R1_001.fastq.gz \
E4_ds.2a75e0d79f9348b7885474632b5445d3/E4_S23_L001_R2_001.fastq.gz \
E4_ds.40faa9c57d8744c399ba95f3bcb39a25/E4_S23_L001_R1_001.fastq.gz \
E4_ds.40faa9c57d8744c399ba95f3bcb39a25/E4_S23_L001_R2_001.fastq.gz \
E4_ds.a075d87297154fa4958815427c72c2bc/E4_S23_L001_R1_001.fastq.gz \
E4_ds.a075d87297154fa4958815427c72c2bc/E4_S23_L001_R2_001.fastq.gz \
E4_ds.ca45ba8ef7754ced87ceb2abbea0d724/E4_S23_L001_R1_001.fastq.gz \
E4_ds.ca45ba8ef7754ced87ceb2abbea0d724/E4_S23_L001_R2_001.fastq.gz


#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_S_3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/S_3

cd /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/01_fastq

# run kallisto
kallisto quant \
-i /home/arc78/genomes/species/scla/refseq_v1/results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx \
-t 25 -b 100 -o /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/S_3 \
DS_ds.13aa695f6b394c33b795658f21e73420/DS_S9_L001_R1_001.fastq.gz \
DS_ds.13aa695f6b394c33b795658f21e73420/DS_S9_L001_R2_001.fastq.gz \
DS_ds.1b904fd7a5ca4901a9b5cc3753bd1d47/DS_S9_L001_R1_001.fastq.gz \
DS_ds.1b904fd7a5ca4901a9b5cc3753bd1d47/DS_S9_L001_R2_001.fastq.gz \
DS_ds.9405817ff6f843ee8af0c2eb69d4a3a5/DS_S9_L001_R1_001.fastq.gz \
DS_ds.9405817ff6f843ee8af0c2eb69d4a3a5/DS_S9_L001_R2_001.fastq.gz \
DS_ds.d7b863fa32d147a5b0b36434b8d3bb5c/DS_S9_L001_R1_001.fastq.gz \
DS_ds.d7b863fa32d147a5b0b36434b8d3bb5c/DS_S9_L001_R2_001.fastq.gz


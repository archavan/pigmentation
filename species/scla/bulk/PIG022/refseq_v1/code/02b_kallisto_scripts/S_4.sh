#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_S_4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/S_4

cd /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/01_fastq

# run kallisto
kallisto quant \
-i /home/arc78/genomes/species/scla/refseq_v1/results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx \
-t 25 -b 100 -o /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/S_4 \
ES_ds.185e399ee66844a99ce3f17ffbfac6ed/ES_S15_L001_R1_001.fastq.gz \
ES_ds.185e399ee66844a99ce3f17ffbfac6ed/ES_S15_L001_R2_001.fastq.gz \
ES_ds.9e9f7b9378894ee08ca820433f98ca30/ES_S15_L001_R1_001.fastq.gz \
ES_ds.9e9f7b9378894ee08ca820433f98ca30/ES_S15_L001_R2_001.fastq.gz \
ES_ds.d0db64fc55be4c4c982d2071e302017e/ES_S15_L001_R1_001.fastq.gz \
ES_ds.d0db64fc55be4c4c982d2071e302017e/ES_S15_L001_R2_001.fastq.gz \
ES_ds.efe8d7ed0f9d4be4bcf8d9e8f469eeff/ES_S15_L001_R1_001.fastq.gz \
ES_ds.efe8d7ed0f9d4be4bcf8d9e8f469eeff/ES_S15_L001_R2_001.fastq.gz


#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_I2_4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I2_4

cd /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/01_fastq

# run kallisto
kallisto quant \
-i /home/arc78/genomes/species/scla/refseq_v1/results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx \
-t 25 -b 100 -o /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I2_4 \
E2_ds.2cd6d077ddea46719ba5a17aeb4721cd/E2_S13_L001_R1_001.fastq.gz \
E2_ds.2cd6d077ddea46719ba5a17aeb4721cd/E2_S13_L001_R2_001.fastq.gz \
E2_ds.6b3fa67e9c41406888700f1103daf7e4/E2_S13_L001_R1_001.fastq.gz \
E2_ds.6b3fa67e9c41406888700f1103daf7e4/E2_S13_L001_R2_001.fastq.gz \
E2_ds.a1259afef5bb451186a7da921cf1bec3/E2_S13_L001_R1_001.fastq.gz \
E2_ds.a1259afef5bb451186a7da921cf1bec3/E2_S13_L001_R2_001.fastq.gz \
E2_ds.df7c8c5d319d4725a2bc852226a7d805/E2_S13_L001_R1_001.fastq.gz \
E2_ds.df7c8c5d319d4725a2bc852226a7d805/E2_S13_L001_R2_001.fastq.gz


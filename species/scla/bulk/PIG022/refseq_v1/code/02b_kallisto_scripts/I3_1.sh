#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_I3_1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I3_1

cd /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/01_fastq

# run kallisto
kallisto quant \
-i /home/arc78/genomes/species/scla/refseq_v1/results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx \
-t 25 -b 100 -o /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I3_1 \
A3_ds.04e7af81120442c5a5989730004fe703/A3_S31_L001_R1_001.fastq.gz \
A3_ds.04e7af81120442c5a5989730004fe703/A3_S31_L001_R2_001.fastq.gz \
A3_ds.4a788bdc7b264f07af3870c2fb8b4e92/A3_S31_L001_R1_001.fastq.gz \
A3_ds.4a788bdc7b264f07af3870c2fb8b4e92/A3_S31_L001_R2_001.fastq.gz \
A3_ds.5a70a89bde5a4c0a841a50c336405f6f/A3_S31_L001_R1_001.fastq.gz \
A3_ds.5a70a89bde5a4c0a841a50c336405f6f/A3_S31_L001_R2_001.fastq.gz \
A3_ds.f6e639d7026b4758a534a5e84d3948e8/A3_S31_L001_R1_001.fastq.gz \
A3_ds.f6e639d7026b4758a534a5e84d3948e8/A3_S31_L001_R2_001.fastq.gz


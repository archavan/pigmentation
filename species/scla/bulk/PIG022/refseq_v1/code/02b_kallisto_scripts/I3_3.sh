#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_I3_3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I3_3

cd /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/01_fastq

# run kallisto
kallisto quant \
-i /home/arc78/genomes/species/scla/refseq_v1/results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx \
-t 25 -b 100 -o /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I3_3 \
D3_ds.0f1602db59224e5eb16aff4cacb3cdda/D3_S19_L001_R1_001.fastq.gz \
D3_ds.0f1602db59224e5eb16aff4cacb3cdda/D3_S19_L001_R2_001.fastq.gz \
D3_ds.328258d691f048e29d63c6f2cd3a2ae1/D3_S19_L001_R1_001.fastq.gz \
D3_ds.328258d691f048e29d63c6f2cd3a2ae1/D3_S19_L001_R2_001.fastq.gz \
D3_ds.f069b8f2eebe4a8eb189003a20f4b4b1/D3_S19_L001_R1_001.fastq.gz \
D3_ds.f069b8f2eebe4a8eb189003a20f4b4b1/D3_S19_L001_R2_001.fastq.gz \
D3_ds.fbfdee6fe9de4dd3846ac59339ab32b8/D3_S19_L001_R1_001.fastq.gz \
D3_ds.fbfdee6fe9de4dd3846ac59339ab32b8/D3_S19_L001_R2_001.fastq.gz


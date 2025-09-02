#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_I5_1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I5_1

cd /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/01_fastq

# run kallisto
kallisto quant \
-i /home/arc78/genomes/species/scla/refseq_v1/results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx \
-t 25 -b 100 -o /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I5_1 \
A5_ds.7258dc3598904959b2a3fe16b0441bf7/A5_S35_L001_R1_001.fastq.gz \
A5_ds.7258dc3598904959b2a3fe16b0441bf7/A5_S35_L001_R2_001.fastq.gz \
A5_ds.9f054169ad4940e19b142f0ff52d60ae/A5_S35_L001_R1_001.fastq.gz \
A5_ds.9f054169ad4940e19b142f0ff52d60ae/A5_S35_L001_R2_001.fastq.gz \
A5_ds.eb543fdfeabd4e1fb1207359cf9a0ef9/A5_S35_L001_R1_001.fastq.gz \
A5_ds.eb543fdfeabd4e1fb1207359cf9a0ef9/A5_S35_L001_R2_001.fastq.gz \
A5_ds.f5ebd7f7e3ea4b87ab66b15a2e4deed0/A5_S35_L001_R1_001.fastq.gz \
A5_ds.f5ebd7f7e3ea4b87ab66b15a2e4deed0/A5_S35_L001_R2_001.fastq.gz


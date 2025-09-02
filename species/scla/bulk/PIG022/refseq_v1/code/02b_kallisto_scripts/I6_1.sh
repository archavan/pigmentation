#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_I6_1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I6_1

cd /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/01_fastq

# run kallisto
kallisto quant \
-i /home/arc78/genomes/species/scla/refseq_v1/results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx \
-t 25 -b 100 -o /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I6_1 \
A6_ds.0184a01ea87041a29c5d7251c79c39d0/A6_S37_L001_R1_001.fastq.gz \
A6_ds.0184a01ea87041a29c5d7251c79c39d0/A6_S37_L001_R2_001.fastq.gz \
A6_ds.383f5e73ba6745ff9f68dea0ceab1da6/A6_S37_L001_R1_001.fastq.gz \
A6_ds.383f5e73ba6745ff9f68dea0ceab1da6/A6_S37_L001_R2_001.fastq.gz \
A6_ds.50795f9a85a6415d8c7b468776720632/A6_S37_L001_R1_001.fastq.gz \
A6_ds.50795f9a85a6415d8c7b468776720632/A6_S37_L001_R2_001.fastq.gz \
A6_ds.a142f31132174496ad7ddca87a402ae7/A6_S37_L001_R1_001.fastq.gz \
A6_ds.a142f31132174496ad7ddca87a402ae7/A6_S37_L001_R2_001.fastq.gz


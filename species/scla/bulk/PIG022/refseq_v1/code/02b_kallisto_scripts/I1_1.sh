#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_I1_1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I1_1

cd /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/01_fastq

# run kallisto
kallisto quant \
-i /home/arc78/genomes/species/scla/refseq_v1/results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx \
-t 25 -b 100 -o /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I1_1 \
A1_ds.70f5c8ab61e44b25a12dfcdafcc78fd8/A1_S27_L001_R1_001.fastq.gz \
A1_ds.70f5c8ab61e44b25a12dfcdafcc78fd8/A1_S27_L001_R2_001.fastq.gz \
A1_ds.b3af62b93de34d458db7a8bd909ccc5c/A1_S27_L001_R1_001.fastq.gz \
A1_ds.b3af62b93de34d458db7a8bd909ccc5c/A1_S27_L001_R2_001.fastq.gz \
A1_ds.c960888370484715bac1979be772b8b0/A1_S27_L001_R1_001.fastq.gz \
A1_ds.c960888370484715bac1979be772b8b0/A1_S27_L001_R2_001.fastq.gz \
A1_ds.e487edf9a1b6449b8ba6c0410b665600/A1_S27_L001_R1_001.fastq.gz \
A1_ds.e487edf9a1b6449b8ba6c0410b665600/A1_S27_L001_R2_001.fastq.gz


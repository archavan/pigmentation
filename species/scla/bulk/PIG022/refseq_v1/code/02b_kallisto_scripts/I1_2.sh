#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_I1_2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I1_2

cd /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/01_fastq

# run kallisto
kallisto quant \
-i /home/arc78/genomes/species/scla/refseq_v1/results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx \
-t 25 -b 100 -o /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I1_2 \
B1_ds.4b6ed60d729d43f39ffcecaa28e0b851/B1_S41_L001_R1_001.fastq.gz \
B1_ds.4b6ed60d729d43f39ffcecaa28e0b851/B1_S41_L001_R2_001.fastq.gz \
B1_ds.a0193aeef71f4f6589d1b887a6a449c3/B1_S41_L001_R1_001.fastq.gz \
B1_ds.a0193aeef71f4f6589d1b887a6a449c3/B1_S41_L001_R2_001.fastq.gz \
B1_ds.a646ead2136042479fd54a9591c1588c/B1_S41_L001_R1_001.fastq.gz \
B1_ds.a646ead2136042479fd54a9591c1588c/B1_S41_L001_R2_001.fastq.gz \
B1_ds.eb25e8a8c666467a98308f301f10f493/B1_S41_L001_R1_001.fastq.gz \
B1_ds.eb25e8a8c666467a98308f301f10f493/B1_S41_L001_R2_001.fastq.gz


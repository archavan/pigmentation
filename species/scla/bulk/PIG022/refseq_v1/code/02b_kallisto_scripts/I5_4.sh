#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_I5_4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I5_4

cd /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/01_fastq

# run kallisto
kallisto quant \
-i /home/arc78/genomes/species/scla/refseq_v1/results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx \
-t 25 -b 100 -o /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I5_4 \
E5_ds.3f9eb4dccbbe42f4945544f06ccf0165/E5_S11_L001_R1_001.fastq.gz \
E5_ds.3f9eb4dccbbe42f4945544f06ccf0165/E5_S11_L001_R2_001.fastq.gz \
E5_ds.43eff27635124e8fabd9467c952d5613/E5_S11_L001_R1_001.fastq.gz \
E5_ds.43eff27635124e8fabd9467c952d5613/E5_S11_L001_R2_001.fastq.gz \
E5_ds.aac716be974f4b7786aaa80989e4d397/E5_S11_L001_R1_001.fastq.gz \
E5_ds.aac716be974f4b7786aaa80989e4d397/E5_S11_L001_R2_001.fastq.gz \
E5_ds.fcc14a1ff5144c43b445dc8c2315940d/E5_S11_L001_R1_001.fastq.gz \
E5_ds.fcc14a1ff5144c43b445dc8c2315940d/E5_S11_L001_R2_001.fastq.gz


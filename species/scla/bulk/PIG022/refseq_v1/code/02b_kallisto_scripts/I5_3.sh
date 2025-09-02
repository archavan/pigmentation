#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_I5_3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I5_3

cd /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/01_fastq

# run kallisto
kallisto quant \
-i /home/arc78/genomes/species/scla/refseq_v1/results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx \
-t 25 -b 100 -o /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I5_3 \
D5_ds.401ad5e0e98c434abe9b7252f651acdf/D5_S21_L001_R1_001.fastq.gz \
D5_ds.401ad5e0e98c434abe9b7252f651acdf/D5_S21_L001_R2_001.fastq.gz \
D5_ds.b0de079491d24e8f9fc59e00d94be1d0/D5_S21_L001_R1_001.fastq.gz \
D5_ds.b0de079491d24e8f9fc59e00d94be1d0/D5_S21_L001_R2_001.fastq.gz \
D5_ds.c3cc51ef26be45e49a3a77e845eb4382/D5_S21_L001_R1_001.fastq.gz \
D5_ds.c3cc51ef26be45e49a3a77e845eb4382/D5_S21_L001_R2_001.fastq.gz \
D5_ds.dd0cb1159f8a47289593db8efacfa0fe/D5_S21_L001_R1_001.fastq.gz \
D5_ds.dd0cb1159f8a47289593db8efacfa0fe/D5_S21_L001_R2_001.fastq.gz


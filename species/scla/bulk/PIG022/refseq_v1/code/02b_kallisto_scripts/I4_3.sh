#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_I4_3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I4_3

cd /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/01_fastq

# run kallisto
kallisto quant \
-i /home/arc78/genomes/species/scla/refseq_v1/results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx \
-t 25 -b 100 -o /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I4_3 \
D4_ds.738e454c93cd49689eb4c9fa83075cc3/D4_S20_L001_R1_001.fastq.gz \
D4_ds.738e454c93cd49689eb4c9fa83075cc3/D4_S20_L001_R2_001.fastq.gz \
D4_ds.82a78902a5e84c0b943235919d3d585f/D4_S20_L001_R1_001.fastq.gz \
D4_ds.82a78902a5e84c0b943235919d3d585f/D4_S20_L001_R2_001.fastq.gz \
D4_ds.c1aebfeb2aeb47eaa4b88a8e6e235681/D4_S20_L001_R1_001.fastq.gz \
D4_ds.c1aebfeb2aeb47eaa4b88a8e6e235681/D4_S20_L001_R2_001.fastq.gz \
D4_ds.e29fd121b0b1446a941d5047fb24e5b1/D4_S20_L001_R1_001.fastq.gz \
D4_ds.e29fd121b0b1446a941d5047fb24e5b1/D4_S20_L001_R2_001.fastq.gz


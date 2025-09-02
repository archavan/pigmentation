#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_I3_4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I3_4

cd /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/01_fastq

# run kallisto
kallisto quant \
-i /home/arc78/genomes/species/scla/refseq_v1/results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx \
-t 25 -b 100 -o /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I3_4 \
E3_ds.221302e276374d6f991b048ad7a325e3/E3_S12_L001_R1_001.fastq.gz \
E3_ds.221302e276374d6f991b048ad7a325e3/E3_S12_L001_R2_001.fastq.gz \
E3_ds.2e86cab742314788913c22875653561f/E3_S12_L001_R1_001.fastq.gz \
E3_ds.2e86cab742314788913c22875653561f/E3_S12_L001_R2_001.fastq.gz \
E3_ds.b251d31bfc2c45569a9bc68b56d382c5/E3_S12_L001_R1_001.fastq.gz \
E3_ds.b251d31bfc2c45569a9bc68b56d382c5/E3_S12_L001_R2_001.fastq.gz \
E3_ds.f0f606c9547e4571a8b580f66be7fd85/E3_S12_L001_R1_001.fastq.gz \
E3_ds.f0f606c9547e4571a8b580f66be7fd85/E3_S12_L001_R2_001.fastq.gz


#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_S_2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/S_2

cd /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/01_fastq

# run kallisto
kallisto quant \
-i /home/arc78/genomes/species/scla/refseq_v1/results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx \
-t 25 -b 100 -o /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/S_2 \
BS_ds.2c087f26e1ab4b5dbe745a6720f689ed/BS_S39_L001_R1_001.fastq.gz \
BS_ds.2c087f26e1ab4b5dbe745a6720f689ed/BS_S39_L001_R2_001.fastq.gz \
BS_ds.47a7c788056b4a6ab88827963fbb0579/BS_S39_L001_R1_001.fastq.gz \
BS_ds.47a7c788056b4a6ab88827963fbb0579/BS_S39_L001_R2_001.fastq.gz \
BS_ds.7f09d693053149f58a0a8861a9d6155d/BS_S39_L001_R1_001.fastq.gz \
BS_ds.7f09d693053149f58a0a8861a9d6155d/BS_S39_L001_R2_001.fastq.gz \
BS_ds.bd465fa5baf84cceb4e7f785bbd9a368/BS_S39_L001_R1_001.fastq.gz \
BS_ds.bd465fa5baf84cceb4e7f785bbd9a368/BS_S39_L001_R2_001.fastq.gz


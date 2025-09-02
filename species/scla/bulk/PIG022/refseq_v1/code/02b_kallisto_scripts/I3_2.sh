#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_I3_2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I3_2

cd /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/01_fastq

# run kallisto
kallisto quant \
-i /home/arc78/genomes/species/scla/refseq_v1/results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx \
-t 25 -b 100 -o /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I3_2 \
B3_ds.0e7cc458ff2249f28fb95b26f74a9432/B3_S45_L001_R1_001.fastq.gz \
B3_ds.0e7cc458ff2249f28fb95b26f74a9432/B3_S45_L001_R2_001.fastq.gz \
B3_ds.1ed9f1708a01479ebf5aa31dbfcdef8f/B3_S45_L001_R1_001.fastq.gz \
B3_ds.1ed9f1708a01479ebf5aa31dbfcdef8f/B3_S45_L001_R2_001.fastq.gz \
B3_ds.72d3074e48a043fa8bffd1341a956e4f/B3_S45_L001_R1_001.fastq.gz \
B3_ds.72d3074e48a043fa8bffd1341a956e4f/B3_S45_L001_R2_001.fastq.gz \
B3_ds.c279245c86c94282a934e6f18912da1c/B3_S45_L001_R1_001.fastq.gz \
B3_ds.c279245c86c94282a934e6f18912da1c/B3_S45_L001_R2_001.fastq.gz


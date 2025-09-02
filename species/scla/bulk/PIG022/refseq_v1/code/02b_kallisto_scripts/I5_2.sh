#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_I5_2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I5_2

cd /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/01_fastq

# run kallisto
kallisto quant \
-i /home/arc78/genomes/species/scla/refseq_v1/results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx \
-t 25 -b 100 -o /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I5_2 \
B5_ds.0ba5d7ca61b9499e956b8a0b6acdf9df/B5_S17_L001_R1_001.fastq.gz \
B5_ds.0ba5d7ca61b9499e956b8a0b6acdf9df/B5_S17_L001_R2_001.fastq.gz \
B5_ds.3d98cb8d3c6b4a8381a1ae9b7ed4cb26/B5_S17_L001_R1_001.fastq.gz \
B5_ds.3d98cb8d3c6b4a8381a1ae9b7ed4cb26/B5_S17_L001_R2_001.fastq.gz \
B5_ds.5cb451e425ee4404b811ab1c0c8f4436/B5_S17_L001_R1_001.fastq.gz \
B5_ds.5cb451e425ee4404b811ab1c0c8f4436/B5_S17_L001_R2_001.fastq.gz \
B5_ds.6605166ac3984177accaf01ba95278fa/B5_S17_L001_R1_001.fastq.gz \
B5_ds.6605166ac3984177accaf01ba95278fa/B5_S17_L001_R2_001.fastq.gz


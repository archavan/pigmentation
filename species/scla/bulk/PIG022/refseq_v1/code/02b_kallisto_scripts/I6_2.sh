#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_I6_2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I6_2

cd /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/01_fastq

# run kallisto
kallisto quant \
-i /home/arc78/genomes/species/scla/refseq_v1/results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx \
-t 25 -b 100 -o /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I6_2 \
B6_ds.4edee964131b4e7684b7f09c1f4f6759/B6_S16_L001_R1_001.fastq.gz \
B6_ds.4edee964131b4e7684b7f09c1f4f6759/B6_S16_L001_R2_001.fastq.gz \
B6_ds.5710bb5f4c2f471e99338b95ef93db94/B6_S16_L001_R1_001.fastq.gz \
B6_ds.5710bb5f4c2f471e99338b95ef93db94/B6_S16_L001_R2_001.fastq.gz \
B6_ds.733668c54ea444f8822148d9d7f20808/B6_S16_L001_R1_001.fastq.gz \
B6_ds.733668c54ea444f8822148d9d7f20808/B6_S16_L001_R2_001.fastq.gz \
B6_ds.94c67fb14465498b815c6cebcf1bedb5/B6_S16_L001_R1_001.fastq.gz \
B6_ds.94c67fb14465498b815c6cebcf1bedb5/B6_S16_L001_R2_001.fastq.gz


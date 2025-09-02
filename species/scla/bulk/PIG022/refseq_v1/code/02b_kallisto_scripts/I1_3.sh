#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_I1_3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I1_3

cd /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/01_fastq

# run kallisto
kallisto quant \
-i /home/arc78/genomes/species/scla/refseq_v1/results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx \
-t 25 -b 100 -o /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I1_3 \
D1_ds.0193fce0d3b14580a75890ea682ca95b/D1_S8_L001_R1_001.fastq.gz \
D1_ds.0193fce0d3b14580a75890ea682ca95b/D1_S8_L001_R2_001.fastq.gz \
D1_ds.97cb3b79b17448eab84b0890883306ba/D1_S8_L001_R1_001.fastq.gz \
D1_ds.97cb3b79b17448eab84b0890883306ba/D1_S8_L001_R2_001.fastq.gz \
D1_ds.da64306ef5e1491a83b76164d38a093e/D1_S8_L001_R1_001.fastq.gz \
D1_ds.da64306ef5e1491a83b76164d38a093e/D1_S8_L001_R2_001.fastq.gz \
D1_ds.df21984a0db44e2b9b4098ae6191c555/D1_S8_L001_R1_001.fastq.gz \
D1_ds.df21984a0db44e2b9b4098ae6191c555/D1_S8_L001_R2_001.fastq.gz


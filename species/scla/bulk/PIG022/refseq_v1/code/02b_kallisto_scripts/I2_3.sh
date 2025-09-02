#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_I2_3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I2_3

cd /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/01_fastq

# run kallisto
kallisto quant \
-i /home/arc78/genomes/species/scla/refseq_v1/results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx \
-t 25 -b 100 -o /home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto/I2_3 \
D2_ds.40f0669888f84e858e2e2c7aa6aefdb0/D2_S18_L001_R1_001.fastq.gz \
D2_ds.40f0669888f84e858e2e2c7aa6aefdb0/D2_S18_L001_R2_001.fastq.gz \
D2_ds.68a42d3415f247d985bd92636ffef0fb/D2_S18_L001_R1_001.fastq.gz \
D2_ds.68a42d3415f247d985bd92636ffef0fb/D2_S18_L001_R2_001.fastq.gz \
D2_ds.841195eab0d146d9a9a75ed800685dce/D2_S18_L001_R1_001.fastq.gz \
D2_ds.841195eab0d146d9a9a75ed800685dce/D2_S18_L001_R2_001.fastq.gz \
D2_ds.b638bae1314145d2bb4ce8ae2359c8f6/D2_S18_L001_R1_001.fastq.gz \
D2_ds.b638bae1314145d2bb4ce8ae2359c8f6/D2_S18_L001_R2_001.fastq.gz


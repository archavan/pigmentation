#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_cr24_gut.j
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6 
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

# load modules
module load Python
module load kallisto

# make a new directory for the sample 
mkdir -p /home/arc78/endoderm/species/crob/bulk/ARCA123/KY21/data/02_kallisto/cr24_gut.j

# cd into the data directory
cd /home/arc78/endoderm/species/crob/bulk/ARCA123/KY21/data/01_fastq

# run kallisto
kallisto quant \
-i ~/endoderm/species/crob/genome/KY21/results/kallisto/KY21_transcripts_with-mito.idx \
-t 6 -b 100 -o /home/arc78/endoderm/species/crob/bulk/ARCA123/KY21/data/02_kallisto/cr24_gut.j \
000012213/Sample_cr24_gut_j_381_196/cr24_gut_j_381_196_S43_L001_R1_001.fastq.gz \
000012213/Sample_cr24_gut_j_381_196/cr24_gut_j_381_196_S43_L001_R2_001.fastq.gz


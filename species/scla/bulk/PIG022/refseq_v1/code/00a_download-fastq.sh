#!/usr/bin/env bash

# run interactively on the cluster

cd /home/arc78/palmer_scratch/data/bulk-rna/PIG022/01_fastq

salloc --ntasks 1 --cpus-per-task 1 --mem-per-cpu 25G --time 4:00:00

module load BaseSpaceCLI

bs list biosamples

# | AS            | 775223100 |      
# | A1            | 775223102 |      
# | A2            | 775223104 |      
# | A3            | 775223106 |      
# | A4            | 775223108 |      
# | A5            | 775223110 |      
# | A6            | 775223112 |      
# | BS            | 775223114 |      
# | B1            | 775223116 |      
# | B2            | 775223118 |      
# | B3            | 775223120 |      
# | B4            | 775223122 |      
# | B5            | 775223092 |      
# | B6            | 775223091 |      
# | DS            | 775223084 |      
# | D1            | 775223083 |      
# | D2            | 775223093 |      
# | D3            | 775223094 |      
# | D4            | 775223095 |      
# | D5            | 775223096 |      
# | D6            | 775223097 |      
# | ES            | 775223090 |      
# | E1            | 775223089 |      
# | E2            | 775223088 |      
# | E3            | 775223087 |      
# | E4            | 775223098 |      
# | E5            | 775223086 |      
# | E6            | 775223085 |      

bs download biosample --name AS --extension fastq.gz      
bs download biosample --name A1 --extension fastq.gz      
bs download biosample --name A2 --extension fastq.gz      
bs download biosample --name A3 --extension fastq.gz      
bs download biosample --name A4 --extension fastq.gz      
bs download biosample --name A5 --extension fastq.gz      
bs download biosample --name A6 --extension fastq.gz      
bs download biosample --name BS --extension fastq.gz      
bs download biosample --name B1 --extension fastq.gz      
bs download biosample --name B2 --extension fastq.gz      
bs download biosample --name B3 --extension fastq.gz      
bs download biosample --name B4 --extension fastq.gz      
bs download biosample --name B5 --extension fastq.gz      
bs download biosample --name B6 --extension fastq.gz      
bs download biosample --name DS --extension fastq.gz      
bs download biosample --name D1 --extension fastq.gz      
bs download biosample --name D2 --extension fastq.gz      
bs download biosample --name D3 --extension fastq.gz      
bs download biosample --name D4 --extension fastq.gz      
bs download biosample --name D5 --extension fastq.gz      
bs download biosample --name D6 --extension fastq.gz      
bs download biosample --name ES --extension fastq.gz      
bs download biosample --name E1 --extension fastq.gz      
bs download biosample --name E2 --extension fastq.gz      
bs download biosample --name E3 --extension fastq.gz      
bs download biosample --name E4 --extension fastq.gz      
bs download biosample --name E5 --extension fastq.gz      
bs download biosample --name E6 --extension fastq.gz    

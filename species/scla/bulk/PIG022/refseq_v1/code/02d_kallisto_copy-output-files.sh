#!/usr/bin/env bash

# set download locations ======================================================
cluster=arc78@mccleary.ycrc.yale.edu:/home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/02_kallisto
computer=/Users/arun/Dropbox/research/projects/endoderm/species/scla/bulk/PIG022/refseq_v1/data/

# download ====================================================================

rsync -avz --exclude "*.h5" $cluster $computer
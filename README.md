## Organization of digestive functions in the chordate gut

Code associated with:

Chavan, A. R.\*, Yaghoubi, P.\*, Langner, O. C., Menges, F. S., Gibson, K. H., Jia, K., Musser, J. M., Miller, S. J., Dishaw, L. J., & Medzhitov, R. (2026). *Organization of digestive functions in the chordate gut* [Preprint]. *bioRxiv*. https://doi.org/10.64898/2026.11.24.768686  
(\*equal contribution)

### Code

Code for figures in the manuscript is under the directory `figures`. Code for data processing and exploratory analysis is under the directory for each species under `species`.

### RNA-seq data

Associated bulk RNA-sequencing data is available in NCBI GEO under the SuperSeries [GSE354417](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE354417) (*Ciona robusta*: [GSE354357](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE354357); *Styela clava*: [GSE354416](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE354416)). For both species, the final cleaned counts and tpm files used for analysis can be downloaded from the processed data section of the GEO record.

### Genomic data


Code for preparation of genomic data and processed lists of functional module genes are in the [genomes](https://github.com/archavan/genomes) repository release [pigmentation_v1](https://github.com/archavan/genomes/releases/tag/pigmentation_v1) (commit [ddbf124](https://github.com/archavan/genomes/commit/ddbf1240fe16d055ee2c8691f5ddd868e7c82e0a)). To retrieve the version of files from this release, as they were used for this manuscript, you can use `git clone --branch pigmentation_v1 https://github.com/archavan/genomes.git`. The lists of functional module genes can also be directly downloaded from the links below.

- *Mus musculus*: https://github.com/archavan/genomes/blob/pigmentation_v1/species/mmus/gene-modules/functional-modules_gut.tsv
- *Ciona robusta*: https://github.com/archavan/genomes/blob/pigmentation_v1/species/crob/KY21/results/genes/functional-modules_gut.tsv
- *Styela clava*: https://github.com/archavan/genomes/blob/pigmentation_v1/species/scla/refseq_v1/results/genes/functional-modules_gut.tsv

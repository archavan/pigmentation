# setup =======================================================================
library(tidyverse)
library(here)
library(glue)
library(fs)

basedir <- here("species/scla/bulk/PIG022/refseq_v1")

md <- read_csv(here(basedir, "data/samples.csv"))

# functions ===================================================================

generate_kallisto_script <- function(
  sample_id,
  fastq.files,
  index = "/home/arc78/genomes/species/scla/refseq_v1/results/kallisto/GCF_013122585.1_ASM1312258v2_rna.idx",
  basedir = "/home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1"
) {
  script <- glue(
    "
#!/usr/bin/env bash
#SBATCH --partition=pi_medzhitov
#SBATCH --job-name=kallisto_{sample_id}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25 
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arun.chavan@yale.edu

module load kallisto

mkdir -p {basedir}/data/02_kallisto/{sample_id}

cd {basedir}/data/01_fastq

# run kallisto
kallisto quant \\\\\n
-i {index} \\\\\n
-t 25 -b 100 -o {basedir}/data/02_kallisto/{sample_id} \\\\\n
{fastq.files}

  "
  )

  return(script)
}

# generate fastq filenames ====================================================
# read filenames
files <- fs::dir_info(
  here(basedir, "data/01_fastq"),
  recurse = TRUE,
  type = "file",
  glob = "*.fastq.gz"
) %>%
  select(path) %>%
  mutate(path = as.character(path)) %>%
  mutate(path = gsub(
    "/home/arc78/endoderm/species/scla/bulk/PIG022/refseq_v1/data/01_fastq/", 
    "", 
    path)) # get path relative to fastq dir

files <- files %>% 
  filter(!str_detect(path, "^Undetermined")) 

# extract features
files$filename <- gsub("(.*)(/)(.+.fastq.gz)", "\\3", files$path)
files$seqfile_id <- gsub(
  "(.+?)_(S\\d+_L\\d+_R\\d+_\\d+\\.fastq.gz)",
  "\\1",
  files$filename
)
files$mate <- gsub("(.*)(R[1-2])(.*)", "\\2", files$filename)

# make sample ids consistent
files <- left_join(files, md, by = "seqfile_id")

# list R1 R2 file pairs:
# We want a list (as a character string) of all filenames for a sample, but
# explicitly ensure that R1 and R2 mate files from each pair of files come after
# one another, first we will create strings "pairA_R1.fastq.gz
# pairA_R2.fastq.gz" for each pair.

files$unique_filebase <- gsub("(.*)(R[1-2]_001.fastq.gz)", "\\1", files$path)

file_pairs <- files %>%
  select(unique_filebase) %>%
  unique()

file_pairs$matepairs <- NA_character_

for (i in 1:nrow(file_pairs)) {
  file_pairs$matepairs[i] <- str_c(
    files$path[which(
      files$unique_filebase == file_pairs$unique_filebase[i] &
        files$mate == "R1"
    )],
    files$path[which(
      files$unique_filebase == file_pairs$unique_filebase[i] &
        files$mate == "R2"
    )],
    sep = " \\\n" # to break up long lines in the final bash script
  )
}

# generate final file_list per sample_id by contatenating pairs by sample_id
file_pairs <- file_pairs %>% # add sample_id
  left_join(
    (files %>%
      select(unique_filebase, sample_id) %>%
      distinct()),
    by = "unique_filebase"
  )

file_list <- file_pairs %>%
  group_by(sample_id) %>%
  summarize(file_list = list(matepairs)) %>%
  mutate(
    file_list = map_chr(.x = file_list, .f = str_flatten, collapse = " \\\n")
  )

# test: writeLines(file_list$file_list)

# generate scripts ============================================================
file_list$kallisto <- NA_character_

for (i in 1:nrow(file_list)) {
  file_list$kallisto[i] <- generate_kallisto_script(
    sample_id = file_list$sample_id[i],
    fastq.files = file_list$file_list[i]
  )
}

# writeLines(file_list$kallisto[4])

# write script files ==========================================================
fs::dir_create(here(basedir, "code/02b_kallisto_scripts"))
for (i in 1:nrow(file_list)) {
  writeLines(
    file_list$kallisto[i],
    here(
      basedir,
      glue("code/02b_kallisto_scripts/{file_list$sample_id[i]}.sh")
    )
  )
}

# write run_all script ========================================================
runall <- glue("sbatch {file_list$sample_id}.sh")

runall <- str_flatten(runall, collapse = "\n")

runall <- str_c(
  "#!/usr/bin/env bash",
  "cd ~/endoderm/species/scla/bulk/PIG022/refseq_v1/code/02b_kallisto_scripts",
  runall,
  sep = "\n\n"
)

writeLines(runall, here(basedir, "code/02c_kallisto_run-all-scripts.sh"))

# end  ========================================================================

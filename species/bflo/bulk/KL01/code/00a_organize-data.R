# Arun Chavan
# Started: 2023-12-27

# background ==================================================================

# Here, we organize the data and metadata for amphioxus bulk RNA-seq data
# received from Kun-Lung Li and Sky Yu.

# setup =======================================================================
library(tidyverse)
library(here)

basedir <- here("species/bflo/bulk/KL01/")
genomedir <- here("species/bflo/genome/VNyyK")

# data ========================================================================
cts <- read_tsv(here(basedir, "data/from-KL/Fcount_v231120-re3.txt"))
coldata <- read_csv(here(basedir, "data/from-KL/colData.csv"))

glen <- read_tsv(here(genomedir, "results/genes/gene-length-median.tsv"))

# organize metadata ===========================================================
md <- coldata %>% 
  rename(sample_name_old = `...1`, 
         rep_id = BioRep) %>% 
  mutate(tissue_full = case_when(
    Organ == "En" ~ "endostyle",
    Organ == "Gi" ~ "gills",
    Organ == "Mi" ~ "midgut", # midgut = stomach + ileocolon ring
    Organ == "Lv" ~ "midgut diverticulum",
    Organ == "Hd" ~ "hindgut"
  )) %>% 
  mutate(tissue = case_when(
    tissue_full == "endostyle"           ~ "estl",
    tissue_full == "gills"               ~ "gill",
    tissue_full == "midgut"              ~ "mgut",
    tissue_full == "midgut diverticulum" ~ "mdiv",
    tissue_full == "hindgut"             ~ "hgut"
  )) %>% 
  mutate(sample_name = paste0(tissue, "_", rep_id))

# organize data ===============================================================
# rename samples 
rename_samples <- function(x) {
  new_names <- rep(NA_character_, length(x))
  
  for (i in seq_along(x)) {
    if (x[i] %in% md$sample_name_old) {
      new_names[i] <- md$sample_name[md$sample_name_old == x[i]]
    } else {
      new_names[i] <- x[i]
    }
  }
  
  return(new_names)
}

cts <- cts %>% 
  rename_with(rename_samples) %>% 
  rename(gene_id = Geneid) 

# convert to integer just to be sure
cts <- cts %>% 
  mutate(across(where(is.double), as.integer))

# calculate tpm ===============================================================
calculate_tpm <- function(feature_counts, feature_lengths) {
  ratio <- feature_counts/feature_lengths
  ratio_sum = sum(ratio)
  tpm <- (ratio/ratio_sum) * 10^6
  return(tpm)
}

tpm_long <- cts %>% 
  pivot_longer(-gene_id, names_to = "sample_name", values_to = "counts") %>% 
  left_join(glen, by = "gene_id") %>% 
  nest(.by = sample_name) %>% 
  mutate(data = map(data, 
                    ~ .x %>% 
                      mutate(tpm = calculate_tpm(counts, length)))) %>% 
  unnest(data)

tpm <- tpm_long %>% 
  pivot_wider(id_cols = "gene_id", 
              names_from = "sample_name", 
              values_from = "tpm")

tpm_mean <- tpm_long %>% 
  separate(sample_name, c("tissue", "rep_id")) %>% 
  group_by(tissue, gene_id) %>% 
  summarize(tpm = mean(tpm)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = "gene_id", names_from = "tissue", values_from = "tpm")

# select only gut samples =====================================================
select_gut_samples <- function(x) {
  select(x, gene_id, starts_with("mgut"), starts_with("mdiv"), starts_with("hgut"))
}

cts_gut <- select_gut_samples(cts)
tpm_gut <- select_gut_samples(tpm)
tpm_mean_gut <- select_gut_samples(tpm_mean)

md_gut <- md %>% 
  filter(tissue %in% c("mgut", "mdiv", "hgut"))

# write files =================================================================
write_tsv(md, here(basedir, "data/00_metadata/samples.tsv"))
write_tsv(md_gut, here(basedir, "data/00_metadata/samples_gut.tsv"))

write_tsv(cts, here(basedir, "data/03_compiled/counts.tsv"))
write_tsv(cts_gut, here(basedir, "data/03_compiled/counts_gut.tsv"))

write_tsv(tpm, here(basedir, "data/03_compiled/tpm.tsv"))
write_tsv(tpm_gut, here(basedir, "data/03_compiled/tpm_gut.tsv"))

write_tsv(tpm_mean, here(basedir, "data/03_compiled/tpm_mean.tsv"))
write_tsv(tpm_mean_gut, here(basedir, "data/03_compiled/tpm_mean_gut.tsv"))

# end =========================================================================

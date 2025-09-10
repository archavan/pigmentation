# Arun Chavan
# Started: 2025-09-09

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)
library(glue)
library(conflicted)
library(ggrepel)

conflicts_prefer(
  dplyr::filter,
  dplyr::select,
  purrr::map
)

scdir <- here("species/scla/bulk/PIG022/refseq_v1")
crdir <- here("species/crob/bulk/ARCA123/KY21")

scgenomedir <- here("species/scla/genome/refseq_v1")
crgenomedir <- here("species/crob/genome/KY21")

orthdir <- here("~/Dropbox/research/projects/genomes/orthofinder/2025-08-29")

resdir <- here(scdir, "results/compare/crob/all-samples/hogs")
fs::dir_create(resdir)

source(here(scdir, "code/utilities_ggplot.R"))

# data ========================================================================
md_sc <- read_csv(here(scdir, "data/samples_clean.csv"))
tpm_sc <- read_csv(here(scdir, "data/03_compiled/tpm_clean.csv"))

md_cr <- read_csv(here(crdir, "data/00_metadata/samples.csv"))
tpm_cr <- read_csv(here(crdir, "data/03_compiled/tpm_no-testis.csv"))

# genomic data ================================================================
genes_sc <- read_tsv(here(scgenomedir, "results/genes/gene-display-names.tsv"))
genes_cr <- read_tsv(here(crgenomedir, "results/genes/gene-display-names.tsv"))

sc_cr_ologs <- read_tsv(here(
  orthdir,
  "results/pairwise-orthologs/scla_vs_crob.tsv"
))

ologs <- sc_cr_ologs |>
  rename(gene_id_sc = scla, gene_id_cr = crob)

###############################################################################
# prepare data ================================================================
###############################################################################

## rename crob samples to match scla naming style =============================
md_cr <- md_cr |>
  mutate(
    sample = case_when(
      tissue == "int.1.w" ~ "I1",
      tissue == "int.2.y" ~ "I2",
      tissue == "int.3.w" ~ "I3",
      tissue == "int.4.y" ~ "I4",
      tissue == "stm.1.s" ~ "S"
    )
  ) |>
  rename(rep = rep_id) |>
  mutate(sample_id = paste0(sample, "_", rep))

rename_crob_samples <- function(oldname) {
  md_cr$sample_id[which(md_cr$sample_name == oldname)]
}

tpm_cr <- tpm_cr |>
  rename_with(rename_crob_samples, .cols = where(is.numeric))

## metadata ===================================================================
md <- bind_rows(
  (md_cr |>
    mutate(species = "crob") |>
    select(species, sample_id, sample_id, sample, rep)),
  (md_sc |>
    mutate(species = "scla") |>
    select(species, sample_id, sample, rep))
) |>
  mutate(
    spp_sample_id = paste0(substr(species, 1, 2), "__", sample_id),
    spp_sample = paste0(substr(species, 1, 2), "__", sample)
  )

## define sample order ========================================================
sample_order_cr <- c("S", "I1", "I2", "I3", "I4")
sample_order_sc <- c("S", "I1", "I2", "I3", "I4", "I5", "I6")
sample_order <- c(
  paste0("cr__", sample_order_cr),
  paste0("sc__", sample_order_sc)
)
md <- mutate(md, spp_sample = fct(spp_sample, sample_order))

sample_id_order <- arrange(md, spp_sample) |> pull(spp_sample_id)


## tpm ========================================================================

tpm_sc_ologs <- inner_join(
  x = ologs |>
    select(hog, gene_id_sc) |>
    distinct(),
  y = tpm_sc,
  by = c("gene_id_sc" = "gene_id")
) |>
  rename_with(
    function(x) paste0("sc__", x),
    where(is.numeric)
  ) |>
  select(-gene_id_sc) |>
  pivot_longer(
    cols = where(is.numeric),
    names_to = "spp_sample_id",
    values_to = "tpm"
  ) |>
  group_by(spp_sample_id, hog) |>
  summarize(tpm = sum(tpm)) |>
  ungroup() |>
  pivot_wider(id_cols = hog, names_from = spp_sample_id, values_from = tpm)

tpm_cr_ologs <- inner_join(
  x = ologs |>
    select(hog, gene_id_cr) |>
    distinct(),
  y = tpm_cr,
  by = c("gene_id_cr" = "gene_id")
) |>
  rename_with(
    function(x) paste0("cr__", x),
    where(is.numeric)
  ) |>
  select(-gene_id_cr) |>
  pivot_longer(
    cols = where(is.numeric),
    names_to = "spp_sample_id",
    values_to = "tpm"
  ) |>
  group_by(spp_sample_id, hog) |>
  summarize(tpm = sum(tpm)) |>
  ungroup() |>
  pivot_wider(id_cols = hog, names_from = spp_sample_id, values_from = tpm)


tpm <- inner_join(x = tpm_sc_ologs, y = tpm_cr_ologs, by = "hog")


# scaling =====================================================================
barescale <- function(x, ...) scale(x, ...)[, 1]

expressed_in_both <- tpm |>
  pivot_longer(
    cols = where(is.numeric),
    names_to = "spp_sample_id",
    values_to = "tpm"
  ) |>
  left_join(select(md, spp_sample_id, species)) |>
  group_by(species, hog) |>
  summarize(tpm_sum = sum(tpm)) |>
  ungroup() |>
  pivot_wider(
    id_cols = "hog",
    names_from = "species",
    values_from = "tpm_sum"
  ) |>
  filter(crob > 0 & scla > 0) |>
  pull(hog)


ztpm_wss <- tpm |>
  filter(hog %in% expressed_in_both) |>
  pivot_longer(
    cols = where(is.numeric),
    names_to = "spp_sample_id",
    values_to = "tpm"
  ) |>
  left_join(select(md, spp_sample_id, species), by = "spp_sample_id") |>
  group_by(species, hog) |>
  mutate(ztpm = barescale(tpm)) |>
  ungroup() |>
  pivot_wider(id_cols = hog, values_from = ztpm, names_from = spp_sample_id)

ztpm_wss_mean <- ztpm_wss |>
  pivot_longer(
    cols = where(is.numeric),
    values_to = "ztpm",
    names_to = "spp_sample_id"
  ) |>
  left_join(select(md, spp_sample_id, spp_sample), by = "spp_sample_id") |>
  group_by(spp_sample, hog) |>
  summarize(ztpm = mean(ztpm)) |>
  ungroup() |>
  pivot_wider(id_cols = hog, values_from = ztpm, names_from = spp_sample)

###############################################################################
# correlations ================================================================
###############################################################################

corrs_wss_full_pearson <- ztpm_wss |>
  select(where(is.numeric)) |>
  cor(method = "pearson") |>
  as.data.frame() |>
  rownames_to_column("v1") |>
  pivot_longer(cols = where(is.numeric), names_to = "v2", values_to = "corr") |>
  filter(substr(v1, 1, 2) != substr(v2, 1, 2)) |>
  filter(substr(v1, 1, 2) != "cr") |>
  mutate(v1 = fct(v1, sample_id_order), v2 = fct(v2, sample_id_order)) |>
  arrange(v1, v2) |>
  pivot_wider(names_from = "v2", values_from = "corr") |>
  column_to_rownames("v1") |>
  pheatmap::pheatmap(cluster_rows = FALSE, cluster_cols = FALSE)

pdf(here(resdir, "corrs_wss_full_pearson.pdf"), 5.25, 6)
print(corrs_wss_full_pearson)
dev.off()


corrs_wss_mean_pearson <- ztpm_wss_mean |>
  select(where(is.numeric)) |>
  cor(method = "pearson") |>
  as.data.frame() |>
  rownames_to_column("v1") |>
  pivot_longer(cols = where(is.numeric), names_to = "v2", values_to = "corr") |>
  filter(substr(v1, 1, 2) != substr(v2, 1, 2)) |>
  filter(substr(v1, 1, 2) != "cr") |>
  mutate(v1 = fct(v1, sample_order), v2 = fct(v2, sample_order)) |>
  arrange(v1, v2) |>
  pivot_wider(names_from = "v2", values_from = "corr") |>
  column_to_rownames("v1") |>
  pheatmap::pheatmap(cluster_rows = FALSE, cluster_cols = FALSE)

pdf(here(resdir, "corrs_wss_mean_pearson.pdf"), 3, 3)
print(corrs_wss_mean_pearson)
dev.off()

# end =========================================================================

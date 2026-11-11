# Arun Chavan
# Started: 2026-11-11

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)

basedir <- here("geo/scla")
scdir <- here("species/scla/bulk/PIG022/refseq_v1")

outdir <- here(basedir, "results/processed-data")
fs::dir_create(here(outdir))

###############################################################################
# data ========================================================================
###############################################################################

fq <- read_csv(here(basedir, "results/fq_all_tidy.csv"))
md <- read_csv(here(scdir, "data/samples.csv")) |>
  filter(sample_id %in% fq$sample_id)

df <- tribble(
  ~id   , ~filename          , ~outname     ,
  "cts" , "counts_clean.csv" , "counts.csv" ,
  "tpm" , "tpm_clean.csv"    , "tpm.csv"
) |>
  mutate(obj = map(filename, ~ read_csv(here(scdir, "data/03_compiled", .x))))

###############################################################################
# reorder samples ==============================================================
###############################################################################

sample_order <- c("S", paste0("I", 1:6))

id_order <- md |>
  mutate(sample = fct(sample, sample_order)) |>
  arrange(sample, rep) |>
  pull(sample_id)

reorder_samples <- function(.df) {
  stopifnot(setequal(setdiff(names(.df), "gene_id"), id_order))
  select(.df, gene_id, all_of(id_order))
}

df <- df |>
  mutate(obj = map(obj, reorder_samples))

###############################################################################
# write files =================================================================
###############################################################################

walk2(df$obj, df$outname, ~ write_csv(.x, here(outdir, .y)))

###############################################################################
# get md5 sums ================================================================
###############################################################################

df |>
  mutate(md5 = map_chr(outname, ~ tools::md5sum(here(outdir, .x)))) |>
  select(file = outname, md5) |>
  write_csv(here(basedir, "results/processed-data-md5.csv"))

# end #########################################################################

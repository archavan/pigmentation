# Arun Chavan
# Started: 2026-11-11

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)

basedir <- here("geo/crob")
crdir <- here("species/crob/bulk/ARCA123/KY21")

outdir <- here(basedir, "results/processed-data")
fs::dir_create(here(outdir))

###############################################################################
# data ========================================================================
###############################################################################

fq <- read_csv(here(basedir, "results/fq_all_tidy.csv"))
md <- read_csv(here(crdir, "data/00_metadata/samples.csv")) |>
  filter(tissue_name %in% fq$sample_id)

df <- tribble(
  ~id             , ~filename              , ~outname                            ,
  "cts"           , "counts.csv"           , "counts.csv"                        ,
  "tpm"           , "tpm.csv"              , "tpm.csv"                           ,
  "cts_no_testis" , "counts_no-testis.csv" , "counts_excluding-testis-genes.csv" ,
  "tpm_no_testis" , "tpm_no-testis.csv"    , "tpm_excluding-testis-genes.csv"
) |>
  mutate(obj = map(filename, ~ read_csv(here(crdir, "data/03_compiled", .x))))

###############################################################################
# rename samples ==============================================================
###############################################################################

## rename samples to match their names of GEO metadata

sample_order <- c("S", paste0("I", 1:4))

dict <- md |>
  mutate(tissue_short = fct(tissue_short, sample_order)) |>
  arrange(tissue_short, rep_id) |>
  select(tissue_name, sample_name) |>
  deframe()

rename_samples <- function(.df) {
  stopifnot(setequal(setdiff(names(.df), "gene_id"), dict))
  select(.df, gene_id, all_of(dict))
}

df <- df |>
  mutate(obj = map(obj, rename_samples))

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

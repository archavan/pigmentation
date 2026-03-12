# Arun Chavan
# Started: 2026-03-11

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)
library(glue)

basedir <- here("species/crob/bulk/ARCA123/KY21")

###############################################################################
# data ========================================================================
###############################################################################

tpm <- read_csv(here(basedir, "data/03_compiled/tpm_no-testis.csv"))
samples <- read_csv(here(basedir, "data/00_metadata/samples.csv")) |>
  filter(tissue_short != "S")

tpm <- tpm |>
  select(all_of(c("gene_id", samples$sample_name))) |>
  rename_with(
    function(x) samples$tissue_name[match(x, samples$sample_name)],
    .cols = !starts_with("gene_id")
  )

tpm_long <- tpm |>
  pivot_longer(cols = -gene_id, names_to = "tissue_name", values_to = "tpm") |>
  separate(tissue_name, into = c("tissue_short", "rep_id"), sep = "_") |>
  mutate(tissue_short = fct(tissue_short, c("I1", "I2", "I3", "I4")))

###############################################################################
# L1 normalization / unit sum scaling =========================================
###############################################################################

# Since all tpm values are positive, anyway, L1 normalization in this case is
# identical to unit sum scaling.

l1_normalize <- function(x) {
  x / sum(abs(x))
}

## remove lowly expressed genes ====
tpm_thresh <- 5

# expressed over threhsold in at least one sample in all replicates
expressed_genes <- tpm_long |>
  summarize(max_tpm = max(tpm), .by = c(gene_id, rep_id)) |>
  arrange(desc(max_tpm)) |>
  summarize(all_over_thresh = all(max_tpm > tpm_thresh), .by = gene_id) |>
  filter(all_over_thresh) |>
  pull(gene_id)

tpm_long <- tpm_long |> filter(gene_id %in% expressed_genes)

## normalization ====
uss_long <- tpm_long |>
  group_by(gene_id, rep_id) |>
  mutate(uss = l1_normalize(tpm)) |>
  ungroup()

uss_long |>
  ggplot(aes(x = uss)) +
  stat_ecdf()

stopifnot(
  "USS values do not sum to 1 for all gene_id/rep_id groups" = uss_long |>
    summarise(sum_uss = sum(uss), .by = c(gene_id, rep_id)) |>
    pull(sum_uss) |>
    near(1) |>
    all()
)

###############################################################################
# save data ===================================================================
###############################################################################

uss_wide <- uss_long |>
  pivot_wider(
    id_cols = gene_id,
    names_from = c("tissue_short", "rep_id"),
    names_sep = "_",
    values_from = "uss"
  )

write_csv(uss_wide, here(basedir, "data/03_compiled/uss_intestine.csv"))

# end #########################################################################

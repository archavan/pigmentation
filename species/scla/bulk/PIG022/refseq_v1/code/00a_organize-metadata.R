# Arun Chavan
# Started: 2025-08-18

# background ==================================================================

# setup =======================================================================
library(tidyverse)
library(here)
library(glue)

basedir <- here("species/scla/bulk/PIG022/refseq_v1")

# data ========================================================================
md <- tibble(
  seqfile_id = outer(
    c("A", "B", "D", "E"),
    c("S", "1", "2", "3", "4", "5", "6"),
    paste0
  ) |>
    t() |>
    as.vector()
)

md <- md |>
  separate(seqfile_id, into = c("animal", "sample"), remove = FALSE, sep = 1) |>
  mutate(sample = ifelse(sample %in% 1:6, paste0("I", sample), sample)) |>
  mutate(rep = as.integer(as.factor(animal))) |>
  mutate(sample_name = paste0(sample, "_", animal)) |>
  mutate(sample_id = paste0(sample, "_", rep)) |>
  mutate(exp_id = "PIG022") |>
  relocate(exp_id)

write_csv(md, here(basedir, "data/samples.csv"))

# end =========================================================================

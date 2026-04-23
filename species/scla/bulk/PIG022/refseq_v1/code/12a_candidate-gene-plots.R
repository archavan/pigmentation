# Arun Chavan
# Started: 2026-04-23

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)
library(glue)
library(ggrepel)
library(ggtext)

basedir <- here("species/scla/bulk/PIG022/refseq_v1")
genomedir <- here("species/scla/genome/refseq_v1")
resdir <- here(basedir, "results/12a_candidate-gene-plots")
fs::dir_create(resdir)

source(here(basedir, "code/utilities_ggplot.R"))

###############################################################################
# data ========================================================================
###############################################################################

tpm <- read_csv(here(basedir, "data/03_compiled/tpm_clean.csv"))
uss <- read_csv(here(basedir, "data/03_compiled/uss_intestine.csv"))
samples <- read_csv(here(basedir, "data/samples_clean.csv"))

tpm_long <- tpm |>
  pivot_longer(cols = -gene_id, names_to = "sample_id", values_to = "tpm") |>
  left_join(samples, by = "sample_id") |>
  mutate(sample = fct(sample, c("S", paste0("I", 1:6))))

genes <- read_tsv(here(genomedir, "results/genes/gene-display-names.tsv"))

###############################################################################
# plots =======================================================================
###############################################################################

get_hog_label <- function(.gene_id) {
  genes |>
    filter(gene_id %in% .gene_id) |>
    distinct(gene_id, hog_label) |>
    deframe() |>
    str_replace("(.+)_(.+)", "\\2 (\\1)")
}

plot_genes <- function(.gene_ids) {
  genes |>
    filter(gene_id %in% .gene_ids) |>
    inner_join(tpm_long, by = "gene_id") |>
    ggplot(aes(sample, tpm)) +
    geom_point(color = "grey75") +
    geom_line(
      data = function(x) {
        x |> summarize(tpm = median(tpm), .by = c(gene_id, sample))
      },
      mapping = aes(sample, tpm, group = gene_id)
    ) +
    geom_point(
      data = function(x) {
        x |> summarize(tpm = median(tpm), .by = c(gene_id, sample))
      },
      mapping = aes(sample, tpm)
    ) +
    facet_wrap(
      vars(gene_id),
      scale = "free_y",
      labeller = labeller(gene_id = function(x) {
        glue("{x}<br>{get_hog_label(x)}")
      })
    ) +
    theme_half_open_grid() +
    theme(
      strip.clip = "off",
      strip.text.x = element_markdown()
    )
}

## parahox genes ==============================================================
genes |>
  filter(str_detect(hog_label, "Pdx1") | str_detect(hog_label, "Cdx")) |>
  pull(gene_id) |>
  plot_genes()

ggsave(here(resdir, "parahox-genes.pdf"), width = 6, height = 3)

## hox genes ==================================================================
genes |>
  filter(str_detect(hog_label, regex("Hox", FALSE))) |>
  pull(gene_id) |>
  plot_genes()

ggsave(here(resdir, "hox-genes.pdf"), width = 8, height = 5)

# end #########################################################################

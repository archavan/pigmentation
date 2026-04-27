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
samples <- read_csv(here(basedir, "data/samples_clean.csv"))

tpm_long <- tpm |>
  pivot_longer(cols = -gene_id, names_to = "sample_id", values_to = "tpm") |>
  left_join(samples, by = "sample_id") |>
  mutate(sample = fct(sample, c("S", paste0("I", 1:6))))

genes <- read_tsv(here(genomedir, "results/genes/gene-display-names.tsv"))

###############################################################################
# plots =======================================================================
###############################################################################

plot_genes <- function(.gene_ids) {
  genes |>
    filter(gene_id %in% .gene_ids) |>
    distinct(gene_id, hog_label) |>
    mutate(
      facet_label = glue(
        "{gene_id}<br>{str_replace(hog_label, '(.+)_(.+)', '\\\\2 (\\\\1)')}"
      )
    ) |>
    inner_join(tpm_long, by = "gene_id") |>
    ggplot(aes(sample, tpm)) +
    geom_point(color = "grey75") +
    stat_summary(fun = median, geom = "line", aes(group = facet_label)) +
    stat_summary(fun = median, geom = "point") +
    facet_wrap(vars(facet_label), scales = "free_y") +
    theme_half_open_grid() +
    theme(
      strip.clip = "off",
      strip.text.x = element_markdown()
    )
}

## parahox genes ==============================================================
p_parahox <- genes |>
  filter(str_detect(hog_label, "Pdx1") | str_detect(hog_label, "Cdx")) |>
  pull(gene_id) |>
  plot_genes()

ggsave(here(resdir, "parahox-genes.pdf"), p_parahox, width = 6, height = 3)

## hox genes ==================================================================
p_hox <- genes |>
  filter(str_detect(hog_label, "Hox")) |>
  pull(gene_id) |>
  plot_genes()

ggsave(here(resdir, "hox-genes.pdf"), p_hox, width = 8, height = 5)

## alkaline phosphatase =======================================================

p_alkphos <- genes |>
  filter(str_detect(hog_label, "m4s7_Akp3.Alp\\(i,l,pl2\\)")) |>
  pull(gene_id) |>
  plot_genes()

ggsave(
  here(resdir, "alkaline-phosphatase.pdf"),
  p_alkphos,
  width = 8,
  height = 6.5
)

# end #########################################################################

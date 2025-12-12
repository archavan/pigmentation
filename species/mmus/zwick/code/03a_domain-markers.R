# Arun Chavan
# Started: 2025-12-11

# setup =======================================================================
library(tidyverse)
library(here)
library(glue)
library(ggrepel)
library(Seurat)

basedir <- here("species/mmus/zwick")
resdir <- here(basedir, "results/markers")
fs::dir_create(resdir)

source(here(basedir, "code/utilities_ggplot.R"))
source(here(basedir, "code/utilities_seurat-plots.R"))

plot_dotplot2 <- function(.genes) {
  dp_extract_dotplot_data(ec, .genes) |>
    dp_define_feature_order() |>
    dp_define_cluster_order(LETTERS[1:5]) |>
    dp_split_dotplot_data() |>
    map(dp_plot)
}

# data ========================================================================
ec <- read_rds(here(basedir, "data/enterocytes.rds"))

# find markers ================================================================
ec <- SetIdent(ec, value = "domain")

markers_raw <- FindAllMarkers(ec, only.pos = TRUE)

markers_top <- markers_raw |>
  filter(p_val_adj < 0.05, avg_log2FC > 1) |>
  slice_head(n = 50, by = cluster) |>
  as_tibble() |>
  rename(domain = cluster)

for (i in unique(markers_top$domain)) {
  markers_top |>
    filter(domain == i) |>
    pull(gene) |>
    unique() |>
    plot_dotplot2() |>
    save_dotplot(
      here(basedir, "results/markers", paste0("top-markers_", i, ".pdf"))
    )
}

write_csv(markers_top, here(resdir, "top-markers.csv"))

# end =========================================================================

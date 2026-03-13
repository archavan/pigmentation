# Arun Chavan
# Started: 2026-03-12

# setup =======================================================================
library(tidyverse)
library(here)
library(glue)
library(ggrepel)

basedir <- here("species/crob/bulk/ARCA123/KY21")
genomedir <- here("species/crob/genome/KY21")
resdir <- here(
  "species/crob/bulk/ARCA123/KY21/results/14b_functional-modules_int"
)
fs::dir_create(resdir)

source(here(basedir, "code/utilities_ggplot.R"))

###############################################################################
# data ========================================================================
###############################################################################

tpm <- read_csv(here(basedir, "data/03_compiled/tpm_no-testis.csv"))
uss <- read_csv(here(basedir, "data/03_compiled/uss_intestine.csv"))
samples <- read_csv(here(basedir, "data/00_metadata/samples.csv")) |>
  filter(tissue_short != "S")

tpm <- tpm |>
  select(-starts_with("stm")) |>
  rename_with(
    function(x) samples$tissue_name[match(x, samples$sample_name)],
    .cols = -gene_id
  )

###############################################################################
# genomic data ================================================================
###############################################################################

genes <- read_tsv(here(genomedir, "results/genes/gene-display-names.tsv"))
fnmod <- read_tsv(here(genomedir, "results/genes/functional-modules_gut.tsv"))

fn <- fnmod |>
  select(module.1, module.2, module.3, module.4, ortholog) |>
  filter(!is.na(ortholog)) |>
  mutate(
    mod = pmap_chr(
      list(module.1, module.2, module.3, module.4),
      ~ c(...) |> na.omit() |> paste0(collapse = " - ")
    )
  ) |>
  summarize(gene_ids = list(ortholog) |> unique(), .by = mod)

###############################################################################
# prepare for plotting ========================================================
###############################################################################

tpm_long <- tpm |>
  pivot_longer(cols = -gene_id, names_to = "tissue_name", values_to = "tpm") |>
  left_join(samples, by = "tissue_name")

uss_long <- uss |>
  pivot_longer(cols = -gene_id, names_to = "tissue_name", values_to = "uss") |>
  left_join(samples, by = "tissue_name")

display_names <- deframe(select(genes, gene_id, display_name))

###############################################################################
# faceted gene plots: tpm =====================================================
###############################################################################

fs::dir_create(here(resdir, "faceted-genes_tpm"))

order_genes_by_expression <- function(.genes) {
  tpm_long |>
    filter(gene_id %in% .genes) |>
    group_by(gene_id) |>
    summarize(tpm = mean(tpm)) |>
    ungroup() |>
    arrange(-tpm) |>
    pull(gene_id)
}

extract_facet_layout <- function(p) {
  q <- ggplot_build(p)
  n_facet_rows <- max(summarise_layout(q)$row)
  n_facet_cols <- max(summarise_layout(q)$col)
  list(nrows = n_facet_rows, ncols = n_facet_cols)
}

generate_module_display_name <- function(x) {
  str_replace_all(x, " - ", " → \n  ")
}

plot_genes_faceted_tpm <- function(.genes, .module = NULL, .save = FALSE) {
  .title <- generate_module_display_name(.module)
  gene_order <- order_genes_by_expression(.genes)
  pdat <- tpm_long |>
    filter(gene_id %in% .genes) |>
    mutate(gene_id = fct(gene_id, gene_order))

  if (nrow(pdat) == 0) {
    return(NULL)
  }

  pp <- ggplot(pdat, aes(tissue_short)) +
    geom_jitter(aes(y = tpm), width = 0.1, color = "grey") +
    stat_summary(aes(y = tpm, group = 1), fun = mean, geom = "line") +
    stat_summary(aes(y = tpm), fun = mean, geom = "point") +
    facet_wrap(
      vars(gene_id),
      scales = "free_y",
      labeller = labeller(gene_id = str_trunc(display_names, 70))
    ) +
    labs(title = .title, x = "segment") +
    theme_bw(base_size = 8) +
    theme(
      strip.text = element_text(size = 5, hjust = 0),
      plot.title = element_text(size = 7, face = 1),
      panel.grid.minor = element_blank()
    )

  if (.save) {
    nfacets <- extract_facet_layout(pp)
    ggsave(
      filename = here(
        resdir,
        "faceted-genes_tpm",
        paste0(.module, ".pdf")
      ),
      plot = pp,
      width = nfacets$ncols * 2.2,
      height = nfacets$nrows * 1.6 + 0.6,
    )
  } else {
    return(pp)
  }
}

walk2(
  fn$gene_ids,
  fn$mod,
  ~ plot_genes_faceted_tpm(.x, .y, .save = TRUE)
)

###############################################################################
# faceted gene plots: uss =====================================================
###############################################################################

fs::dir_create(here(resdir, "faceted-genes_uss"))

plot_genes_faceted_uss <- function(.genes, .module = NULL, .save = FALSE) {
  .title <- generate_module_display_name(.module)
  gene_order <- order_genes_by_expression(.genes)
  pdat <- uss_long |>
    filter(gene_id %in% .genes) |>
    mutate(gene_id = fct(gene_id, gene_order))

  if (nrow(pdat) == 0) {
    return(NULL)
  }

  pp <- ggplot(pdat, aes(tissue_short)) +
    geom_jitter(aes(y = uss), width = 0.1, color = "grey") +
    stat_summary(aes(y = uss, group = 1), fun = mean, geom = "line") +
    stat_summary(aes(y = uss), fun = mean, geom = "point") +
    facet_wrap(
      vars(gene_id),
      scales = "free_y",
      labeller = labeller(gene_id = str_trunc(display_names, 70))
    ) +
    labs(title = .title, x = "segment") +
    theme_bw(base_size = 8) +
    theme(
      strip.text = element_text(size = 5, hjust = 0),
      plot.title = element_text(size = 7, face = 1),
      panel.grid.minor = element_blank()
    )

  if (.save) {
    nfacets <- extract_facet_layout(pp)
    ggsave(
      filename = here(
        resdir,
        "faceted-genes_uss",
        paste0(.module, ".pdf")
      ),
      plot = pp,
      width = nfacets$ncols * 2.2,
      height = nfacets$nrows * 1.6 + 0.6,
    )
  } else {
    return(pp)
  }
}

safely_plot_uss <- possibly(
  function(.genes, .module) {
    plot_genes_faceted_uss(.genes, .module, .save = TRUE)
  },
  otherwise = NULL
)

walk2(fn$gene_ids, fn$mod, safely_plot_uss)

###############################################################################
# module heatmaps =============================================================
###############################################################################
fs::dir_create(here(resdir, "module-heatmaps_uss"))

find_expressed_genes <- function(.genes) {
  tpm_long |>
    filter(gene_id %in% .genes) |>
    group_by(gene_id) |>
    summarize(max_tpm = max(tpm)) |>
    ungroup() |>
    filter(max_tpm > 2) |>
    pull(gene_id)
}

order_genes_by_expression_pattern <- function(.pdat) {
  .pdat |>
    group_by(gene_id) |>
    slice_max(uss, n = 1, with_ties = FALSE) |>
    ungroup() |>
    arrange(tissue_short) |>
    group_by(tissue_short) |>
    arrange(desc(uss), .by_group = TRUE) |>
    ungroup() |>
    pull(gene_id) |>
    as.character() |>
    rev()
}

plot_module_heatmap_uss <- function(.genes, .module = NULL, .save = FALSE) {
  .title <- generate_module_display_name(.module)
  pdat <- uss_long |>
    filter(gene_id %in% .genes)
  pdat <- pdat |>
    mutate(gene_id = fct(gene_id, order_genes_by_expression_pattern(pdat)))

  pp <- ggplot(pdat, aes(x = tissue_short, y = gene_id, fill = uss)) +
    geom_tile(color = "white") +
    scale_y_discrete(
      labels = function(x) str_trunc(display_names[x], 90, "right"),
      position = "right"
    ) +
    scale_x_discrete(name = "segment") +
    scale_fill_gradient(low = "white", high = "blue") +
    coord_fixed() +
    labs(title = .title) +
    theme_bw(base_size = 8) +
    theme(
      legend.position = "left",
      axis.title.y = element_blank(),
      axis.text.y = element_text(hjust = 0),
      plot.title = element_text(size = 8, face = 1),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank()
    )

  if (.save) {
    ngenes <- length(unique(pdat$gene_id))
    ggsave(
      filename = here(
        resdir,
        "module-heatmaps_uss",
        paste0(.module, ".pdf")
      ),
      plot = pp,
      width = 6,
      height = (ngenes * 0.15) + 1
    )
  } else {
    return(pp)
  }
}

walk2(
  fn$gene_ids,
  fn$mod,
  ~ plot_module_heatmap_uss(.x, .y, .save = TRUE)
)

###############################################################################
# plot aggregated heatmap =====================================================
###############################################################################

calculate_aggregated_uss <- function(.genes) {
  uss_long |>
    filter(gene_id %in% .genes) |>
    group_by(tissue_short) |>
    summarize(median_uss = median(uss)) |>
    ungroup()
}

fn <- fn |> mutate(agg_uss = map(gene_ids, calculate_aggregated_uss))

quartz(
  type = "pdf",
  file = here(resdir, "module-heatmap_aggregated-uss_z-scaled.pdf"),
  width = 6.75,
  height = 7.5
)
print(
  fn |>
    unnest(agg_uss) |>
    select(mod, tissue_short, median_uss) |>
    pivot_wider(
      id_cols = "mod",
      names_from = "tissue_short",
      values_from = "median_uss"
    ) |>
    column_to_rownames("mod") |>
    pheatmap::pheatmap(cluster_cols = FALSE, cluster_rows = TRUE, scale = "row")
)
dev.off()

quartz(
  type = "pdf",
  file = here(resdir, "module-heatmap_aggregated-uss_unscaled.pdf"),
  width = 6.75,
  height = 7.5
)
print(
  fn |>
    unnest(agg_uss) |>
    select(mod, tissue_short, median_uss) |>
    pivot_wider(
      id_cols = "mod",
      names_from = "tissue_short",
      values_from = "median_uss"
    ) |>
    column_to_rownames("mod") |>
    pheatmap::pheatmap(
      cluster_cols = FALSE,
      cluster_rows = TRUE,
      scale = "none"
    )
)
dev.off()

# end =========================================================================

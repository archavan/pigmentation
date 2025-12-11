# Arun Chavan
# Started: 2025-12-01

# setup =======================================================================
library(tidyverse)
library(here)
library(glue)
library(ggrepel)

basedir <- here("species/crob/bulk/ARCA123/KY21")
genomedir <- here("species/crob/genome/KY21")
resdir <- here("species/crob/bulk/ARCA123/KY21/results/14a_functional-modules")
fs::dir_create(resdir)

source(here(basedir, "code/utilities_ggplot.R"))

###############################################################################
# data ========================================================================
###############################################################################

tpm <- read_csv(here(basedir, "data/03_compiled/tpm_no-testis.csv"))
samples <- read_csv(here(basedir, "data/00_metadata/samples.csv"))

###############################################################################
# genomic data ================================================================
###############################################################################

genes <- read_tsv(here(genomedir, "results/genes/gene-display-names.tsv"))
fnmod <- read_tsv(here(genomedir, "results/genes/functional-modules_gut.tsv"))

fnmod$module <- NA_character_
for (i in seq_len(nrow(fnmod))) {
  fnmod$module[i] <- fnmod[i, 1:4] |>
    unlist() |>
    na.omit() |>
    paste0(collapse = " - ")
}

fnmod <- relocate(fnmod, module, .after = module.4)

fnmod_nest <- fnmod |>
  filter(!is.na(ortholog)) |>
  select(module, gene_mm, ortholog, display_name) |>
  distinct() |>
  nest_by(module)

###############################################################################
# prepare for plotting ========================================================
###############################################################################

samples <- samples |>
  mutate(
    tissue_short = case_when(
      tissue == "stm.1.s" ~ "S",
      tissue == "int.1.w" ~ "W1",
      tissue == "int.2.y" ~ "Y1",
      tissue == "int.3.w" ~ "W2",
      tissue == "int.4.y" ~ "Y2"
    )
  ) |>
  mutate(tissue_short = fct(tissue_short, c("S", "W1", "Y1", "W2", "Y2")))

tpm_long <- tpm |>
  pivot_longer(cols = -gene_id, names_to = "sample_name", values_to = "tpm") |>
  left_join(samples, by = "sample_name")

display_names <- deframe(select(genes, gene_id, display_name))

###############################################################################
# faceted gene plots ==========================================================
###############################################################################

fs::dir_create(here(resdir, "faceted-genes"))

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

plot_genes_faceted <- function(.genes, .module = NULL, .save = FALSE) {
  .title <- generate_module_display_name(.module)
  gene_order <- order_genes_by_expression(.genes)
  pdat <- tpm_long |>
    filter(gene_id %in% .genes) |>
    mutate(gene_id = fct(gene_id, gene_order))

  pp <- ggplot(pdat, aes(tissue_short)) +
    geom_jitter(aes(y = tpm), width = 0.1, color = "grey") +
    stat_summary(aes(y = tpm, group = 1), fun = mean, geom = "line") +
    stat_summary(aes(y = tpm), fun = mean, geom = "point") +
    facet_wrap(
      vars(gene_id),
      scales = "free_y",
      labeller = labeller(gene_id = display_names)
    ) +
    labs(title = .title, x = "segment") +
    theme_bw(base_size = 8) +
    theme(
      strip.text = element_text(size = 5, hjust = 0),
      plot.title = element_text(size = 8, face = 1),
      panel.grid.minor = element_blank()
    )

  if (.save) {
    nfacets <- extract_facet_layout(pp)
    ggsave(
      filename = here(
        resdir,
        "faceted-genes",
        paste0(.module, ".png")
      ),
      plot = pp,
      width = nfacets$ncols * 2.2,
      height = nfacets$nrows * 1.6 + 0.6,
      dpi = 300
    )
  } else {
    return(pp)
  }
}

walk2(
  fnmod_nest$data,
  fnmod_nest$module,
  ~ plot_genes_faceted(.x$ortholog, .y, .save = TRUE)
)

###############################################################################
# module heatmaps =============================================================
###############################################################################
fs::dir_create(here(resdir, "module-heatmaps"))

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
    slice_max(ztpm, n = 1, with_ties = FALSE) |>
    ungroup() |>
    arrange(tissue_short) |>
    group_by(tissue_short) |>
    arrange(desc(ztpm), .by_group = TRUE) |>
    ungroup() |>
    pull(gene_id) |>
    as.character() |>
    rev()
}

plot_module_heatmap <- function(.genes, .module = NULL, .save = FALSE) {
  .title <- generate_module_display_name(.module)
  pdat <- tpm_long |>
    filter(gene_id %in% find_expressed_genes(.genes)) |>
    group_by(gene_id) |>
    mutate(ztpm = scale(tpm)[, 1]) |>
    ungroup()
  pdat <- pdat |>
    mutate(gene_id = fct(gene_id, order_genes_by_expression_pattern(pdat)))

  pp <- ggplot(pdat, aes(x = tissue_short, y = gene_id, fill = ztpm)) +
    geom_tile(color = "black") +
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
        "module-heatmaps",
        paste0(.module, ".png")
      ),
      plot = pp,
      width = 6,
      height = (ngenes * 0.15) + 1,
      dpi = 300
    )
  } else {
    return(pp)
  }
}

walk2(
  fnmod_nest$data,
  fnmod_nest$module,
  ~ plot_module_heatmap(.x$ortholog, .y, .save = TRUE)
)

###############################################################################
# plot aggregated (median z tpm for all genes in module) ======================
###############################################################################

calculate_aggregated_ztpm <- function(.genes) {
  tpm_long |>
    filter(gene_id %in% find_expressed_genes(.genes)) |>
    group_by(gene_id) |>
    mutate(ztpm = scale(tpm)[, 1]) |>
    ungroup() |>
    group_by(tissue_short) |>
    summarize(ztpm_agg = mean(ztpm)) |>
    ungroup()
}

fnmod_nest$agg_z <- vector("list", nrow(fnmod_nest))

for (i in seq_len(nrow(fnmod_nest))) {
  fnmod_nest$agg_z[[i]] <- fnmod_nest$data[[i]] |>
    filter(!is.na(ortholog)) |>
    pull(ortholog) |>
    unique() |>
    calculate_aggregated_ztpm()
}

fnmod_nest |>
  select(module, agg_z) |>
  unnest(col = agg_z) |>
  ggplot(aes(tissue_short, module, fill = ztpm_agg)) +
  geom_tile() +
  scale_y_discrete(position = "right") +
  scale_x_discrete(name = "segment") +
  # scale_fill_gradient(low = "white", high = "blue") +
  scale_fill_viridis_c(direction = -1, name = "mean ztpm") +
  coord_fixed() +
  theme(
    legend.position = "left",
    axis.title.y = element_blank(),
    axis.text.y = element_text(hjust = 0),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank()
  )

ggsave(
  here(resdir, "all-aggregated.png"),
  width = 8,
  height = 6,
  dpi = 300
)

# end =========================================================================

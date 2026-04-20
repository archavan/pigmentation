# Arun Chavan
# Started: 2026-04-16

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)
library(glue)
library(Seurat)
library(ggtext)
library(patchwork)

basedir <- here("figures")
genomedir <- here("species/mmus/genome")
dir_mm <- here("species/mmus/zwick")

resdir <- here(basedir, "results/03a_mmus_fn-mods")
fs::dir_create(resdir)

source(here(basedir, "code/utils_colors.R"), local = TRUE)
source(here(basedir, "code/utils_themes.R"), local = TRUE)
source(here(basedir, "code/utils_seurat-plots.R"), local = TRUE)

###############################################################################
# data ========================================================================
###############################################################################

## single cell data ===========================================================
ec <- read_rds(here(dir_mm, "data/enterocytes.rds"))
ec$domain <- factor(ec$domain, LETTERS[1:5])
ec <- SetIdent(ec, value = "domain")

## genomic data ===============================================================
fn4 <- read_tsv(here(genomedir, "gene-modules/functional-modules_gut.tsv"))

## crob module enrichment =====================================================
# We want to only plot the modules that are enriched in crob pigmented and
# unpigmented segments, and specifically just the ones that highlight the
# division of labor between uptake and metabolism.
crob_enrichment <- read_tsv(here(
  basedir,
  "results/01a_crob_fn-mods/module-expression-pattern.tsv"
))

keep_mods <- c(
  "BB hydrolysis: peptide",
  "Uptake: AA",
  "Uptake: lipid",
  "Endocytosis",
  "Uptake: carb",
  "BB hydrolysis: carb",
  "Storage: glycogen",
  "Synthesis: AA",
  "Catabolism: AA",
  "Catabolism: FA",
  "Synthesis: cholesterol",
  "Catabolism: carb"
)

fn4_categories <- crob_enrichment |>
  filter(label %in% keep_mods) |>
  mutate(
    category = case_when(
      enriched_in == "pigmented" ~ "Uptake",
      enriched_in == "unpigmented" ~ "Metabolism",
      .default = NA_character_
    )
  ) |>
  mutate(category = fct(category, c("Uptake", "Metabolism")))

###############################################################################
# module dotplot ==============================================================
###############################################################################

## prepare data ===============================================================

#' Collapse Seurat object from gene level to geneset level by aggregating counts.
#'   This is achieved by transposing the Seurat object so the genes are "cells"
#'   so that we can AggregateExpression() across genes grouped by geneset, and
#'   then transpose back so columns are cells, rows are genesets, and counts are
#'   aggregated.
#'
#' @param .obj Seurat object
#' @param .ids dataframe with one column named .group_col. Rownames should be
#'   gene ids of .obj
#' @param .group_col name of the column in .ids to group by
collapse_genesets <- function(.obj, .ids, .group_col) {
  # Helper to protect/restore underscores from Seurat mangling
  protect <- function(x) gsub("_", "..", x)
  restore <- function(x) gsub("\\.\\.", "_", x)

  shared <- intersect(rownames(.obj), rownames(.ids))
  if (length(shared) == 0) {
    stop("No shared gene ids between .obj and .ids")
  }

  x <- subset(.obj, features = shared)
  tcts <- t(LayerData(x, "counts"))
  .ids <- .ids[colnames(tcts), , drop = FALSE]

  # Protect all names from _ -> - mangling
  rownames(tcts) <- protect(rownames(tcts))
  colnames(tcts) <- protect(colnames(tcts))
  rownames(.ids) <- protect(rownames(.ids))
  .ids[[.group_col]] <- protect(.ids[[.group_col]])

  y <- CreateSeuratObject(tcts, meta.data = .ids)

  agg <- AggregateExpression(
    y,
    assays = DefaultAssay(y),
    return.seurat = TRUE,
    group.by = .group_col
  ) |>
    LayerData("counts") |>
    t()

  # Restore barcodes only (geneset names stay protected for final CreateSeuratObject)
  colnames(agg) <- restore(colnames(agg))

  if (!setequal(colnames(agg), colnames(x))) {
    stop("Mismatched barcode ids")
  }
  agg <- agg[, colnames(x), drop = FALSE]
  stopifnot(identical(colnames(agg), rownames(x[[]])))

  z <- CreateSeuratObject(agg, meta.data = x[[]])
  return(z)
}

fn4_for_collapsing <- fn4_categories |>
  inner_join(fn4 |> select(id, genes), by = "id") |>
  select(genes, id) |>
  separate_longer_delim(genes, "; ") |>
  filter(!is.na(genes)) |>
  filter(!duplicated(genes)) |>
  column_to_rownames("genes")

ec_fn4 <- collapse_genesets(ec, fn4_for_collapsing, "id") |>
  NormalizeData() |>
  SetIdent(value = "domain")

## plot =======================================================================

p_mod_dotplot <- dp_extract_dotplot_data(ec_fn4, rownames(ec_fn4)) |>
  left_join(fn4_categories, by = c("features.plot" = "id")) |>
  mutate(label = fct(label, levels = keep_mods)) |>
  ggplot() +
  geom_point(
    aes(id, label, size = pct.exp, color = avg.exp.scaled),
    shape = 16,
    stroke = 0
  ) +
  facet_wrap(
    vars(category),
    strip.position = "top",
    scales = "free_y",
    ncol = 1,
    labeller = labeller(category = function(x) glue("**{x}**"))
  ) +
  scale_size(
    range = c(0, 4),
    name = "% Cells with Expression",
    guide = guide_legend(
      theme = theme(
        legend.key.height = unit(0.6, "lines"),
        legend.key.width = unit(0.5, "lines")
      )
    )
  ) +
  scale_color_gradient(
    low = "white",
    high = clrs_gradient2$high,
    name = "z(Average Expression)",
    guide = guide_colorbar(
      theme = theme(
        legend.key.height = unit(3, "lines"),
        legend.key.width = unit(0.5, "lines")
      )
    )
  ) +
  scale_y_discrete(labels = function(y) str_trunc(y, 50), position = "right") +
  labs(x = "Domain") +
  theme_figure() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_line(linewidth = 0.1, color = "grey"),
    axis.line = element_line(linewidth = lwd_thin, color = "black"),
    axis.title.y = element_blank(),
    legend.position = "left",
    legend.title.position = "left",
    legend.title = element_markdown(
      hjust = 0.5,
      angle = 90,
      size = ts_small
    ),
    strip.clip = "off",
    strip.text = element_markdown(
      hjust = 0,
      vjust = 0,
      margin = margin(b = 1),
    )
  )

save_quartz(
  here(resdir, "module-dotplot.pdf"),
  .width = 2.2,
  .height = 2.1,
  .plot = p_mod_dotplot
)

###############################################################################
# PCA and crytp villus markers on PCA =========================================
###############################################################################

d_pca <- FetchData(ec, c("PC_1", "PC_2", "domain"))

cv_markers <- c(Crypt = "Rpl41", Villus = "Clca4a")
d_markers <- FetchData(ec, c("PC_1", "PC_2", cv_markers))

p_pca <- d_pca |>
  ggplot(aes(PC_1, PC_2, color = domain)) +
  geom_point(size = 0.4, shape = 16, stroke = 0) +
  scale_color_discrete(
    name = "Domain",
    guide = guide_legend(
      override.aes = list(size = 2),
      theme = theme(
        legend.key.size = unit(0.1, "lines"),
        legend.text.position = "bottom"
      )
    )
  )

p_markers <- d_markers |>
  rownames_to_column("barcode") |>
  pivot_longer(
    cols = all_of(unname(cv_markers)),
    names_to = "gene",
    values_to = "expression"
  ) |>
  mutate(gene = fct(gene, cv_markers)) |>
  ggplot(aes(PC_1, PC_2, z = expression)) +
  stat_summary_hex(bins = 60) +
  facet_grid(
    cols = vars(gene),
    labeller = as_labeller(function(x) {
      glue("*{x}* ({names(cv_markers[cv_markers == x])})")
    })
  ) +
  scale_fill_viridis_c(
    direction = -1,
    name = "Norm. Expression",
    guide = guide_colorbar(
      theme = theme(
        legend.key.height = unit(0.4, "lines"),
        legend.key.width = unit(3, "lines")
      )
    )
  )

wrap_plots(
  p_pca,
  p_markers,
  axes = "collect_y",
  axis_titles = "collect",
  widths = c(1, 2)
) &
  scale_x_continuous(breaks = scales::breaks_width(20)) &
  labs(x = "PC1 (Crypt to Villus)", y = "PC2 (Posterior to Anterior)") &
  theme_figure() &
  theme(
    legend.position = "bottom",
    panel.border = element_blank(),
    plot.background = element_blank(),
    plot.margin = margin(1, 1, 1, 1),
    axis.line = element_line(linewidth = lwd_thin, color = "black"),
    axis.title = element_text(size = ts_tiny),
    axis.text = element_text(size = ts_tiny),
    strip.text = element_markdown(size = ts_tiny, margin = margin(t = 1)),
    legend.title.position = "top",
    legend.title = element_text(
      hjust = 0.5,
      margin = margin(b = 2),
      size = ts_tiny
    ),
    legend.text = element_text(margin = margin(t = 2), size = ts_tiny),
    legend.box.spacing = unit(0.1, "mm"),
    legend.box.margin = margin()
  )

save_quartz(
  here(resdir, "pca_and_cv-markers-on-pca.pdf"),
  .width = 2.26,
  .height = 1.3,
  .plot = last_plot()
)

###############################################################################
# module against PC1 ==========================================================
###############################################################################

fn4_for_pc1_corr <- fn4_categories |>
  filter(str_detect(label, "Uptake|Catabolism")) |>
  separate_wider_delim(
    label,
    ": ",
    names = c(NA, "macromolecule"),
    cols_remove = FALSE
  ) |>
  select(id, label, category, macromolecule) |>
  mutate(
    macromolecule = case_when(
      macromolecule == "FA" ~ "lipid",
      .default = macromolecule
    )
  ) |>
  mutate(macromolecule = fct(macromolecule, c("AA", "lipid", "carb"))) |>
  arrange(macromolecule, category) |>
  mutate(label = fct_inorder(label))


d_pc1_corr <- inner_join(
  FetchData(ec, "PC_1") |> rownames_to_column("barcode"),
  FetchData(ec_fn4, fn4_for_pc1_corr$id) |> rownames_to_column("barcode"),
  by = "barcode"
) |>
  pivot_longer(
    cols = starts_with("gut"),
    names_to = "id",
    values_to = "expression"
  ) |>
  left_join(fn4_for_pc1_corr, by = "id") |>
  filter(expression > 0)

p_pc1_corr <- ggplot(d_pc1_corr, aes(PC_1, expression)) +
  geom_point(shape = 16, stroke = 0, size = 0.4, color = "grey", alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE, aes(color = category)) +
  facet_wrap(vars(label), scales = "free_y", nrow = 3) +
  scale_x_continuous(breaks = scales::breaks_width(20)) +
  scale_color_discrete(name = NULL) +
  labs(x = "PC1 (Crypt to Villus)", y = "Normalized Expression") +
  theme_figure() +
  theme(
    panel.border = element_blank(),
    axis.title.x = element_text(margin = margin(t = 3)),
    axis.title.y = element_text(margin = margin(r = 3)),
    axis.line = element_line(linewidth = lwd_thin, color = "black"),
    strip.text = element_text(margin = margin(b = 2)),
    legend.position = "bottom",
    legend.box.spacing = unit(0, "lines")
  )

save_quartz(
  here(resdir, "modules-vs-pc1.pdf"),
  .width = 2.2,
  .height = 3.5,
  .plot = p_pc1_corr
)

# end #########################################################################

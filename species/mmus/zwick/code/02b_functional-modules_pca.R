# Arun Chavan
# Started: 2026-04-07

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)
library(glue)
library(Seurat)
library(ggrepel)
library(ggtext)

basedir <- here("species/mmus/zwick")
genomedir <- here("species/mmus/genome")

resdir <- here(basedir, "results/02b_functional-modules_pca")
fs::dir_create(resdir)

source(here(basedir, "code/utilities_ggplot.R"))
source(here(basedir, "code/utilities_seurat-plots.R"))

###############################################################################
# data ========================================================================
###############################################################################

## single cell data ===========================================================

ec <- read_rds(here(basedir, "data/enterocytes.rds"))
ec$domain <- factor(ec$domain, LETTERS[1:5])
ec <- SetIdent(ec, value = "domain")

## genomic data ===============================================================
fn4 <- read_tsv(here(genomedir, "gene-modules/functional-modules_gut.tsv"))

fn4_mapping <- fn4 |>
  select(id, label) |>
  distinct()

# used in crob
fn4_order_manual <- tribble(
  ~label                   , ~group                           ,
  "BB hydrolysis: peptide" , "Enriched in<br>**pigmented**"   ,
  "BB hydrolysis: carb"    , "Enriched in<br>**pigmented**"   ,
  "Uptake: AA"             , "Enriched in<br>**pigmented**"   ,
  "Uptake: carb"           , "Enriched in<br>**pigmented**"   ,
  "Uptake: lipid"          , "Enriched in<br>**pigmented**"   ,
  "Endocytosis"            , "Enriched in<br>**pigmented**"   ,
  "Storage: lipid droplet" , "Enriched in<br>**pigmented**"   ,
  "Synthesis: cholesterol" , "Enriched in<br>**unpigmented**" ,
  "Bile acid metabolism"   , "Enriched in<br>**unpigmented**" ,
  "Catabolism: carb"       , "Enriched in<br>**unpigmented**" ,
  "Storage: glycogen"      , "Enriched in<br>**unpigmented**" ,
  "Synthesis: AA"          , "Enriched in<br>**unpigmented**" ,
  "Catabolism: FA"         , "Enriched in<br>**unpigmented**" ,
  "Catabolism: AA"         , "Enriched in<br>**unpigmented**" ,
  "Endocrine"              , "Enriched in<br>**unpigmented**"
)


###############################################################################
# plotting setup ==============================================================
###############################################################################

plot_dotplot2 <- function(
  .features,
  .title,
  .reorder.genes = TRUE,
  .wrap.title = TRUE
) {
  if (.wrap.title) {
    .title <- str_replace_all(.title, " - ", " → \n  ")
  }

  dd <- dp_extract_dotplot_data(ec, .features)
  if (.reorder.genes) {
    dd <- dp_define_feature_order(dd)
  }
  dp_split_dotplot_data(dd) |>
    purrr::map(dp_plot) |>
    purrr::map(dp_add_title, .title) |>
    purrr::map(function(p) {
      p +
        theme(plot.title = element_text(size = 8, face = 1))
    })
}

###############################################################################
# PCA =========================================================================
###############################################################################

ecpc <- FetchData(ec, c("domain", "PC_1", "PC_2"))

ggplot(ecpc, aes(domain, PC_1)) +
  geom_violin(fill = "grey") +
  theme_bw_2()

ggsave(here(resdir, "domain-vs-pc1.pdf"), width = 3, height = 2.5)

ggplot(ecpc, aes(domain, PC_2)) +
  geom_violin(fill = "grey") +
  theme_bw_2()

ggsave(here(resdir, "domain-vs-pc2.pdf"), width = 3, height = 2.5)

DimPlot(ec, reduction = "pca", group.by = "domain") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

ggsave(here(resdir, "PCA_domain.pdf"), width = 5, height = 3.5)

FeaturePlot(ec, "Apoa4", reduction = "pca")
ggsave(here(resdir, "PCA_Apoa4.pdf"), width = 5, height = 3.5)

###############################################################################
# plot aggregated module score on PCA =========================================
###############################################################################

## prepare ====================================================================

#' Collapse Seurat object from gene level to geneset level by aggregating counts
#'
#' @param .obj Seurat object
#' @param .ids dataframe with one column named .group_col. Rownames should be gene ids of .obj
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


fn4_for_collapsing <- fn4 |>
  select(genes, id, label) |>
  separate_longer_delim(genes, "; ") |>
  filter(!is.na(genes)) |>
  filter(!duplicated(genes)) |>
  select(genes, id) |>
  column_to_rownames("genes")

ec_fn4 <- collapse_genesets(ec, fn4_for_collapsing, "id") |>
  NormalizeData() |>
  SetIdent(value = "domain")

## plot on pca ================================================================
pc_pdat <- inner_join(
  ecpc |> rownames_to_column("barcode"),
  FetchData(ec_fn4, rownames(ec_fn4)) |> rownames_to_column("barcode"),
  by = "barcode"
) |>
  pivot_longer(
    cols = starts_with("gut"),
    names_to = "fn4_id",
    values_to = "expression"
  ) |>
  left_join(fn4_mapping, by = c("fn4_id" = "id"))

plot_module_on_pc <- function(.id) {
  pc_pdat |>
    filter(fn4_id == .id) |>
    arrange(expression) |>
    mutate(barcode = fct_inorder(barcode)) |>
    ggplot(aes(PC_1, PC_2, z = expression)) +
    stat_summary_hex(fun = median, bins = 70) +
    scale_fill_viridis_c(direction = -1, name = "Aggregated\nExpression") +
    labs(title = fn4_mapping$label[fn4_mapping$id == .id]) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 8, face = 1)
    )
}

quartz(
  type = "pdf",
  file = here(resdir, "aggregated-fnmod-on-pca.pdf"),
  width = 4.75,
  height = 3
)
for (i in rownames(ec_fn4)) {
  print(plot_module_on_pc(i))
}
dev.off()

## plot against PC1 and PC2 ===================================================
plot_module_against_pc <- function(.id, .pc) {
  pc_pdat |>
    filter(fn4_id == .id) |>
    filter(expression > 0) |>
    ggplot(aes(.data[[paste0("PC_", .pc)]], expression)) +
    geom_point(alpha = 0.5, shape = 1) +
    geom_smooth(method = "loess") +
    labs(title = fn4_mapping$label[fn4_mapping$id == .id]) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 8, face = 1)
    )
}

quartz(
  type = "pdf",
  file = here(resdir, "aggregated-fnmod-against-PC1.pdf"),
  width = 4.5,
  height = 3.5
)
for (i in rownames(ec_fn4)) {
  print(plot_module_against_pc(i, 1))
}
dev.off()

quartz(
  type = "pdf",
  file = here(resdir, "aggregated-fnmod-against-PC2.pdf"),
  width = 4.5,
  height = 3.5
)
for (i in rownames(ec_fn4)) {
  print(plot_module_against_pc(i, 2))
}
dev.off()

###############################################################################
# Aggregated module as dotplot ================================================
###############################################################################

dp_extract_dotplot_data(ec_fn4, rownames(ec_fn4)) |>
  left_join(fn4_mapping, by = c("features.plot" = "id")) |>
  filter(label %in% fn4_order_manual$label) |>
  left_join(fn4_order_manual) |>
  mutate(label = fct(label, levels = fn4_order_manual$label)) |>
  ggplot() +
  geom_point(
    aes(id, label, size = pct.exp, color = avg.exp.scaled),
    shape = 16,
    stroke = 0
  ) +
  facet_wrap(vars(group), strip.position = "top", scales = "free", ncol = 1) +
  scale_size(range = c(0, 6)) +
  scale_color_gradient(low = "white", high = pal_okabe[[5]]) +
  scale_y_discrete(labels = function(y) str_trunc(y, 50), position = "right") +
  theme_half_open_grid() +
  theme(
    strip.text.x.top = element_markdown(),
    axis.title = element_blank(),
    legend.position = "left"
  )

quartz(
  type = "pdf",
  file = here(resdir, "dotplot_aggregated_fn4-mods-crob-enriched.pdf"),
  width = 4,
  height = 5
)
print(last_plot())
dev.off()

###############################################################################
# plot crob enriched modules: individual genes ================================
###############################################################################

cr_enriched_fn4 <- fn4 |>
  inner_join(fn4_order_manual, by = "label") |>
  select(id, label, group, genes) |>
  separate_longer_delim(genes, "; ") |>
  mutate(label = fct(label, levels(fn4_order_manual$label)))

dp_extract_dotplot_data(.obj = ec, .features = unique(cr_enriched_fn4$genes)) |>
  rename(domain = id) |>
  left_join(cr_enriched_fn4, by = c("features.plot" = "genes")) |>
  mutate(max_pct = max(pct.exp), .by = features.plot) |>
  filter(max_pct > 20) |>
  ggplot() +
  geom_point(
    aes(domain, features.plot, size = pct.exp, color = avg.exp.scaled),
    shape = 16,
    stroke = 0
  ) +
  facet_wrap(
    vars(label),
    scales = "free",
    space = "free_y",
    strip.position = "top",
    ncol = 1
  ) +
  scale_size(range = c(0, 6)) +
  scale_color_gradient(low = "white", high = pal_okabe[[5]]) +
  scale_y_discrete(labels = function(y) str_trunc(y, 50), position = "right") +
  theme_half_open_grid() +
  theme(
    strip.text.x.top = element_markdown(),
    axis.title = element_blank(),
    legend.position = "left",
    strip.clip = "off"
  )

quartz(
  type = "pdf",
  file = here(resdir, "dotplot_genes_fn4-mods-crob-enriched.pdf"),
  width = 3.5,
  height = 16
)
print(last_plot())
dev.off()

# end =========================================================================

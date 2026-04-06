# Arun Chavan
# Started: 2026-04-06

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)
library(glue)
library(Seurat)
library(ggrepel)

basedir <- here("species/bflo/sc/gut/VNyyK")
genomedir <- here("species/bflo/genome/VNyyK")
scdir <- here("~/Dropbox/research/projects/endoderm/species/bflo/10x/gut/VNyyK")

resdir <- here(basedir, "results/01c_fn-mods_hgut-clusters_epith")
fs::dir_create(resdir)

source(here(basedir, "code/utilities_ggplot.R"))
source(here(basedir, "code/utilities_seurat-plots.R"))

###############################################################################
# functions ===================================================================
###############################################################################

define_cluster_order <- function(.obj) {
  hcc <- .obj %>%
    AggregateExpression(
      return.seurat = TRUE,
      features = VariableFeatures(.obj)
    ) %>%
    ScaleData() %>%
    GetAssayData(layer = "scale.data") %>%
    cor() %>%
    dist() %>%
    hclust()

  hcc$labels[hcc$order]
}

plot_dotplot2 <- function(
  .features,
  .title = NULL,
  .cluster_order = NULL,
  .reorder_genes = TRUE,
  .obj
) {
  xx <- dp_extract_dotplot_data(.obj = .obj, .features = .features)

  if (!is.null(.cluster_order)) {
    xx <- dp_define_cluster_order(xx, .cluster_order = .cluster_order)
  }

  if (.reorder_genes) {
    xx <- dp_define_feature_order(p.dat = xx)
  }

  xx |>
    dp_split_dotplot_data() |>
    purrr::map(dp_plot) |>
    purrr::map(dp_replace_gene_id_with_display_name, display_names, 80) |>
    purrr::map(dp_add_title, .title = .title)
}

#' Collapse Seurat object from gene level to HOG level by aggregating counts
#'
#' @param .obj Seurat object
#' @param .ids dataframe with one column "hog_name". Rownames should be in gene ids of .obj
collapse_hogs <- function(.obj, .ids) {
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
  .ids$hog_name <- protect(.ids$hog_name)

  y <- CreateSeuratObject(tcts, meta.data = .ids)

  agghog <- AggregateExpression(
    y,
    assays = DefaultAssay(y),
    return.seurat = TRUE,
    group.by = "hog_name"
  ) |>
    LayerData("counts") |>
    t()

  # Restore barcodes only (HOG names stay protected for final CreateSeuratObject)
  colnames(agghog) <- restore(colnames(agghog))

  if (!setequal(colnames(agghog), colnames(x))) {
    stop("Mismatched barcode ids")
  }
  agghog <- agghog[, colnames(x), drop = FALSE]
  stopifnot(identical(colnames(agghog), rownames(x[[]])))

  z <- CreateSeuratObject(agghog, meta.data = x[[]])
  return(z)
}

define_hog_order <- function(p.dat) {
  .hog_order <- p.dat %>%
    mutate(combined.exp = pct.exp * avg.exp.scaled) %>%
    group_by(hog_label) %>%
    slice_max(combined.exp, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(id) %>%
    group_by(id) %>%
    arrange(desc(combined.exp), .by_group = TRUE) %>%
    ungroup() %>%
    pull(hog_label) %>%
    as.character() %>%
    rev()
  p.dat %>%
    mutate(hog_label = fct(x = as.character(hog_label), levels = .hog_order))
}


###############################################################################
# data ========================================================================
###############################################################################

## genomic data ===============================================================
genes <- read_tsv(here(genomedir, "results/genes/gene-display-names.tsv"))
display_names <- genes |>
  select(gene_id, display_name) |>
  deframe()

fn4 <- read_tsv(here(genomedir, "results/genes/functional-modules_gut.tsv"))

fn4_hogs <- fn4 |>
  filter(!is.na(hog_label)) |>
  select(ortholog, hog_name = hog_label) |>
  distinct() |>
  column_to_rownames("ortholog")

fn4_order_manual <- tribble(
  ~label                   , ~group                           ,
  "BB hydrolysis: peptide" , "Enriched in<br>**pigmented**"   ,
  "BB hydrolysis: carb"    , "Enriched in<br>**pigmented**"   ,
  "Uptake: AA"             , "Enriched in<br>**pigmented**"   ,
  "Uptake: carb"           , "Enriched in<br>**pigmented**"   ,
  "Uptake: lipid"          , "Enriched in<br>**pigmented**"   ,
  "Endocytosis"            , "Enriched in<br>**pigmented**"   ,
  "Storage: lipid droplet" , "Enriched in<br>**pigmented**"   ,
  "Storage: glycogen"      , "Enriched in<br>**unpigmented**" ,
  "Bile acid metabolism"   , "Enriched in<br>**unpigmented**" ,
  "Catabolism: AA"         , "Enriched in<br>**unpigmented**" ,
  "Catabolism: FA"         , "Enriched in<br>**unpigmented**" ,
  "Catabolism: carb"       , "Enriched in<br>**unpigmented**" ,
  "Synthesis: AA"          , "Enriched in<br>**unpigmented**" ,
  "Synthesis: cholesterol" , "Enriched in<br>**unpigmented**" ,
  "Endocrine"              , "Enriched in<br>**unpigmented**"
)

## single-cell data ===========================================================
bf <- read_rds(here(scdir, "data/02_processed/03_clustered/clustered_v1.rds"))

tissue_of_origin <- read_csv(here(
  scdir,
  "results/03_clustering/r1/tissue-of-origin.csv"
))

bf_fn4 <- subset(bf, features = fn4$ortholog)

## subset to hgut clusters ====================================================
hgut_clusters <- tissue_of_origin |>
  filter(tissue == "gut") |>
  pull(cluster_v1) |>
  as.character()
hgut_epith <- hgut_clusters[
  hgut_clusters %in%
    c("2", "4", "7", "9", "11", "14", "16", "23", "28", "31", "32")
]

bf_fn4 <- subset(bf_fn4, idents = hgut_epith)

bf_fn4_order <- define_cluster_order(bf_fn4) |> str_replace("g", "")

###############################################################################
# module genes ================================================================
###############################################################################

fs::dir_create(here(resdir, "functional-modules/dotplots"))

fn4_nest <- fn4 |>
  group_by(id, label) |>
  nest() |>
  ungroup()

for (i in seq_len(nrow(fn4_nest))) {
  orthologs <- fn4_nest$data[[i]]$ortholog |>
    unique() |>
    na.omit() |>
    as.character()

  if (length(orthologs) == 0) {
    next
  }

  make_safe_filename <- function(id, label) {
    glue("{id}_{label}.pdf") |>
      str_replace_all(": ", "_") |>
      str_replace_all(" ", ".")
  }

  plot_dotplot2(
    .features = orthologs,
    .title = paste0(fn4_nest$id[[i]], " ", fn4_nest$label[[i]]),
    .cluster_order = bf_fn4_order,
    .reorder_genes = FALSE,
    .obj = bf_fn4
  ) |>
    save_dotplot(
      here(
        resdir,
        "functional-modules/dotplots",
        make_safe_filename(fn4_nest$id[[i]], fn4_nest$label[[i]])
      ),
      height = 10,
      width = 8.75
    )
  rm(orthologs)
}

###############################################################################
# plot individual HOGs ========================================================
###############################################################################

bf_fn4_hog <- collapse_hogs(bf_fn4, fn4_hogs) |>
  SetIdent(value = "cluster_v1") |>
  NormalizeData()

bf_fn4_hog_order <- define_cluster_order(bf_fn4_hog) |>
  str_replace("g", "")

bf_fn4_hog_dp_dat <- dp_extract_dotplot_data(
  bf_fn4_hog,
  Features(bf_fn4_hog)
) |>
  dp_define_cluster_order(bf_fn4_hog_order) |>
  dp_define_feature_order() |>
  as_tibble() |>
  rename(hog_label = features.plot, cluster_int = id) |>
  mutate(hog_label = gsub("\\.\\.", "_", hog_label))

fn4_hog_dp_data <- fn4 |>
  filter(!is.na(hog_label)) |>
  select(id, label, hog_label) |>
  distinct() |>
  left_join(bf_fn4_hog_dp_dat, by = 'hog_label') |>
  filter(!is.na(avg.exp))

fn4_hog_dp_data <- define_hog_order(fn4_hog_dp_data)

p_fn4_hogs <- ggplot(fn4_hog_dp_data) +
  geom_point(
    aes(cluster_int, hog_label, size = pct.exp, color = avg.exp.scaled),
    shape = 16,
    stroke = 0
  ) +
  facet_wrap(vars(label), space = "free_y", scale = "free") +
  scale_size(range = c(0, 6)) +
  scale_color_gradient(low = "white", high = pal_okabe[[5]]) +
  scale_y_discrete(labels = function(y) str_trunc(y, 50), position = "right") +
  theme_half_open_grid() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_blank(),
    panel.grid.major = element_line(linewidth = 0.1),
    legend.position = "left"
  )

quartz(
  type = "pdf",
  file = here(resdir, "functional-modules/fn4-hogs-by-cluster.pdf"),
  width = 7.75,
  height = 50
)
print(p_fn4_hogs)
dev.off()

###############################################################################
# plot individual genes =======================================================
###############################################################################

fn4_genes_with_hog <- fn4 |>
  select(id, label, gene_id = ortholog, hog_label) |>
  filter(!is.na(hog_label)) |>
  distinct()

bf_fn4 <- SetIdent(bf_fn4, value = "cluster_v1")
fn4_genes_with_hog_pdat <- dp_extract_dotplot_data(
  bf_fn4,
  unique(fn4_genes_with_hog$gene_id)
) |>
  dp_define_feature_order() |>
  as_tibble() |>
  rename(gene_id = features.plot, cluster_int = id)

p_fn4_genes <- fn4_genes_with_hog |>
  left_join(fn4_genes_with_hog_pdat, by = "gene_id") |>
  filter(!is.na(avg.exp)) |>
  filter(
    gene_id %in%
      (fn4_genes_with_hog_pdat |>
        summarize(max_pct_exp = max(pct.exp), .by = gene_id) |>
        filter(max_pct_exp > 20) |>
        pull(gene_id))
  ) |>
  filter(label %in% fn4_order_manual$label) |>
  mutate(label = fct(label, fn4_order_manual$label)) |>
  mutate(gene_id = fct(gene_id, levels(fn4_genes_with_hog_pdat$gene_id))) |>
  arrange(gene_id) |>
  mutate(display_name = paste0(gene_id, "_", hog_label)) |>
  mutate(display_name = fct_inorder(display_name)) |>
  mutate(cluster_int = fct(as.character(cluster_int), bf_fn4_order)) |>
  ggplot() +
  geom_point(
    aes(cluster_int, display_name, size = pct.exp, color = avg.exp.scaled),
    shape = 16,
    stroke = 0
  ) +
  facet_wrap(vars(label), space = "free_y", scale = "free") +
  scale_size(range = c(0, 6)) +
  scale_color_gradient(low = "white", high = pal_okabe[[5]]) +
  scale_y_discrete(labels = function(y) str_trunc(y, 50), position = "right") +
  theme_half_open_grid() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_blank(),
    panel.grid.major = element_line(linewidth = 0.1),
    legend.position = "left",
    strip.clip = "off"
  )

quartz(
  type = "pdf",
  file = here(resdir, "functional-modules/fn4-genes-by-cluster.pdf"),
  width = 7.75,
  height = 15
)
print(p_fn4_genes)
dev.off()

# end #########################################################################

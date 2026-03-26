# Arun Chavan
# Started: 2026-03-23

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)
library(glue)
library(Seurat)
library(ggrepel)
library(harmony)
library(ggsankey)
library(patchwork)
library(UCell)

basedir <- here("species/crob/sc/run-02/KY21")
genomedir <- here("species/crob/genome/KY21")
datadir <- here(
  "~/Dropbox/research/projects/endoderm/species/crob/10x/run-02/KY21"
)
bulkdir <- here("species/crob/bulk/ARCA123/KY21")

resdir <- here(basedir, "results/01a_intestine-subclustering")
fs::dir_create(resdir)

source(here(basedir, "code/utilities_ggplot.R"))
source(here(basedir, "code/utilities_clustering.R"))
source(here(basedir, "code/utilities_seurat-plots.R"))

barebones_seurat <- function(.obj) {
  cts <- GetAssayData(.obj, assay = "RNA", layer = "counts")
  md <- .obj[[]]
  new_obj <- CreateSeuratObject(counts = cts, meta.data = md)
  return(new_obj)
}

cpt_preprocess <- function(obj_merged) {
  obj_merged <- obj_merged |>
    NormalizeData() |>
    FindVariableFeatures(selection.method = "vst", nfeatures = 3000) |>
    ScaleData() |>
    RunPCA(npcs = 30)
  return(obj_merged)
}

add_prefix_to_idents <- function(.obj, .cluster.name, .prefix) {
  old_levels <- levels(.obj[[.cluster.name, drop = TRUE]])
  .obj[[.cluster.name]] <- paste0(.prefix, .obj[[.cluster.name, drop = TRUE]])
  .obj[[.cluster.name]] <- factor(
    .obj[[.cluster.name, drop = TRUE]],
    paste0(.prefix, old_levels)
  )
  .obj <- SetIdent(.obj, value = .cluster.name)
  return(.obj)
}

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
  .cluster_order = int_clust_order,
  .reorder_genes = TRUE,
  .obj = int
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
    purrr::map(dp_replace_gene_id_with_display_name, display_names, 90) |>
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
    LayerData("counts") %>%
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

###############################################################################
# data ========================================================================
###############################################################################

## genomic data ===============================================================
genes <- read_tsv(here(genomedir, "results/genes/gene-display-names.tsv"))
display_names <- genes |>
  select(gene_id, display_name) |>
  deframe()

testis_genes <- read_csv(here(
  bulkdir,
  "results/04b_testis-genes/testis-enriched-genes.csv"
))

fn4 <- read_tsv(here(genomedir, "results/genes/functional-modules_gut.tsv"))

fn4_hogs <- fn4 |>
  filter(!is.na(hog_label)) |>
  select(ortholog, hog_name = hog_label) |>
  distinct() |>
  column_to_rownames("ortholog")

## single cell data ===========================================================
cr <- read_rds(here(datadir, "data/processed/03_clustered/clustered_v1.rds"))

cr <- PercentageFeatureSet(
  object = cr,
  features = testis_genes$gene_id,
  col.name = "pct.testis"
)

FeaturePlot(cr, "pct.testis", order = TRUE)
VlnPlot(cr, "pct.testis")

## bulk rna-seq data ==========================================================
tpm <- read_csv(here(bulkdir, "data/03_compiled/tpm_no-testis.csv"))
uss <- read_csv(here(bulkdir, "data/03_compiled/uss_intestine.csv"))
bulk_samples <- read_csv(here(bulkdir, "data/00_metadata/samples.csv"))

tpm <- tpm |>
  rename_with(
    \(x) bulk_samples$tissue_name[match(x, bulk_samples$sample_name)],
    .cols = -gene_id
  )

segment_de <- tibble(
  segment = c("I1", "I2", "I3", "I4"),
  tissue = c("int.1.w", "int.2.y", "int.3.w", "int.4.y")
) |>
  mutate(
    de_up = map(tissue, \(t) {
      read_csv(here(bulkdir, "results/09a_segment-specific", t, "genes/up.csv"))
    })
  )

###############################################################################
# subclustering: round 1 ======================================================
###############################################################################

## cluster ====================================================================

int0 <- subset(cr, subset = tissue == "int")
int0 <- barebones_seurat(int0)
int0 <- cpt_preprocess(int0)
int0 <- RunHarmony(int0, group.by.vars = "animal_id")
int0 <- RunUMAP(int0, reduction = "harmony", dims = 1:30, n.neighbors = 20L)
int0 <- FindNeighbors(int0, reduction = "harmony", dims = 1:30, k.param = 20L)
int0 <- FindClusters(int0, algorithm = 2, cluster.name = "cluster_int0")
int0 <- add_prefix_to_idents(int0, "cluster_int0", "ii")

int0_clust_order <- define_cluster_order(int0)

## find markers ===============================================================
markers_raw0 <- FindAllMarkers(int0, only.pos = TRUE)
markers_raw0 |>
  filter(p_val_adj < 0.05, avg_log2FC > 1) |>
  group_by(cluster) |>
  slice_head(n = 20) |>
  ungroup() |>
  pull(gene) |>
  unique() |>
  plot_dotplot2(.obj = int0, .cluster_order = int0_clust_order)

markers_raw0 |>
  filter(p_val_adj < 0.05, avg_log2FC > 1) |>
  filter(cluster == "ii9") |>
  slice_head(n = 50) |>
  pull(gene) |>
  unique() |>
  plot_dotplot2(.obj = int0, .cluster_order = int0_clust_order)

VlnPlot(int0, "pct.testis")
FeaturePlot(int0, "pct.testis", order = TRUE)
DimPlot(int0, label = TRUE)

FetchData(int0, c("umap_1", "umap_2", "pct.testis", "cluster_int0")) |>
  as.data.frame() |>
  rownames_to_column("barcode") |>
  as_tibble() |>
  ggplot(aes(cluster_int0, pct.testis)) +
  geom_violin()

# clusters 1 are testis cells. Also remove cell with high fraction of
# testis genes. Removes only a small numnber of cells from other clusters.
clean_barcodes0 <- int0[[]] |>
  rownames_to_column("barcode") |>
  as_tibble() |>
  filter(
    pct.testis < 10 &
      !cluster_int0 %in% c("ii1")
  ) |>
  pull(barcode)

###############################################################################
# suclustering: round 2 =======================================================
###############################################################################

int <- subset(cr, cells = clean_barcodes0)

nonepithelial_or_stm_cells <- int[[]] |>
  filter(
    cluster_v1 %in%
      c(
        1,
        16,
        6,
        19,
        10,
        11,
        14,
        15,
        17,
        18,
        20,
        22,
        25,
        26,
        27,
        28,
        32,
        34,
        35,
        36
      )
  ) |>
  rownames()


# int_epi_cells <- int[[]] |>
#   filter(cluster_v1 %in% c(3, 9, 24)) |>
#   rownames()

# int <- subset(int, cells = int_epi_cells)
int <- subset(int, cells = setdiff(Cells(int), nonepithelial_or_stm_cells))

FeaturePlot(int, "pct.testis", order = TRUE)
VlnPlot(int, "pct.testis")

## cluster ===================================================================
int <- barebones_seurat(int)
int <- cpt_preprocess(int)
int <- RunHarmony(int, group.by.vars = "animal_id")
int <- RunUMAP(int, reduction = "harmony", dims = 1:30, n.neighbors = 20L)
int <- FindNeighbors(int, reduction = "harmony", dims = 1:30, k.param = 20L)
int <- FindClusters(int, algorithm = 2, cluster.name = "cluster_int")
int <- add_prefix_to_idents(int, "cluster_int", "i")

int_clust_order <- define_cluster_order(int)

## save dim plots =============================================================
fs::dir_create(here(resdir, "clustering"))

DimPlot(int, label = TRUE, pt.size = 0.5)
ggsave(here(resdir, "clustering/umap.pdf"), width = 6, height = 4.5)

DimPlot(int, label = TRUE, pt.size = 0.5, split.by = "animal_id") +
  theme_bw_2()
ggsave(here(resdir, "clustering/umap_split.pdf"), width = 8, height = 4)

DimPlot(
  int,
  group.by = c("cluster_v1", "cluster_int"),
  label = TRUE,
  pt.size = 0.5
)
ggsave(
  here(resdir, "clustering/umap_colored-by-full-data-clusters.png"),
  width = 12,
  height = 6
)

## markers ====================================================================
fs::dir_create(here(resdir, "markers/all"))
fs::dir_create(here(resdir, "markers/top50"))

markers_raw <- FindAllMarkers(int, only.pos = TRUE)
markers_top <- markers_raw |>
  as_tibble() |>
  filter(p_val_adj < 0.05, avg_log2FC > 1)

for (i in unique(markers_top$cluster)) {
  markers_top |>
    filter(cluster == i) |>
    slice_head(n = 50) |>
    pull(gene) |>
    plot_dotplot2(.reorder_genes = FALSE) |>
    save_dotplot(
      here(resdir, "markers/top50", paste0("cluster-", i, ".pdf"))
    )
}

for (i in unique(markers_top$cluster)) {
  markers_top |>
    filter(cluster == i) |>
    pull(gene) |>
    plot_dotplot2(.reorder_genes = FALSE) |>
    save_dotplot(
      here(resdir, "markers/all", paste0("cluster-", i, ".pdf"))
    )
}

markers_top |>
  group_by(cluster) |>
  slice_head(n = 4) |>
  ungroup() |>
  pull(gene) |>
  unique() |>
  plot_dotplot2() |>
  save_dotplot(here(resdir, "markers", "top4.pdf"), height = 9.5)

markers_top |>
  rename(gene_id = gene) |>
  left_join(genes |> select(gene_id, display_name), by = "gene_id") |>
  write_tsv(here(resdir, "markers/markers_top.tsv"))

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

  plot_dotplot2(
    .features = orthologs,
    .title = paste0(fn4_nest$id[[i]], " ", fn4_nest$label[[i]]),
    .reorder_genes = FALSE
  ) |>
    save_dotplot(
      here(
        resdir,
        "functional-modules/dotplots",
        paste0(fn4_nest$id[[i]], " ", fn4_nest$label[[i]], ".pdf")
      ),
      height = 10,
      width = 9
    )
  rm(orthologs)
}

###############################################################################
# segment specific genes from bulk RNA-seq data ===============================
###############################################################################

# manually selected from yellow vs white analysis

segment_genes <- list(
  W1 = paste0(
    "KY21.Chr",
    c(
      "3.1527",
      "3.670",
      "4.262",
      "10.314"
    )
  ),
  W2 = paste0(
    "KY21.Chr",
    c(
      "8.1322",
      "9.489",
      "9.712",
      "6.67",
      "1.1867",
      "9.128",
      "6.393",
      "1.29",
      "13.312",
      "12.802",
      "10.579",
      "14.873",
      "3.69",
      "9.927",
      "5.529",
      "9.14"
    )
  ),
  Y1 = paste0(
    "KY21.Chr",
    c(
      "11.1215",
      "10.1002",
      "11.820",
      "11.819",
      "7.1034",
      "6.543",
      "1.1648"
    )
  ),
  Y2 = c(
    paste0(
      "KY21.Chr",
      c(
        "10.784",
        "1.556",
        "1.649",
        "12.696",
        "1.156",
        "11.1051",
        "13.377",
        "1.1917",
        "11.706",
        "2.1450",
        "13.350",
        "9.461",
        "5.768"
      )
    ),
    "KY21.UAContig20.2"
  )
)

fs::dir_create(here(resdir, "segment-genes"))

iwalk(segment_genes, function(x, y) {
  plot_dotplot2(x, .title = paste0(y, " genes")) |>
    save_dotplot(
      here(resdir, "segment-genes", paste0(y, "-genes.pdf")),
      width = 9
    )
})

###############################################################################
# pseudobulk correlation with bulk RNA-seq segments ===========================
###############################################################################

fs::dir_create(here(resdir, "pseudobulk-correlation"))

# 1. Pseudobulk from scRNA-seq
pseudobulk_obj <- AggregateExpression(
  int,
  assays = "RNA",
  return.seurat = TRUE,
  verbose = FALSE
)
pseudobulk <- GetAssayData(pseudobulk_obj, layer = "data") |> as.matrix()

# 2. Bulk mean TPM per segment (log1p to match scRNA-seq log-normalization)
bulk_mean <- tpm |>
  pivot_longer(-gene_id, names_to = "tissue_name", values_to = "tpm") |>
  separate(tissue_name, into = c("tissue_short", "rep_id"), sep = "_") |>
  summarize(tpm = mean(tpm), .by = c(gene_id, tissue_short)) |>
  pivot_wider(names_from = tissue_short, values_from = tpm) |>
  column_to_rownames("gene_id") |>
  as.matrix() |>
  log1p()

# 3. Shared variable features
shared_genes <- intersect(VariableFeatures(int), rownames(bulk_mean))

# 4. Spearman correlation (cluster × segment)
cor_mat <- cor(
  pseudobulk[shared_genes, ],
  bulk_mean[shared_genes, ],
  method = "spearman"
)

# 5. Heatmap
segment_order <- c("S", "I1", "I2", "I3", "I4")

cor_long <- cor_mat |>
  as.data.frame() |>
  rownames_to_column("cluster") |>
  pivot_longer(-cluster, names_to = "segment", values_to = "correlation") |>
  mutate(
    cluster = factor(cluster, int_clust_order),
    segment = factor(segment, segment_order)
  )

p_pseudobulk_cor <- ggplot(cor_long, aes(segment, cluster, fill = correlation)) +
  geom_tile() +
  geom_text(aes(label = round(correlation, 2)), size = 3) +
  scale_fill_viridis_c(direction = -1) +
  theme_bw_2() +
  labs(
    title = "Pseudobulk\u2013bulk Spearman correlation",
    x = "Bulk segment",
    y = "scRNA-seq cluster"
  )

ggsave(
  here(resdir, "pseudobulk-correlation/heatmap.pdf"),
  p_pseudobulk_cor,
  width = 6,
  height = 5
)

write_csv(
  cor_mat |> as.data.frame() |> rownames_to_column("cluster"),
  here(resdir, "pseudobulk-correlation/correlation-matrix.csv")
)

###############################################################################
# segment scoring with UCell ==================================================
###############################################################################

fs::dir_create(here(resdir, "segment-scoring"))

# 1. Filter DE genes by USS specificity
#    For each segment, keep only upregulated genes where mean USS >= 0.5
uss_long <- uss |>
  pivot_longer(-gene_id, names_to = "tissue_name", values_to = "uss") |>
  separate(tissue_name, into = c("tissue_short", "rep_id"), sep = "_") |>
  summarize(mean_uss = mean(uss), .by = c(gene_id, tissue_short))

segment_signatures <- segment_de |>
  mutate(
    sig_genes = map2(de_up, segment, \(de, seg) {
      uss_seg <- uss_long |> filter(tissue_short == seg)
      de |>
        inner_join(uss_seg, by = "gene_id") |>
        filter(mean_uss >= 0.5) |>
        pull(gene_id)
    })
  )

sig_list <- set_names(segment_signatures$sig_genes, segment_signatures$segment)

# 2. Score each cell with UCell
int <- AddModuleScore_UCell(int, features = sig_list)

# 3. Visualize: Feature plots
ucell_features <- paste0(names(sig_list), "_UCell")

feature_plots <- FeaturePlot(
  int,
  features = ucell_features,
  order = TRUE,
  combine = FALSE
)

feature_patch <- wrap_plots(feature_plots, ncol = 2)

ggsave(
  here(resdir, "segment-scoring/ucell-feature-plots.pdf"),
  feature_patch,
  width = 10,
  height = 10
)

# 4. Visualize: Violin plots per cluster
violin_data <- FetchData(int, c(ucell_features, "cluster_int")) |>
  pivot_longer(-cluster_int, names_to = "segment", values_to = "score") |>
  mutate(
    segment = gsub("_UCell", "", segment),
    cluster_int = factor(cluster_int, int_clust_order)
  )

p_ucell_violin <- ggplot(violin_data, aes(cluster_int, score, fill = segment)) +
  geom_violin(scale = "width", linewidth = 0.3) +
  facet_wrap(vars(segment), ncol = 1) +
  theme_bw_2() +
  labs(x = "Cluster", y = "UCell score") +
  theme(legend.position = "none")

ggsave(
  here(resdir, "segment-scoring/ucell-violin-plots.pdf"),
  p_ucell_violin,
  width = 7,
  height = 8
)

# 5. Summary heatmap: mean UCell score per cluster × segment
ucell_summary <- violin_data |>
  summarize(mean_score = mean(score), .by = c(cluster_int, segment)) |>
  mutate(segment = factor(segment, c("I1", "I2", "I3", "I4")))

p_ucell_heatmap <- ggplot(ucell_summary, aes(segment, cluster_int, fill = mean_score)) +
  geom_tile() +
  geom_text(aes(label = round(mean_score, 3)), size = 3) +
  scale_fill_gradient(low = "white", high = "red") +
  theme_bw_2() +
  labs(
    title = "Mean UCell segment score per cluster",
    x = "Bulk segment",
    y = "scRNA-seq cluster"
  )

ggsave(
  here(resdir, "segment-scoring/ucell-heatmap.pdf"),
  p_ucell_heatmap,
  width = 5,
  height = 5
)

###############################################################################
# per-cell correlation with bulk segments =====================================
###############################################################################

fs::dir_create(here(resdir, "per-cell-correlation"))

# 1. Get per-cell normalized expression (log1p CPM from "data" layer)
cell_exp <- GetAssayData(int, layer = "data") |> as.matrix()

# 2. Bulk mean TPM per segment (already computed as bulk_mean above)
#    Use the same shared variable features
shared_genes_bulk <- intersect(rownames(cell_exp), rownames(bulk_mean))

# 3. Spearman correlation of each cell against each bulk segment
#    Only use variable features for speed and informativeness
shared_var <- intersect(VariableFeatures(int), shared_genes_bulk)

cell_bulk_cor <- cor(
  cell_exp[shared_var, ],
  bulk_mean[shared_var, ],
  method = "spearman"
)

# 4. Add correlations and best-matching segment to metadata
int$best_segment <- colnames(cell_bulk_cor)[apply(cell_bulk_cor, 1, which.max)]
for (seg in colnames(cell_bulk_cor)) {
  int[[paste0("cor_", seg)]] <- cell_bulk_cor[, seg]
}

# 5. UMAP colored by best segment
p_best_seg <- DimPlot(int, group.by = "best_segment", pt.size = 0.5) +
  ggtitle("Best-matching bulk segment per cell")

ggsave(
  here(resdir, "per-cell-correlation/umap-best-segment.pdf"),
  p_best_seg,
  width = 7, height = 5
)

# 6. Feature plots of per-cell correlations
cor_features <- paste0("cor_", colnames(cell_bulk_cor))

cor_feature_plots <- FeaturePlot(
  int,
  features = cor_features,
  order = TRUE,
  combine = FALSE
)

cor_patch <- wrap_plots(cor_feature_plots, ncol = 2)

ggsave(
  here(resdir, "per-cell-correlation/umap-correlations.pdf"),
  cor_patch,
  width = 10, height = 10
)

# 7. Proportion of cells per cluster assigned to each segment
prop_data <- int[[]] |>
  count(cluster_int, best_segment) |>
  group_by(cluster_int) |>
  mutate(prop = n / sum(n)) |>
  ungroup() |>
  mutate(
    cluster_int = factor(cluster_int, int_clust_order),
    best_segment = factor(best_segment, c("S", "I1", "I2", "I3", "I4"))
  )

p_seg_prop <- ggplot(prop_data, aes(cluster_int, prop, fill = best_segment)) +
  geom_col() +
  theme_bw_2() +
  labs(x = "Cluster", y = "Proportion of cells", fill = "Segment")

ggsave(
  here(resdir, "per-cell-correlation/cluster-segment-proportions.pdf"),
  p_seg_prop,
  width = 7, height = 4
)

# 8. Heatmap: mean correlation per cluster × segment
cor_summary <- int[[]] |>
  select(cluster_int, starts_with("cor_")) |>
  pivot_longer(-cluster_int, names_to = "segment", values_to = "correlation") |>
  mutate(segment = gsub("cor_", "", segment)) |>
  summarize(mean_cor = mean(correlation), .by = c(cluster_int, segment)) |>
  mutate(
    cluster_int = factor(cluster_int, int_clust_order),
    segment = factor(segment, c("S", "I1", "I2", "I3", "I4"))
  )

p_cor_heatmap <- ggplot(cor_summary, aes(segment, cluster_int, fill = mean_cor)) +
  geom_tile() +
  geom_text(aes(label = round(mean_cor, 3)), size = 3) +
  scale_fill_viridis_c(direction = -1) +
  theme_bw_2() +
  labs(
    title = "Mean per-cell Spearman correlation with bulk segments",
    x = "Bulk segment", y = "scRNA-seq cluster"
  )

ggsave(
  here(resdir, "per-cell-correlation/correlation-heatmap.pdf"),
  p_cor_heatmap,
  width = 6, height = 5
)

## biological sanity check: fn4 modules by assigned segment ===================
fs::dir_create(here(resdir, "per-cell-correlation/fn4-by-segment"))

# Set best_segment as identity for dotplot extraction
int_by_seg <- SetIdent(int, value = "best_segment")
seg_order <- c("S", "I1", "I2", "I3", "I4")

# Individual gene dotplot grouped by segment
fn4_seg_pdat <- dp_extract_dotplot_data(
  int_by_seg,
  unique(na.omit(fn4$ortholog))
) |>
  dp_define_cluster_order(seg_order) |>
  dp_define_feature_order() |>
  as_tibble() |>
  rename(gene_id = features.plot, segment = id)

fn4_seg_plot_data <- fn4 |>
  select(id, label, gene_id = ortholog, hog_label) |>
  filter(!is.na(hog_label)) |>
  distinct() |>
  left_join(fn4_seg_pdat, by = "gene_id") |>
  filter(!is.na(avg.exp)) |>
  filter(
    gene_id %in%
      (fn4_seg_pdat |>
        summarize(max_pct_exp = max(pct.exp), .by = gene_id) |>
        filter(max_pct_exp > 20) |>
        pull(gene_id))
  ) |>
  filter(label %in% fn4_order_manual$label) |>
  mutate(label = fct(label, fn4_order_manual$label)) |>
  mutate(gene_id = fct(gene_id, levels(fn4_seg_pdat$gene_id))) |>
  arrange(gene_id) |>
  mutate(display_name = paste0(gene_id, "_", hog_label)) |>
  mutate(display_name = fct_inorder(display_name)) |>
  mutate(segment = fct(as.character(segment), seg_order))

p_fn4_gene_seg <- ggplot(fn4_seg_plot_data) +
  geom_point(
    aes(segment, display_name, size = pct.exp, color = avg.exp.scaled),
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

ggsave(
  here(resdir, "per-cell-correlation/fn4-by-segment/fn4-genes-by-segment.pdf"),
  p_fn4_gene_seg,
  width = 6, height = 12
)

# HOG-level dotplot grouped by segment
int_hog_by_seg <- collapse_hogs(int_by_seg, fn4_hogs) |>
  SetIdent(value = "best_segment") |>
  NormalizeData()

int_hog_seg_pdat <- dp_extract_dotplot_data(
  int_hog_by_seg,
  Features(int_hog_by_seg)
) |>
  dp_define_cluster_order(seg_order) |>
  dp_define_feature_order() |>
  as_tibble() |>
  rename(hog_label = features.plot, segment = id) |>
  mutate(hog_label = gsub("\\.\\.", "_", hog_label))

fn4_hog_seg_data <- fn4 |>
  filter(!is.na(hog_label)) |>
  select(id, label, hog_label) |>
  distinct() |>
  left_join(int_hog_seg_pdat, by = "hog_label") |>
  filter(!is.na(avg.exp)) |>
  filter(label %in% fn4_order_manual$label) |>
  mutate(label = fct(label, fn4_order_manual$label))

fn4_hog_seg_data <- define_hog_order(fn4_hog_seg_data)

p_fn4_hog_seg <- ggplot(fn4_hog_seg_data) +
  geom_point(
    aes(segment, hog_label, size = pct.exp, color = avg.exp.scaled),
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

ggsave(
  here(resdir, "per-cell-correlation/fn4-by-segment/fn4-hogs-by-segment.pdf"),
  p_fn4_hog_seg,
  width = 6, height = 40
)

###############################################################################
# plot individual HOGs ========================================================
###############################################################################

int_hog <- collapse_hogs(int, fn4_hogs) |>
  SetIdent(value = "cluster_int") |>
  NormalizeData()

int_hog_cluster_order <- define_cluster_order(int_hog)

int_hog_dp_dat <- dp_extract_dotplot_data(int_hog, Features(int_hog)) |>
  dp_define_cluster_order(int_hog_cluster_order) |>
  dp_define_feature_order() |>
  as_tibble() |>
  rename(hog_label = features.plot, cluster_int = id) |>
  mutate(hog_label = gsub("\\.\\.", "_", hog_label))

fn4_hog_dp_data <- fn4 |>
  filter(!is.na(hog_label)) |>
  select(id, label, hog_label) |>
  distinct() |>
  left_join(int_hog_dp_dat, by = 'hog_label') |>
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

quartz(type = "pdf", file = here(resdir, "functional-modules/fn4-hogs-by-cluster.pdf"), width = 7, height = 50)
print(p_fn4_hogs)
dev.off()

###############################################################################
# plot individual genes =======================================================
###############################################################################

fn4_genes_with_hog <- fn4 |>
  select(id, label, gene_id = ortholog, hog_label) |>
  filter(!is.na(hog_label)) |>
  distinct()

int <- SetIdent(int, value = "cluster_int")
fn4_genes_with_hog_pdat <- dp_extract_dotplot_data(
  int,
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
  mutate(cluster_int = fct(as.character(cluster_int), int_clust_order)) |>
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

quartz(type = "pdf", file = here(resdir, "functional-modules/fn4-genes-by-cluster.pdf"), width = 6, height = 15)
print(p_fn4_genes)
dev.off()


# save ########################################################################
write_rds(int, here(resdir, "int.rds"), compress = "gz")

# end =========================================================================

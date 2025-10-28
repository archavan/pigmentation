# Arun Chavan
# Started: 2025-09-19

# setup =======================================================================
library(tidyverse)
library(Seurat)
library(here)
library(glue)

basedir <- here("species/mmus/PY_25-06245")
datadir <- here(
  "~/Library/CloudStorage/OneDrive-YaleUniversity/shared/PY+ARC/ghr/sc/25-06245"
)
genomedir <- here("species/mmus/genome")

# data ========================================================================
mm <- read_rds(here(datadir, "data/annotated_v1.rds"))
mm <- subset(mm, genotype == "WT")

modules <- read_tsv(here(genomedir, "gene-modules/functional-modules_gut.tsv"))

# subset to enterocytes =======================================================
ec <- subset(
  mm,
  cells = WhichCells(mm, idents = levels(mm)[grepl("Enterocyte", levels(mm))])
)
ec$ann_fine <- droplevels(ec$ann_fine)

barebones_seurat <- function(.obj) {
  cts <- GetAssayData(.obj, assay = "RNA", layer = "counts")
  md <- .obj[[]]
  new_obj <- CreateSeuratObject(counts = cts, meta.data = md)
  return(new_obj)
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

ec <- ec %>%
  barebones_seurat() %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 20) %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(
    resolution = 0.5,
    algorithm = 2,
    cluster.name = "ec_subcluster"
  ) %>%
  add_prefix_to_idents(.cluster.name = "ec_subcluster", .prefix = "ec")

# plot modules ================================================================

FetchData(
  ec,
  vars = c(str_split_1(modules$genes[[3]], "; "), "PC_1", "PC_3")
) |>
  pivot_longer(
    cols = -c(PC_1, PC_3),
    names_to = "gene",
    values_to = "expression"
  ) |>
  ggplot(aes(x = PC_1, y = PC_3, z = expression)) +
  stat_summary_hex(bins = 50, fun = "mean") +
  facet_wrap(vars(gene)) +
  scale_fill_viridis_c(direction = -1) +
  theme_minimal() +
  coord_fixed()

# end =========================================================================

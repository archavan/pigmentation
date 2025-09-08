# Arun Chavan
# Started: 2025-09-08

# background ==================================================================

# Here we:
# 1. Perform DEseq2 to find genes that are differentially expressed across the 5
#    tissues
# 2. Perform k-means clustering of differentially expressed genes to find
#    clusters of genes that show the same patterns of expression.

###############################################################################
# setup =======================================================================
###############################################################################

library(here)
library(DESeq2)
library(ggrepel)
library(tidyverse)
library(patchwork)
library(pvclust)
library(pheatmap)
library(conflicted)

conflicts_prefer(
  dplyr::select,
  dplyr::filter,
  dplyr::count,
  purrr::map,
  base::as.factor,
  dplyr::rename
)

basedir <- here("species/scla/bulk/PIG022/refseq_v1")
genomedir <- here("species/scla/genome/refseq_v1")
orthdir <- here("~/Dropbox/research/projects/genomes/orthofinder/2025-08-29")

resdir <- here(basedir, "results/cluster-genes/all-samples")
fs::dir_create(resdir)

source(here(basedir, "code/utilities_ggplot.R"))
source(here(basedir, "code/utilities_goseq.R"))

## data =======================================================================
include_samples <- c("S", "I1", "I2", "I3", "I4", "I5", "I6")

md <- read_csv(here(basedir, "data/samples_clean.csv")) |>
  filter(sample %in% include_samples) |>
  mutate(sample = fct(sample, include_samples))

tpm <- read_csv(here(basedir, "data/03_compiled/tpm_clean.csv")) |>
  select(all_of(c("gene_id", md$sample_id)))

cts <- read_csv(here(basedir, "data/03_compiled/counts_clean.csv")) |>
  select(all_of(c("gene_id", md$sample_id)))

## genomic data ===============================================================

genes <- read_tsv(here(genomedir, "results/genes/gene-display-names.tsv"))

mmgenes <- read_tsv(here(orthdir, "results/id-name-mapping/mmus.tsv")) |>
  select(gene_id, symbol) |>
  distinct()

scmm <- read_tsv(here(
  orthdir,
  "results/pairwise-orthologs/scla_vs_mmus.tsv"
)) |>
  left_join(mmgenes, by = c("mmus" = "gene_id"))

###############################################################################
# subset genes to include in k-means clustering ===============================
###############################################################################

# make sure that rows of colData correspond to the columns of countData
md <- md[match(md$sample_id, names(cts)[names(cts) != "gene_id"]), ]

# run DESeq
dds <- DESeqDataSetFromMatrix(
  countData = as.data.frame(cts),
  colData = as.data.frame(md),
  design = ~sample,
  tidy = TRUE
)

dds <- DESeq(dds, test = "LRT", reduced = ~1)

p_cutoff <- 0.1

# `results()` by defaults returns l2fcs between the last level and and the base
# level if contrasts are not specified. Regardless of the comparison, the p
# values for all genes identical because they are for the overall model, not for
# a specific comparison. Here we are interested only in finding genes that have
# some difference across all tissues, so we won't specify a contrast.

res_raw <- results(dds, alpha = p_cutoff)
res <- as.data.frame(res_raw) |>
  filter(!is.na(padj)) |>
  rownames_to_column("gene_id")

res_de <- res |>
  filter(padj < p_cutoff)

hist(res$pvalue)

include_genes <- res_de$gene_id

write_csv(res_de, here(resdir, "included-genes.csv"))

###############################################################################
# k-means clustering ==========================================================
###############################################################################

## prepare data for gene clustering ==========================================

tpm_de <- filter(tpm, gene_id %in% res_de$gene_id)

z_de <- tpm_de |>
  column_to_rownames("gene_id") |>
  # sqrt() |>
  t() |>
  scale() |>
  t() |>
  as.data.frame() |>
  rownames_to_column("gene_id") |>
  as_tibble()

z_de_mean <- z_de |>
  pivot_longer(
    cols = -gene_id,
    names_to = "sample_id",
    values_to = "z_score"
  ) |>
  left_join(select(md, sample_id, sample, rep), by = "sample_id") |>
  group_by(sample, gene_id) |>
  summarize(z_score = mean(z_score)) |>
  ungroup() |>
  pivot_wider(
    id_cols = "gene_id",
    names_from = "sample",
    values_from = "z_score"
  )

pheatmap::pheatmap(
  column_to_rownames(z_de_mean, "gene_id"),
  show_rownames = FALSE,
  cluster_cols = FALSE,
  filename = here(resdir, "included-genes.pdf"),
  width = 4,
  height = 8
)

## functions for clustering and plotting =======================================
cluster_genes <- function(k) {
  set.seed(100)
  k_out <- z_de_mean |>
    column_to_rownames("gene_id") |>
    kmeans(centers = k, iter.max = 25)
  rm(.Random.seed, envir = .GlobalEnv)
  return(k_out)
}

plot_clusters <- function(kclusts) {
  # order clusters
  cluster_order <- rownames(kclusts$centers)

  # get mean values of all clusters for plotting the average
  centers <- kclusts$centers |>
    as.data.frame() |>
    rownames_to_column("cluster") |>
    pivot_longer(
      cols = -cluster,
      names_to = "sample",
      values_to = "z_score"
    ) |>
    mutate(cluster = factor(cluster, levels = cluster_order)) |>
    mutate(sample = factor(sample, levels = include_samples))

  # number of genes in each cluster
  sizes <- kclusts$size |>
    enframe(name = "cluster", value = "size") |>
    mutate(cluster = factor(cluster, levels = cluster_order)) |>
    mutate(x = 5, y = -Inf)

  # gene level data for plotting
  plot_data <- z_de_mean |>
    pivot_longer(
      cols = -gene_id,
      values_to = "z_score",
      names_to = "sample"
    ) |>
    mutate(cluster = kclusts$cluster[.data[["gene_id"]]]) |>
    mutate(cluster = factor(cluster, levels = cluster_order)) |>
    mutate(sample = factor(sample, levels = include_samples))

  # plot
  ggplot(plot_data, aes(sample, z_score)) +
    geom_line(aes(group = gene_id), alpha = 0.2, color = "grey") +
    geom_line(
      data = centers,
      aes(sample, z_score, group = cluster),
      color = "red"
    ) +
    geom_point(data = centers, aes(sample, z_score), color = "red") +
    geom_text(
      data = sizes,
      aes(x = x, y = y, label = paste0(size, " genes")),
      hjust = 1,
      vjust = -1,
      size = 7 / .pt
    ) +
    facet_wrap(facets = vars(cluster), ncol = 4, labeller = label_both) +
    scale_x_discrete(guide = guide_axis(n.dodge = 1)) +
    labs(caption = paste0("k = ", length(unique(kclusts$cluster)))) +
    theme_bw_2()
}

## cluster at different k values to determine optimal k =======================
km <- tibble(k = seq(from = 1, to = 60, by = 1))

km <- mutate(km, km_out = map(k, cluster_genes))

# plot total of within-cluster SS as a measure of tightness of clusters
total_withinss <- km |>
  mutate(withinss = map_dbl(km_out, ~ .x$tot.withinss)) |>
  ggplot(aes(k, withinss)) +
  geom_point() +
  geom_line(aes(group = 1)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(y = "total of within-cluster sums of squares") +
  theme_bw_2()

ggsave(
  here(resdir, "total-within-cluster-ss.png"),
  total_withinss,
  width = 6,
  height = 6,
  units = "in",
  dpi = 300
)

# plot clusters at different k values
plot_k10 <- plot_clusters(km$km_out[km$k == 10][[1]])
plot_k12 <- plot_clusters(km$km_out[km$k == 12][[1]])
plot_k15 <- plot_clusters(km$km_out[km$k == 15][[1]])
plot_k20 <- plot_clusters(km$km_out[km$k == 20][[1]])

ggsave(here(resdir, "cluster_k10.png"), plot_k10, height = 6, width = 8)
ggsave(here(resdir, "cluster_k12.png"), plot_k12, height = 6, width = 8)
ggsave(here(resdir, "cluster_k15.png"), plot_k15, height = 7, width = 8)
ggsave(here(resdir, "cluster_k20.png"), plot_k20, height = 8, width = 8)

###############################################################################
## clusters at k = 12 =========================================================
###############################################################################

k12_out <- km$km_out[km$k == 12][[1]]

k12 <- enframe(k12_out$cluster, name = "gene_id", value = "cluster") |>
  left_join(select(genes, gene_id, symbol, description)) |>
  nest(.by = cluster) |>
  arrange(cluster)

# sort genes by average expression across all samples
tpm_mean_across_all_tissues <- tpm |>
  column_to_rownames("gene_id") |>
  apply(1, mean) |>
  enframe(name = "gene_id", value = "tpm") |>
  arrange(-tpm)

k12 <- k12 |>
  mutate(
    data = map(
      data,
      ~ left_join(.x, tpm_mean_across_all_tissues) |>
        relocate(tpm, .after = gene_id) |>
        rename(tpm_mean = tpm) |>
        arrange(-tpm_mean)
    )
  )


## write k12 output to file ====================================================
fs::dir_create(here(resdir, "k12/genes"))

write_rds(k12_out, here(resdir, "k12/k12_out.rds"))

walk2(
  k12$cluster,
  k12$data,
  ~ write_csv(.y, here(resdir, paste0("k12/genes/cluster-", .x, ".csv")))
)

# plot individual genes to test patterns ======================================
tpm_long <- tpm |>
  pivot_longer(-gene_id, names_to = "sample_id", values_to = "tpm") |>
  left_join(md)

tpm_mean_long <- tpm_long |>
  group_by(sample, gene_id) |>
  summarize(tpm = mean(tpm)) |>
  ungroup()

plot_single_gene <- function(.gene_id) {
  dd <- filter(tpm_long, gene_id == .gene_id)
  dd_mean <- filter(tpm_mean_long, gene_id == .gene_id)
  display_name <- paste0(
    "[",
    genes$gene_id[genes$gene_id == .gene_id],
    "] ",
    genes$symbol[genes$gene_id == .gene_id],
    "\n",
    str_trunc(genes$description[genes$gene_id == .gene_id], 72)
  )
  ggplot() +
    geom_jitter(data = dd, aes(sample, tpm), width = 0.1, color = "grey") +
    geom_line(data = dd_mean, aes(sample, tpm, group = gene_id)) +
    geom_point(data = dd_mean, aes(sample, tpm)) +
    scale_x_discrete(guide = guide_axis(n.dodge = 1)) +
    labs(title = display_name) +
    theme_bw_2() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 6),
      axis.text = element_text(size = 8),
      plot.title = element_text(size = 7)
    )
}

plot_single_gene("ncbi.120326176")
genes |>
  filter(str_detect(display_name, regex("scpep", TRUE))) |>
  pull(gene_id) |>
  map(plot_single_gene) |>
  wrap_plots()

###############################################################################
# GO enrichment ===============================================================
###############################################################################

k12 <- k12 |>
  mutate(
    go_bg_de = map(
      data,
      ~ run_goseq(.x$gene_id, .bg = tpm_de$gene_id) |>
        filter(p_adj < 0.05)
    )
  ) |>
  mutate(
    go_bg_all = map(
      data,
      ~ run_goseq(.x$gene_id, .bg = tpm$gene_id) |>
        filter(p_adj < 0.05)
    )
  )

fs::dir_create(here(resdir, "k12/go/bg_de"))
fs::dir_create(here(resdir, "k12/go/bg_all"))
walk2(
  k12$cluster,
  k12$go_bg_de,
  ~ write_tsv(.y, here(resdir, paste0("k12/go/bg_de/cluster-", .x, ".tsv")))
)
walk2(
  k12$cluster,
  k12$go_bg_all,
  ~ write_tsv(.y, here(resdir, paste0("k12/go/bg_all/cluster-", .x, ".tsv")))
)

## generate input for running metascape =======================================
# get orthologs
get_orthologs <- function(.gene_ids) {
  scmm |>
    filter(scla %in% .gene_ids) |>
    pull(symbol) |>
    unique()
}

generate_metascape_input <- function(.geneset_names, .genesets) {
  meta <- as.data.frame(
    matrix(
      nrow = max(lengths(.genesets)),
      ncol = length(.geneset_names)
    )
  )

  names(meta) <- .geneset_names

  for (i in seq_along(.geneset_names)) {
    meta[[.geneset_names[i]]] <- c(
      .genesets[[i]],
      rep("", nrow(meta) - length(.genesets[[i]]))
    )
  }

  return(meta)
}

k12 <- k12 |>
  mutate(mm_orth = map(data, ~ get_orthologs(.x$gene_id)))

metascape <- generate_metascape_input(
  .geneset_names = paste0("cluster-", k12$cluster),
  .genesets = k12$mm_orth
)

fs::dir_create(here(resdir, "k12/metascape"))
write_csv(metascape, here(resdir, "k12/metascape/input.csv"))

###############################################################################
# heatmaps ====================================================================
###############################################################################
fs::dir_create(here(resdir, "k12/heatmaps"))

save_heatmap <- function(.clust, .w, .h) {
  clust_genes <- k12$data[[which(k12$cluster == .clust)]]$gene_id
  x <- filter(z_de_mean, gene_id %in% clust_genes)
  x <- column_to_rownames(x, "gene_id")
  rownames(x) <- genes$display_name[match(rownames(x), genes$gene_id)]
  pheatmap(
    x,
    cluster_cols = FALSE,
    filename = here(resdir, "k12/heatmaps", paste0("cluster-", .clust, ".pdf")),
    width = .w,
    height = .h
  )
}

save_heatmap(1, 10, 60)
save_heatmap(2, 9.75, 40)
save_heatmap(3, 9.75, 35)
save_heatmap(4, 9.75, 35)
save_heatmap(5, 10, 12)
save_heatmap(6, 9.5, 27)
save_heatmap(7, 10, 32)
save_heatmap(8, 10, 25)
save_heatmap(9, 10.75, 52)
save_heatmap(10, 9.5, 22)
save_heatmap(11, 10, 58)
save_heatmap(12, 9.5, 33)

# end =========================================================================

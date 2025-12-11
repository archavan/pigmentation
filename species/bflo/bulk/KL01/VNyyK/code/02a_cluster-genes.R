# Arun Chavan
# Started: 2023-12-31

# background ==================================================================

# Here we:
# 1. Perform DEseq2 to find genes that are differentially expressed across the 5
#    tissues
# 2. Perform k-means clustering of differentially expressed genes to find 
#    clusters of genes that show the same patterns of expression.

# setup =======================================================================
library(here)
library(DESeq2)
library(ggrepel)
library(tidyverse)
library(patchwork)

basedir <- here("species/bflo/bulk/KL01/")
genomedir <- here("species/bflo/genome/VNyyK")
resdir <- here(basedir, "results/02a_cluster-genes")
fs::dir_create(resdir)

source(here(basedir, "code/utilities_ggplot.R"))
source(here(basedir, "code/utilities_goseq.R"))

# data ========================================================================
tpm <- read_tsv(here(basedir, "data/03_compiled/tpm_gut.tsv"))
counts <- read_tsv(here(basedir, "data/03_compiled/counts_gut.tsv"))
samples <- read_tsv(here(basedir, "data/00_metadata/samples_gut.tsv"))

genes <- read_tsv(here(genomedir, "results/genes/gene-display-names.tsv"))
orth <- read_csv(here(genomedir, "results/emapper/processed/emapper_orthologs_mouse.csv"))

# DESeq =======================================================================
# prepare data
tissue_order <- c("mdiv", "mgut", "hgut")

samples <- mutate(samples, tissue = factor(tissue, levels = tissue_order))

# make sure that rows of colData correspond to the columns of countData
samples <- samples[match(names(counts)[names(counts) != "gene_id"], 
                         samples$sample_name), ]

# run DESeq
dds <- DESeqDataSetFromMatrix(countData = as.data.frame(counts),
                              colData = as.data.frame(samples),
                              design = ~ tissue,
                              tidy = TRUE)

dds <- DESeq(dds, test = "LRT", reduced = ~ 1)

p_cutoff <- 0.05

# `results()` by defaults returns l2fcs between the last level and and the base
# level if contrasts are not specified. Regardless of the comparison, the p
# values for all genes identical because they are for the overall model, not for
# a specific comparison. Here we are interested only in finding genes that have
# some difference across all tissues, so we won't specify a contrast.

res_raw <- results(dds, alpha = p_cutoff)
res <- as.data.frame(res_raw) %>% 
  filter(!is.na(padj)) %>% 
  rownames_to_column("gene_id")

# prepare data for gene clustering ============================================
# extract genes that are different across tissues
res_de <- res %>% 
  as_tibble() %>% 
  filter(padj < p_cutoff)

write_csv(res_de, here(resdir, "de-genes-across-all-tissues.csv"))

# subset tpm to DE genes
tpm_de <- filter(tpm, gene_id %in% res_de$gene_id)

# sqrt transform and scale tpms
z_de <- tpm_de %>% 
  column_to_rownames("gene_id") %>% 
  sqrt() %>% 
  t() %>% 
  scale() %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("gene_id") %>% 
  as_tibble()

z_de_mean <- z_de %>% 
  pivot_longer(cols = -gene_id, names_to = "sample_name", values_to = "z_score") %>% 
  left_join(select(samples, sample_name, tissue, rep_id)) %>% 
  group_by(tissue, gene_id) %>% 
  summarize(z_score = mean(z_score)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = "gene_id", names_from = "tissue", values_from = "z_score")

# functions for clustering and plotting =======================================
cluster_genes <- function(k) {
  set.seed(100)
  k_out <- z_de_mean %>% 
    column_to_rownames("gene_id") %>% 
    kmeans(centers = k, iter.max = 25)
  rm(.Random.seed, envir = .GlobalEnv)
  return(k_out)
}

plot_clusters <- function(kclusts) {
  
  # order clusters 
  cluster_order <- rownames(kclusts$centers)
  
  # get mean values of all clusters for plotting the average
  centers <- kclusts$centers %>% 
    as.data.frame() %>% 
    rownames_to_column("cluster") %>% 
    pivot_longer(cols = -cluster, names_to = "tissue", values_to = "z_score") %>% 
    mutate(cluster = factor(cluster, levels = cluster_order)) %>% 
    mutate(tissue = factor(tissue, levels = tissue_order))
  
  # number of genes in each cluster
  sizes <- kclusts$size %>% 
    enframe(name = "cluster", value = "size") %>% 
    mutate(cluster = factor(cluster, levels = cluster_order)) %>% 
    mutate(x = 3, y = -Inf)
  
  # gene level data for plotting
  plot_data <- z_de_mean %>% 
    pivot_longer(cols = -gene_id, values_to = "z_score", names_to = "tissue") %>% 
    mutate(cluster = kclusts$cluster[.$gene_id]) %>% 
    mutate(cluster = factor(cluster, levels = cluster_order))  %>% 
    mutate(tissue = factor(tissue, levels = tissue_order))
  
  # plot
  ggplot(plot_data, aes(tissue, z_score)) +
    geom_line(aes(group = gene_id), alpha = 0.2, color = "grey") +
    geom_line(data = centers, aes(tissue, z_score, group = cluster), color = "red") +
    geom_point(data = centers, aes(tissue, z_score), color = "red") +
    geom_text(data = sizes, aes(x = x, y = y, label = paste0(size, " genes")),
              hjust = 1, vjust = -1, size = 7/.pt) +
    facet_wrap(facets = vars(cluster), ncol = 4, labeller = label_both) +
    scale_x_discrete(guide = guide_axis(n.dodge = 1)) +
    labs(caption = paste0("k = ", length(unique(kclusts$cluster)))) +
    theme_bw_2()
}

# cluster at different k values to determine optimal k ========================
km <- tibble(k = seq(from = 1, to = 60, by = 1))

km <- mutate(km, km_out = map(k, cluster_genes))

# plot total of within-cluster SS as a measure of tightness of clusters
total_withinss <- km %>% 
  mutate(withinss = map_dbl(km_out, ~ .x$tot.withinss)) %>% 
  ggplot(aes(k, withinss)) +
  geom_point() +
  geom_line(aes(group = 1)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(y = "total of within-cluster sums of squares") +
  theme_bw_2()

ggsave(here(resdir, "total-within-cluster-ss.png"),
       total_withinss, width = 6, height = 6, units = "in", dpi = 300)
       
# plot clusters at different k values
plot_k6 <- plot_clusters(km$km_out[km$k == 6][[1]])
plot_k8 <- plot_clusters(km$km_out[km$k == 8][[1]])
plot_k10 <- plot_clusters(km$km_out[km$k == 10][[1]])

ggsave(here(resdir, "cluster_k6.png"),  plot_k6,  height = 4, width = 6, dpi = 300)
ggsave(here(resdir, "cluster_k8.png"),  plot_k8,  height = 4, width = 6, dpi = 300)
ggsave(here(resdir, "cluster_k10.png"), plot_k10, height = 5, width = 6, dpi = 300)

# clusters at k = 6 ==========================================================
k6_out <- km$km_out[km$k == 6][[1]]
# if restarting analysis, read saved k6_out object
# k6_out <- read_rds(here(resdir, "k6/k6_out.rds"))

k6 <- enframe(k6_out$cluster, name = "gene_id", value = "cluster") %>% 
  left_join(select(genes, gene_id, symbol, description)) %>% 
  nest(.by = cluster) %>% 
  arrange(cluster)

# sort genes by average expression across all samples
tpm_mean_across_all_tissues <- tpm %>% 
  column_to_rownames("gene_id") %>% 
  apply(., 1, mean) %>% 
  enframe(name = "gene_id", value = "tpm") %>% 
  arrange(-tpm)

k6 <- k6 %>% 
  mutate(data = map(data, ~ left_join(.x, tpm_mean_across_all_tissues) %>% 
                      relocate(tpm, .after = gene_id) %>% 
                      rename(tpm_mean = tpm) %>% 
                      arrange(-tpm_mean)))

# write k6 output to file ====================================================
fs::dir_create(here(resdir, "k6"))
fs::dir_create(here(resdir, "k6/genes"))

write_rds(k6_out, here(resdir, "k6/k6_out.rds"))

walk2(k6$cluster, 
      k6$data, 
      ~ write_csv(.y, here(resdir, paste0("k6/genes/cluster-", .x, ".csv" ))))

# plot individual genes to test patterns ======================================
tpm_long <- tpm %>% 
  pivot_longer(-gene_id, names_to = "sample_name", values_to = "tpm") %>% 
  separate(sample_name, c("tissue", "rep_id"), sep = "_", remove = FALSE)

tpm_mean_long <- tpm_long %>% 
  separate(sample_name, c("tissue", NA), sep = "_") %>% 
  group_by(tissue, gene_id) %>% 
  summarize(tpm = mean(tpm)) %>% 
  ungroup()

plot_single_gene <- function(.gene_id) {
  dd <- filter(tpm_long, gene_id == .gene_id) %>% 
    mutate(tissue = factor(tissue, levels = tissue_order))
  dd_mean <- filter(tpm_mean_long, gene_id == .gene_id) %>% 
    mutate(tissue = factor(tissue, levels = tissue_order))
  display_name <- paste0("[", 
                         genes$gene_id[genes$gene_id == .gene_id], 
                         "] ", 
                         genes$symbol[genes$gene_id == .gene_id], 
                         "\n", 
                         str_trunc(genes$description[genes$gene_id == .gene_id], 
                                   72))
  ggplot() +
    geom_jitter(data = dd, aes(tissue, tpm), width = 0.1, color = "grey") +
    geom_line(data = dd_mean, aes(tissue, tpm, group = gene_id)) +
    geom_point(data = dd_mean, aes(tissue, tpm)) +
    scale_x_discrete(guide = guide_axis(n.dodge = 1)) +
    labs(title = display_name) +
    theme_bw_2() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 6),
          axis.text = element_text(size = 8),
          plot.title = element_text(size = 7)) 
}

plot_single_gene("ncbi:118432570")

# GO enrichment ===============================================================

k6 <- k6 %>% 
  mutate(go_bg_de = map(data, ~ run_goseq(.x$gene_id, .bg = tpm_de$gene_id) %>% 
                          filter(p_adj < 0.05))) %>% 
  mutate(go_bg_all = map(data, ~ run_goseq(.x$gene_id, .bg = tpm$gene_id) %>% 
                           filter(p_adj < 0.05)))

fs::dir_create(here(resdir, "k6/go/bg_de"))
fs::dir_create(here(resdir, "k6/go/bg_all"))
walk2(k6$cluster, 
      k6$go_bg_de,
      ~ write_tsv(.y, here(resdir, paste0("k6/go/bg_de/cluster-", .x, ".tsv"))))
walk2(k6$cluster, 
      k6$go_bg_all,
      ~ write_tsv(.y, here(resdir, paste0("k6/go/bg_all/cluster-", .x, ".tsv"))))

# generate input for running metascape ========================================
# get orthologs
get_orthologs <- function(.gene_ids, types = c("one2one")) {
  orth %>% 
    filter(bf_gene_id %in% .gene_ids) %>% 
    filter(orth_type %in% types) %>% 
    pull(mm_gene_id) %>% 
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
    meta[[.geneset_names[i]]] <- c(.genesets[[i]],
                                   rep("", nrow(meta) - length(.genesets[[i]])))
  }
  
  return(meta)
}

k6 <- k6 %>% 
  mutate(one2one = map(data, ~ get_orthologs(.x$gene_id)))

metascape <- generate_metascape_input(.geneset_names = paste0("cluster-", k6$cluster),
                                      .genesets = k6$one2one)

fs::dir_create(here(resdir, "k6/metascape"))
write_csv(metascape, here(resdir, "k6/metascape/input.csv"))

# end =========================================================================





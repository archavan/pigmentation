# Arun Chavan
# Started: 2023-11-08

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
library(pvclust)

basedir <- here("species/crob/bulk/ARCA123/KY21")
genomedir <- here("species/crob/genome/KY21")
fs::dir_create(here(basedir, "results/06a_cluster-genes_all-samples"))
resdir <- here(basedir, "results/06a_cluster-genes_all-samples")

source(here(basedir, "code/utilities_ggplot.R"))
source(here(basedir, "code/utilities_goseq.R"))

# data ========================================================================
tpm <- read_csv(here(basedir, "data/03_compiled/tpm_no-testis.csv"))
counts <- read_csv(here(basedir, "data/03_compiled/counts_no-testis.csv"))
samples <- read_csv(here(basedir, "data/00_metadata/samples.csv"))

genes <- read_tsv(here(genomedir, "results/genes/gene-display-names.tsv"))
orth <- read_tsv(here(genomedir, "results/ghost/gene-info_ghostdb_parsed.tsv"))

# DESeq =======================================================================
# prepare data
tissue_order <- c("stm.1.s", "int.1.w", "int.2.y", "int.3.w", "int.4.y")

samples <- mutate(samples, tissue = factor(tissue, levels = tissue_order))

# make sure that rows of colData correspond to the columns of countData
samples <- samples[match(samples$sample_name, 
                         names(counts)[names(counts) != "gene_id"]), ]

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
    mutate(x = 5, y = -Inf)
  
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
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
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
plot_k10 <- plot_clusters(km$km_out[km$k == 10][[1]])
plot_k15 <- plot_clusters(km$km_out[km$k == 15][[1]])
plot_k20 <- plot_clusters(km$km_out[km$k == 20][[1]])

ggsave(here(resdir, "cluster_k10.png"), plot_k10, height = 6, width = 8, dpi = 300)
ggsave(here(resdir, "cluster_k15.png"), plot_k15, height = 7, width = 8, dpi = 300)
ggsave(here(resdir, "cluster_k20.png"), plot_k20, height = 8, width = 8, dpi = 300)

# clusters at k = 20 ==========================================================
k20_out <- km$km_out[km$k == 20][[1]]
# if restarting analysis, read saved k20_out object
# k20_out <- read_rds(here(resdir, "k20/k20_out.rds"))

k20 <- enframe(k20_out$cluster, name = "gene_id", value = "cluster") %>% 
  left_join(select(genes, gene_id, symbol, description)) %>% 
  nest(.by = cluster) %>% 
  arrange(cluster)

# sort genes by average expression across all samples
tpm_mean_across_all_tissues <- tpm %>% 
  column_to_rownames("gene_id") %>% 
  apply(., 1, mean) %>% 
  enframe(name = "gene_id", value = "tpm") %>% 
  arrange(-tpm)

k20 <- k20 %>% 
  mutate(data = map(data, ~ left_join(.x, tpm_mean_across_all_tissues) %>% 
                      relocate(tpm, .after = gene_id) %>% 
                      rename(tpm_mean = tpm) %>% 
                      arrange(-tpm_mean)))

# similarities of clusters ====================================================

# Here we perform clustering of clusters to get a big picture of how the gene
# clusters are related to each other, whether there are some clusters that are
# so similar to each othe that they can be merged into a single cluster, etc.

pv <- pvclust(t(k20_out$centers), method.hclust = "complete")

fs::dir_create(here(resdir, "k20"))
png(here(resdir, "k20/pvlust-on-cluster-centers.png"), 
    10, 6, res = 300, units = "in")
plot(pv)
dev.off()

# write k20 output to file ====================================================
fs::dir_create(here(resdir, "k20/genes"))

write_rds(k20_out, here(resdir, "k20/k20_out.rds"))

walk2(k20$cluster, 
      k20$data, 
      ~ write_csv(.y, here(resdir, paste0("k20/genes/cluster-", .x, ".csv" ))))

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
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    labs(title = display_name) +
    theme_bw_2() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 6),
          axis.text = element_text(size = 8),
          plot.title = element_text(size = 7)) 
}

 plot_single_gene("KY21.Chr1.1278")

# GO enrichment ===============================================================

k20 <- k20 %>% 
  mutate(go_bg_de = map(data, ~ run_goseq(.x$gene_id, .bg = tpm_de$gene_id) %>% 
                          filter(p_adj < 0.05))) %>% 
  mutate(go_bg_all = map(data, ~ run_goseq(.x$gene_id, .bg = tpm$gene_id) %>% 
                           filter(p_adj < 0.05)))

fs::dir_create(here(resdir, "k20/go/bg_de"))
fs::dir_create(here(resdir, "k20/go/bg_all"))
walk2(k20$cluster, 
      k20$go_bg_de,
      ~ write_tsv(.y, here(resdir, paste0("k20/go/bg_de/cluster-", .x, ".tsv"))))
walk2(k20$cluster, 
      k20$go_bg_all,
      ~ write_tsv(.y, here(resdir, paste0("k20/go/bg_all/cluster-", .x, ".tsv"))))

# generate input for running metascape ========================================
# get orthologs
get_orthologs <- function(.gene_ids) {
  orth %>% 
    filter(gene_id %in% .gene_ids) %>% 
    pull(hsap_homolog_symbol) %>% 
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

k20 <- k20 %>% 
  mutate(hsap_orth = map(data, ~ get_orthologs(.x$gene_id)))

metascape <- generate_metascape_input(.geneset_names = paste0("cluster-", k20$cluster),
                                      .genesets = k20$hsap_orth)

# clusters 12 and 16 is mostly testis genes that dominate the metascape results,
# so we will remove those columns and save the input file.

metascape <- metascape[!names(metascape) %in% c("cluster-12", "cluster-16")]

fs::dir_create(here(resdir, "k20/metascape"))
write_csv(metascape, here(resdir, "k20/metascape/input.csv"))

# end =========================================================================





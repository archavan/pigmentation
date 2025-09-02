# Arun Chavan
# Started: 2024-11-13

# background ==================================================================

# Here, we compare the bulk RNA-seq data with the single-cell RNA-seq data.
# Specifically, we test the expression patterns in the bulk RNA-seq data of the
# marker genes of clusters in single-cell RNA seq data to get a sense of the
# spatial distribution of the cell types.

# setup =======================================================================
library(tidyverse)
library(here)
library(ggrepel)
library(Seurat)
library(glue)
library(patchwork)

basedir <- here("species/crob/bulk/ARCA123/KY21")
genomedir <- here("species/crob/genome/KY21")
resdir <- here(basedir, "results/11a_single-cell-comparison")
fs::dir_create(resdir)

scdir <- here("species/crob/10x/run-02/KY21/")

source(here(basedir, "code/utilities_ggplot.R"))
source(here(scdir, "code/utilities_seurat-plots.R"))

# data ========================================================================
tpm <- read_csv(here(basedir, "data/03_compiled/tpm_no-testis.csv"))
samples <- read_csv(here(basedir, "data/00_metadata/samples.csv"))

genes <- read_tsv(here(genomedir, "results/genes/gene-display-names.tsv"))

sc <- read_rds(here(scdir, "data/processed/03_clustered/clustered_v1.rds"))
markers <- read_csv(here(scdir, "results/v1/01_clustering/markers_filtered.csv"))
tau <- read_csv(here(scdir, "results/v1/01_clustering/tau.csv"))

# for plotting ================================================================
tissue_order <- c("stm.1.s", "int.1.w", "int.2.y", "int.3.w", "int.4.y")
cluster_order <- read_csv(here(scdir, 
                               "results/v1/01_clustering/cluster-order.csv")) %>% 
  pull(cluster_v1)

tissue_origin <- read_csv(here(scdir, "results/v1/01_clustering/tissue-of-origin.csv")) %>% 
  select(cluster_v1, tissue) %>% 
  deframe()

display_names <- deframe(select(genes, gene_id, display_name))

plot_dotplot2 <- function(.features, .title) {
  plot_dotplot(.obj = sc, 
               .features = .features,
               .title = .title, 
               .cluster_order = cluster_order, 
               .display_names = display_names, 
               .tissue_of_origin = tissue_origin)
}

# find cluster specific markers to test =======================================
markers <- markers %>% 
  rename(gene_id = gene) %>% 
  left_join(tau, by = "gene_id") %>% 
  mutate(specificity = tau * abs(avg_log2FC)) %>% 
  nest(.by = cluster) %>% 
  mutate(spec_markers = map(data, ~ slice_max(.x, specificity, n = 10)))

markers <- markers %>% 
  mutate(dp = map2(spec_markers, 
                   cluster, 
                   ~ plot_dotplot2(.x$gene_id, paste0("cluster ", .y))))

# plot the expression of markers along AP axis with bulk data =================
tpm_long <- tpm %>% 
  pivot_longer(cols = -gene_id, names_to = "sample_name", values_to = "tpm") %>% 
  left_join(samples) %>% 
  mutate(tissue = factor(tissue, tissue_order))

plot_genes_faceted <- function(.gene_ids, .title) {
  
  plot_data <- tpm_long %>% 
    filter(gene_id %in% .gene_ids) %>% 
    left_join(genes) %>% 
    mutate(display_name = paste0("[", gene_id, "] ", symbol, "\n", 
                                 str_trunc(description, 72)))
  
  # order genes by average expression
  gene_order <- plot_data %>% 
    summarise(tpm = mean(tpm), .by = display_name) %>% 
    arrange(-tpm) %>% 
    pull(display_name) 
  
  plot_data <- mutate(plot_data, display_name = factor(display_name, gene_order))
  plot_data <- arrange(plot_data, display_name)
  
  ggplot(plot_data, aes(x = tissue)) +
    geom_jitter(aes(y = tpm), width = 0.1, color = "grey") +
    stat_summary(aes(y = tpm, group = 1), fun = mean, geom = "line") +
    stat_summary(aes(y = tpm), fun = mean, geom = "point") +
    facet_wrap(vars(display_name), 
               scales = "free_y",
               nrow = 2, ncol = 5) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    theme_bw_2() +
    theme(strip.text = element_text(size = 6, hjust = 0)) +
    ggtitle(.title)
  
}

markers <- markers %>% 
  mutate(plot_bulk = map2(cluster, 
                          spec_markers,
                          ~ plot_genes_faceted(.y$gene_id,
                                               paste0("cluster ", .x))))

# put markers plots and bulk plots together for each cluster ==================
markers <- markers %>% 
  mutate(patch = map2(dp, 
                      plot_bulk, 
                      ~ wrap_plots(.x[[1]],
                                   .y,
                                   nrow = 2,
                                   heights = c(3, 2))))

fs::dir_create(here(resdir, "marker-plots-and-expression-pattern"))
walk2(markers$patch,
      markers$cluster,
      ~ ggsave(here(resdir, 
                    "marker-plots-and-expression-pattern", 
                    paste0("cluster-", .y, ".pdf")),
               .x,
               width = 13, height = 8, 
               device = cairo_pdf))

# end =========================================================================



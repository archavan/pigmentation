# Arun Chavan
# Started: 2023-12-01

# background ==================================================================

# Here we perform DEseq2 to find DE genes between the yellow and white segments
# of the intestine

# setup =======================================================================
library(here)
library(DESeq2)
library(ggrepel)
library(tidyverse)
library(ggforce)

basedir <- here("species/crob/bulk/ARCA123/KY21")
genomedir <- here("species/crob/genome/KY21")
fs::dir_create(here(basedir, "results/08a_yellow-vs-white"))
resdir <- here(basedir, "results/08a_yellow-vs-white")

source(here(basedir, "code/utilities_ggplot.R"))
source(here(basedir, "code/utilities_goseq.R"))

# data ========================================================================
tpm <- read_csv(here(basedir, "data/03_compiled/tpm_no-testis.csv"))
counts <- read_csv(here(basedir, "data/03_compiled/counts_no-testis.csv"))
samples <- read_csv(here(basedir, "data/00_metadata/samples.csv"))

genes <- read_tsv(here(genomedir, "results/genes/gene-display-names.tsv"))
orth <- read_tsv(here(genomedir, "results/ghost/gene-info_ghostdb_parsed.tsv"))

# select intestine samples
tpm <- select(tpm, -starts_with("stm"))
counts <- select(counts, -starts_with("stm"))
samples <- filter(samples, tissue != "stm.1.s")

# DESeq =======================================================================
# prepare data
tissue_order <- c("int.1.w", "int.2.y", "int.3.w", "int.4.y")
color_order <- c("white", "yellow")

samples <- mutate(samples, 
                  tissue = factor(tissue, levels = tissue_order),
                  colour = factor(colour, levels = color_order))

# make sure that rows of colData correspond to the columns of countData
samples <- samples[match(samples$sample_name, 
                         names(counts)[names(counts) != "gene_id"]), ]

# run DESeq
dds <- DESeqDataSetFromMatrix(countData = as.data.frame(counts),
                              colData = as.data.frame(samples),
                              design = ~ colour,
                              tidy = TRUE)

dds <- DESeq(dds, test = "LRT", reduced = ~ 1)

p_cutoff <- 0.05

# extract results =============================================================
# log fold change shrinkage for visualization
res_shrunken_raw <- lfcShrink(dds, 
                              coef = "colour_yellow_vs_white", 
                              type = "apeglm")

res_shrunken <- as.data.frame(res_shrunken_raw) %>% 
  filter(!is.na(padj)) %>% 
  rownames_to_column("gene_id") %>% 
  rename(l2fc_yellow_vs_white = log2FoldChange)

write_csv(res_shrunken, here(resdir, "deseq_out_shrunken.csv"))

# results
res_raw <- results(dds, 
                   contrast = c("colour", "yellow", "white"), 
                   alpha = p_cutoff)

res <- as.data.frame(res_raw) %>% 
  filter(!is.na(padj)) %>% 
  rownames_to_column("gene_id") %>% 
  rename(l2fc_yellow_vs_white = log2FoldChange)

write_csv(res, here(resdir, "deseq_out.csv"))

res_de <- res %>% 
  as_tibble() %>% 
  filter(padj < p_cutoff)

# subset tpm to DE genes
tpm_de <- filter(tpm, gene_id %in% res_de$gene_id)

# diagnostic plots ============================================================
# p value distribution 
p_val_dist <- ggplot(res, aes(pvalue)) +
  geom_histogram(binwidth = 0.05, boundary = 0, color = "white") +
  theme_bw_2()

ggsave(here(resdir, "pvalue-distribution.png"), p_val_dist,
       width = 4, height = 3, units = "in", dpi = 300)

# MA plot
png(here(resdir, "MA-plots.png"), width = 8, height = 4, units = "in", res = 300)
par(mfrow = c(1, 2))
plotMA(res_raw, main = "raw log2 fold changes")
plotMA(res_shrunken_raw, main = "shrunken log2 fold changes")
dev.off()

# lists of genes ==============================================================
res_de <- left_join(res_de, genes, by = "gene_id")

yellow_up <- filter(res_de, l2fc_yellow_vs_white > 0) %>% 
  arrange(-l2fc_yellow_vs_white)
white_up  <- filter(res_de, l2fc_yellow_vs_white < 0) %>% 
  arrange(l2fc_yellow_vs_white)

fs::dir_create(here(resdir, "gene-lists"))
write_csv(arrange(res_de, padj), here(resdir, "gene-lists/de.csv"))
write_csv(yellow_up, here(resdir, "gene-lists/yellow-up.csv"))
write_csv(white_up, here(resdir, "gene-lists/white-up.csv"))

# volcano plot ================================================================
volcano_data <- as_tibble(res) %>% 
  left_join(genes, by = "gene_id")

volcano_labels <- volcano_data  %>% 
  filter(padj < 0.05) %>% 
  mutate(change = ifelse(l2fc_yellow_vs_white > 0, "yellow_up", "white_up")) %>% 
  nest(.by = change) %>% 
  mutate(data = map(data, 
                    ~ .x %>% 
                      arrange(-abs(l2fc_yellow_vs_white)) %>% 
                      filter(!is.na(symbol)) %>% 
                      slice_head(n = 50))) %>% 
  unnest(data)

volcano_plot <- ggplot(volcano_data, 
       aes(l2fc_yellow_vs_white, -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.8, size = 1.2) +
  geom_text_repel(data = volcano_labels, 
                  mapping = aes(label = symbol),
                  color = "blue",
                  size = 8/.pt) +
  scale_color_manual(values = c("grey90", "grey40")) +
  labs(x = "log2 fold change (yellow vs white)") +
  theme_bw_2()

ggsave(here(resdir, "volcano-plot-labelled.png"), volcano_plot,
       width = 8, height = 6, units = "in", dpi = 300)

# MA plot, labelled ===========================================================
MA_data <- as_tibble(res_shrunken) %>% 
  left_join(genes, by = "gene_id")

MA_labels <- MA_data  %>% 
  filter(padj < 0.05, abs(l2fc_yellow_vs_white) > 1) 

MA_plot_labeled <- ggplot(MA_data, aes(log10(baseMean), l2fc_yellow_vs_white)) +
  geom_hline(yintercept = 0, linetype = 2, size = 0.5) +
  geom_point(aes(color = padj < 0.05), alpha = 1, size = 1.2) +
  geom_text_repel(data = MA_labels, 
                  mapping = aes(label = symbol),
                  color = "blue", size = 8/.pt) +
  scale_color_manual(values = c("grey90", "grey40")) +
  labs(x = "Average Expression (log10)", 
       y = "log2 fold change (yellow vs white)") +
  theme_bw_2()

ggsave(here(resdir, "MA-plot-labelled.png"), MA_plot_labeled,
       width = 8, height = 6, units = "in", dpi = 300)


# plot individual DE genes ====================================================
tpm_long <- tpm %>% 
  pivot_longer(cols = -gene_id, names_to = "sample_name", values_to = "tpm") %>% 
  left_join(samples)

calculate_num_pages <- function(n_genes, genes_per_page) {
  ifelse(n_genes %% genes_per_page == 0, 
         n_genes %/% genes_per_page, 
         (n_genes %/% genes_per_page) + 1) %>% 
    round() %>% 
    as.integer()
}

plot_genes_paginated <- function(.gene_ids, .pg_num) {
  
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
    facet_wrap_paginate(vars(display_name), 
                        scales = "free_y",
                        nrow = 8, ncol = 5, page = .pg_num) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    theme_bw_2() +
    theme(strip.text = element_text(size = 6, hjust = 0))
}

cairo_pdf(filename = here(resdir, "gene-lists/white-up.pdf"), 
          width = 8.25, height = 11,
          onefile = TRUE)
for (i in seq_len(calculate_num_pages(length(white_up$gene_id), 40))) {
  print(plot_genes_paginated(white_up$gene_id, i))
}
dev.off()

cairo_pdf(filename = here(resdir, "gene-lists/yellow-up.pdf"), 
          width = 8.25, height = 11,
          onefile = TRUE)
for (i in seq_len(calculate_num_pages(length(yellow_up$gene_id), 40))) {
  print(plot_genes_paginated(yellow_up$gene_id, i))
}
dev.off()

# GO enrichment ===============================================================
yellow_up_go <- run_goseq(yellow_up$gene_id, tpm$gene_id) %>% 
  filter(p_adj < 0.05)
white_up_go <- run_goseq(white_up$gene_id, tpm$gene_id) %>% 
  filter(p_adj < 0.05)

fs::dir_create(here(resdir, "go"))
write_tsv(yellow_up_go, here(resdir, "go/yellow_up.tsv"))
write_tsv(white_up_go,  here(resdir, "go/white_up.tsv"))

# generate input for metascape ================================================
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

metascape <- list(white_up$gene_id, yellow_up$gene_id) %>% 
  map(get_orthologs) %>% 
  generate_metascape_input(.genesets = ., 
                           .geneset_names = c("white_up", "yellow_up"))

fs::dir_create(here(resdir, "metascape"))
write_csv(metascape, here(resdir, "metascape/input.csv"))

# end =========================================================================


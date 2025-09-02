# Arun Chavan
# Started: 2023-12-07

# background ==================================================================

# Here we perform DESeq2 between the two yellow segments and between the two
# white segment. This is important because even though they appear the same
# colour, the two segments can have different transcriptional signatures.

# setup =======================================================================
library(here)
library(DESeq2)
library(ggrepel)
library(tidyverse)
library(ggforce)

basedir <- here("species/crob/bulk/ARCA123/KY21")
genomedir <- here("species/crob/genome/KY21")
resdir <- here(basedir, "results/08b_segments-of-same-colour")
fs::dir_create(resdir)
fs::dir_create(here(resdir, "yellow"))
fs::dir_create(here(resdir, "white"))

source(here(basedir, "code/utilities_ggplot.R"))
source(here(basedir, "code/utilities_goseq.R"))

# data ========================================================================
tpm <- read_csv(here(basedir, "data/03_compiled/tpm_no-testis.csv"))
counts <- read_csv(here(basedir, "data/03_compiled/counts_no-testis.csv"))
samples <- read_csv(here(basedir, "data/00_metadata/samples.csv"))

# for plotting
tissue_order <- c("stm.1.s", "int.1.w", "int.2.y", "int.3.w", "int.4.y")

genes <- read_tsv(here(genomedir, "results/genes/gene-display-names.tsv"))
orth <- read_tsv(here(genomedir, "results/ghost/gene-info_ghostdb_parsed.tsv"))

# prepare data ================================================================
df <- tibble(analysis = c("white", "yellow"),
             segments = list(c("int.1.w", "int.3.w"), # for setting levels
                             c("int.2.y", "int.4.y")))


df <- df %>% 
  mutate(md = map(analysis, ~ filter(samples, colour == .x))) %>% 
  mutate(md = map2(md, 
                   segments, 
                   ~ mutate(.x, tissue = factor(tissue, .y)))) %>% 
  mutate(cts = map(md, 
                   ~ counts %>% 
                     select(all_of(c("gene_id", .x$sample_name)))
  ))

# make sure that rows of colData correspond to the columns of countData
df <- df %>% 
  mutate(md = map2(md,
                   cts,
                   function(x, y) {
                     x[match(x$sample_name, 
                             names(y)[names(y) != "gene_id"]), ]
                   })) 

# run DESeq ===================================================================
p_cutoff <- 0.05

df$dds <- vector("list", nrow(df))
for (i in seq_len(nrow(df))) {
  df$dds[[i]] <- DESeqDataSetFromMatrix(countData = as.data.frame(df$cts[[i]]),
                                        colData = as.data.frame(df$md[[i]]),
                                        design = ~ tissue,
                                        tidy = TRUE)
  df$dds[[i]] <- DESeq(df$dds[[i]], test = "LRT", reduced = ~ 1)
}

# extract results =============================================================
df$res_shrunken_raw <- vector("list", nrow(df))
df$res_shrunken     <- vector("list", nrow(df))
df$res_raw          <- vector("list", nrow(df))
df$res              <- vector("list", nrow(df))
df$res_de           <- vector("list", nrow(df))

for (i in seq_len(nrow(df))) {
  
  # shunken l2fc for plotting
  df$res_shrunken_raw[[i]] <- lfcShrink(
    df$dds[[i]],
    coef = paste0("tissue_", df$segments[[i]][2], "_vs_", df$segments[[i]][1]),
    type = "apeglm"
  )
  
  df$res_shrunken[[i]] <- df$res_shrunken_raw[[i]] %>% 
    as.data.frame() %>% 
    filter(!is.na(padj)) %>% 
    rownames_to_column("gene_id") %>% 
    as_tibble()
  
  # original l2fc
  df$res_raw[[i]] <- results(df$dds[[i]], 
                             contrast = c("tissue", 
                                          df$segments[[i]][2], 
                                          df$segments[[i]][1]), 
                             alpha = p_cutoff)
  
  df$res[[i]] <- as.data.frame(df$res_raw[[i]]) %>% 
    filter(!is.na(padj)) %>% 
    rownames_to_column("gene_id") %>% 
    as_tibble()
  
  # DE genes only
  df$res_de[[i]] <- df$res[[i]] %>% 
    filter(padj < p_cutoff) %>% 
    arrange(desc(baseMean))
  
}

# separate up and down genes
df <- df %>% 
  mutate(res_up = map(res_de, ~ filter(.x, log2FoldChange > 0)),
         res_down = map(res_de, ~ filter(.x, log2FoldChange < 0)))

# write gene lists ============================================================

for (i in seq_len(nrow(df))) {
  write_csv(df$res_shrunken[[i]],
            here(resdir, df$analysis[[i]], "deseq_out_shrunken.csv"))
  write_csv(df$res[[i]],
            here(resdir, df$analysis[[i]], "deseq_out.csv"))
  
  fs::dir_create(here(resdir, df$analysis[[i]], "genes"))
  write_csv(df$res_de[[i]],   
            here(resdir, df$analysis[[i]], "genes", "de.csv"))
  write_csv(df$res_up[[i]],   
            here(resdir, df$analysis[[i]], "genes", paste0(df$segments[[i]][2], 
                                                           "_up.csv")))
  write_csv(df$res_down[[i]], 
            here(resdir, df$analysis[[i]], "genes", paste0(df$segments[[i]][1], 
                                                           "_up.csv")))
}

# plot p value distributions ==================================================
pval_dist_plot <- df %>% 
  select(analysis, res) %>% 
  unnest(col = res) %>% 
  ggplot(aes(pvalue)) +
  geom_histogram(binwidth = 0.05, boundary = 0, color = "white") +
  facet_grid(cols = vars(analysis)) +
  labs(caption = "white: int.3.w vs int.1.w\nyellow: int.4.y vs int.2.y") +
  theme_bw_2()

ggsave(here(resdir, "pvalue-distribution.pdf"),
       pval_dist_plot,
       width = 5, height = 2.75, device = cairo_pdf)

# plot number of DE genes =====================================================
n_de_plot <- df %>% 
  select(analysis, res_up, res_down) %>% 
  pivot_longer(cols = -analysis, 
               names_to = "change", names_prefix = "res_",
               values_to = "genes") %>% 
  mutate(n_genes = map_int(genes, nrow)) %>% 
  ggplot(aes(change, n_genes, label = n_genes)) +
  geom_col() +
  geom_text(vjust = 0, nudge_y = 20, size = 7/.pt) +
  facet_grid(cols = vars(analysis), scales = "free") +
  labs(caption = "white: int.3.w vs int.1.w\nyellow: int.4.y vs int.2.y") +
  theme_bw_2()

ggsave(here(resdir, "n_de_genes.pdf"), 
       n_de_plot, 
       width = 3, height = 3.5, device = cairo_pdf)

# plot individual DE genes ====================================================
# data
# de genelist in a long format 
de_genes <- df %>% 
  select(analysis, res_up, res_down) %>% 
  pivot_longer(cols = -analysis, 
               names_to = "change", names_prefix = "res_",
               values_to = "genes") %>% 
  mutate(change = case_when(
    analysis == "white"  & change == "up"   ~ "int.3.w_up",
    analysis == "white"  & change == "down" ~ "int.1.w_up",
    analysis == "yellow" & change == "up"   ~ "int.4.y_up",
    analysis == "yellow" & change == "down" ~ "int.2.y_up"
  ))

# tpm long
tpm_long <- tpm %>% 
  pivot_longer(cols = -gene_id, names_to = "sample_name", values_to = "tpm") %>% 
  left_join(samples) %>% 
  mutate(tissue = factor(tissue, tissue_order))

# functions
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

for (i in seq_len(nrow(de_genes))) {
  
  cairo_pdf(here(resdir, de_genes$analysis[[i]], 
                 "genes", paste0(de_genes$change[[i]], ".pdf")),
            width = 8.25, height = 11, onefile = TRUE)
  
  n_pages <- calculate_num_pages(length(de_genes$genes[[i]]$gene_id), 40)
  for (j in seq_len(n_pages)) {
    print(plot_genes_paginated(de_genes$genes[[i]]$gene_id, j))
  }
  
  dev.off()
  rm(n_pages)
  
}

# MA plots ===============================================================

plot_ma <- function(.res_shrunken, .level1, .level2) {
  
  MA_data <- left_join(.res_shrunken, genes, by = "gene_id")
  MA_labels <- MA_data %>% 
    filter(padj < 0.05, abs(log2FoldChange) > 1)
  
  ggplot(MA_data, aes(log10(baseMean), log2FoldChange)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.5) +
    geom_point(aes(color = padj < 0.05), alpha = 1, size = 1.2) +
    geom_text_repel(data = MA_labels, 
                    mapping = aes(label = symbol),
                    color = "blue", size = 8/.pt) +
    scale_color_manual(values = c("grey90", "grey40")) +
    labs(x = "Average Expression (log10)", 
         y = "log2 fold change",
         title = paste0(.level2, " vs ", .level1)) +
    theme_bw_2()
  
}

df$ma <- vector("list", nrow(df))
for (i in seq_len(nrow(df))) {
  df$ma[[i]] <- plot_ma(df$res_shrunken[[i]], 
                        df$segments[[i]][1], 
                        df$segments[[i]][2])
  
  ggsave(here(resdir, df$analysis[[i]], "MA-plot.pdf"),
         df$ma[[i]],
         width = 8, height = 6, device = cairo_pdf)
}

# GO enrichment ===============================================================

de_genes <- de_genes %>% 
  mutate(go = map(genes, 
                  ~ run_goseq(.x$gene_id, tpm$gene_id) %>% 
                    filter(p_adj < 0.05)))

for (i in seq_len(nrow(de_genes))) {
  
  fs::dir_create(here(resdir, de_genes$analysis[[i]], "go"))
  
  write_tsv(de_genes$go[[i]],
            here(resdir, de_genes$analysis[[i]], "go", 
                 paste0(de_genes$change[[i]], ".tsv")))
  
}

# end =========================================================================

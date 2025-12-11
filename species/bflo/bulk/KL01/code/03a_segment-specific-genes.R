# Arun Chavan
# Started: 2024-01-01

# background ==================================================================

# Here we perform DEseq2 analysis to find genes specific to each segment of the
# gut tract.

# setup =======================================================================
library(here)
library(DESeq2)
library(ggrepel)
library(tidyverse)
library(ggforce)

basedir <- here("species/bflo/bulk/KL01")
genomedir <- here("species/bflo/genome/VNyyK")
resdir <- here(basedir, "results/03a_segment-specific")
fs::dir_create(resdir)

source(here(basedir, "code/utilities_ggplot.R"))
source(here(basedir, "code/utilities_goseq.R"))

# data ========================================================================
tpm <- read_tsv(here(basedir, "data/03_compiled/tpm_gut.tsv"))
counts <- read_tsv(here(basedir, "data/03_compiled/counts_gut.tsv"))
samples <- read_tsv(here(basedir, "data/00_metadata/samples_gut.tsv"))

genes <- read_tsv(here(genomedir, "results/genes/gene-display-names.tsv"))
orth <- read_csv(here(genomedir, "results/emapper/processed/emapper_orthologs_mouse.csv"))

# for plotting
tissue_order <- c("mdiv", "mgut", "hgut")

# DESeq2 ======================================================================
df <- tibble(tissue = unique(samples$tissue))

# prepare metadata for each pairwise comparison
df <- df %>% 
  mutate(md = map(tissue, 
                  ~ samples %>% 
                    mutate(segment = ifelse(.data$tissue == .x,
                                            .x,
                                            "others")) %>% 
                    mutate(segment = factor(segment, c("others", .x)))))

# make sure that rows of colData correspond to the columns of countData
df <- df %>% 
  mutate(md = map(md, 
                  function(x) {
                    x[match(names(counts)[names(counts) != "gene_id"], 
                            x$sample_name), ]
                  })) 

# run DESeq2
p_cutoff <- 0.05

df$dds <- vector("list", nrow(df))
for (i in seq_len(nrow(df))) {
  df$dds[[i]] <- DESeqDataSetFromMatrix(countData = as.data.frame(counts),
                                        colData = as.data.frame(df$md[[i]]),
                                        design = ~ segment,
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
    coef = paste0("segment_", df$tissue[[i]], "_vs_others"),
    type = "apeglm"
  )
  
  df$res_shrunken[[i]] <- df$res_shrunken_raw[[i]] %>% 
    as.data.frame() %>% 
    filter(!is.na(padj)) %>% 
    rownames_to_column("gene_id") %>% 
    as_tibble()

  # original l2fc
  df$res_raw[[i]] <- results(df$dds[[i]], 
                             contrast = c("segment", df$tissue[[i]], "others"), 
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
walk(df$tissue, ~ fs::dir_create(here(resdir, .x)))

for (i in seq_len(nrow(df))) {
  write_csv(df$res_shrunken[[i]],
            here(resdir, df$tissue[[i]], "deseq_out_shrunken.csv"))
  write_csv(df$res[[i]],
            here(resdir, df$tissue[[i]], "deseq_out.csv"))
  
  fs::dir_create(here(resdir, df$tissue[[i]], "genes"))
  write_csv(df$res_de[[i]],   here(resdir, df$tissue[[i]], "genes", "de.csv"))
  write_csv(df$res_up[[i]],   here(resdir, df$tissue[[i]], "genes", "up.csv"))
  write_csv(df$res_down[[i]], here(resdir, df$tissue[[i]], "genes", "down.csv"))
}

# plot p value distributions ==================================================
pval_dist_plot <- df %>% 
  select(tissue, res) %>% 
  unnest(col = res) %>% 
  mutate(tissue = factor(tissue, tissue_order)) %>% 
  ggplot(aes(pvalue)) +
  geom_histogram(binwidth = 0.05, boundary = 0, color = "white") +
  facet_grid(cols = vars(tissue)) +
  theme_bw_2()
  
ggsave(here(resdir, "pvalue-distribution.pdf"),
       pval_dist_plot,
       width = 7, height = 3, device = cairo_pdf)

# plot number of DE genes =====================================================
n_de_plot <- df %>% 
  select(tissue, res_up, res_down) %>% 
  pivot_longer(cols = -tissue, 
               names_to = "change", names_prefix = "res_",
               values_to = "genes") %>% 
  mutate(n_genes = map_int(genes, nrow)) %>% 
  mutate(tissue = factor(tissue, tissue_order)) %>% 
  ggplot(aes(change, n_genes, label = n_genes)) +
  geom_col() +
  geom_text(vjust = 0, nudge_y = 20, size = 7/.pt) +
  facet_grid(cols = vars(tissue), scales = "free") +
  theme_bw_2()

ggsave(here(resdir, "n_de_genes.pdf"), 
       n_de_plot, 
       width = 4, height = 2.75, device = cairo_pdf)


# plot individual DE genes ====================================================
# data
# de genelist in a long format 
de_genes <- df %>% 
  select(tissue, res_up, res_down) %>% 
  pivot_longer(cols = -tissue, 
               names_to = "change", names_prefix = "res_",
               values_to = "genes")

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
                        nrow = 8, ncol = 6, page = .pg_num) +
    scale_x_discrete(guide = guide_axis(n.dodge = 1)) +
    theme_bw_2() +
    theme(strip.text = element_text(size = 6, hjust = 0))
}

for (i in seq_len(nrow(de_genes))) {
  
  cairo_pdf(here(resdir, de_genes$tissue[[i]], 
                 "genes", paste0(de_genes$change[[i]], ".pdf")),
            width = 8.25, height = 11, onefile = TRUE)
  
  n_pages <- calculate_num_pages(length(de_genes$genes[[i]]$gene_id), 48)
  for (j in seq_len(n_pages)) {
    print(plot_genes_paginated(de_genes$genes[[i]]$gene_id, j))
  }
  
  dev.off()
  rm(n_pages)
  
}

# MA plots ===============================================================

plot_ma <- function(.res_shrunken, .segment) {
  
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
         title = paste0(.segment, " vs others")) +
    theme_bw_2()
  
}

df <- df %>% 
  mutate(ma = map2(res_shrunken, tissue, plot_ma))

walk2(
  df$ma, 
  df$tissue,
  ~ ggsave(here(resdir, .y, "MA-plot.pdf"), 
           .x,
           width = 8, height = 6, device = cairo_pdf)
)

# GO enrichment ===============================================================

de_genes <- de_genes %>% 
  mutate(go = map(genes, 
                  ~ run_goseq(.x$gene_id, tpm$gene_id) %>% 
                    filter(p_adj < 0.05)))

for (i in seq_len(nrow(de_genes))) {
  
  fs::dir_create(here(resdir, de_genes$tissue[[i]], "go"))
  
  write_tsv(de_genes$go[[i]],
            here(resdir, de_genes$tissue[[i]], "go", 
                 paste0(de_genes$change[[i]], ".tsv")))
  
}

# generate metascape input ====================================================
get_orthologs <- function(.gene_ids, types = c("one2one")) {
  orth %>% 
    filter(cr_gene_id %in% .gene_ids) %>% 
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

de_genes <- de_genes %>% 
  mutate(one2one = map(genes, ~ get_orthologs(.x$gene_id)))

de_genes_up <- filter(de_genes, change == "up")

metascape <- generate_metascape_input(de_genes_up$tissue, de_genes_up$one2one)

fs::dir_create(here(resdir, "metascape"))
write_csv(metascape, here(resdir, "metascape/input.csv"))
# end =========================================================================
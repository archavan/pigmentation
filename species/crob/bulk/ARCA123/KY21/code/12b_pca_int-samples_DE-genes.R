# Arun Chavan
# Started: 2026-03-13

# background ==================================================================

# setup =======================================================================
library(DESeq2)
library(tidyverse)
library(here)
library(glue)
library(factoextra)
library(ggrepel)
library(colorspace)
library(pheatmap)
library(conflicted)

conflicts_prefer(
  dplyr::filter,
  purrr::map,
  dplyr::select,
  dplyr::count,
  base::as.factor
)

basedir <- here("species/crob/bulk/ARCA123/KY21")
genomedir <- here("species/crob/genome/KY21")

source(here(basedir, "code/utilities_ggplot.R"))
source(here(basedir, "code/utilities_jackstraw.R"))
source(here(basedir, "code/utilities_goseq.R"))

resdir <- here(basedir, "results/12b_pca_int-samples_DE-genes")
fs::dir_create(resdir)

# data ========================================================================
include_samples <- c("I1", "I2", "I3", "I4")
sample_order <- c("S", "I1", "I2", "I3", "I4") |>
  (function(x) x[x %in% include_samples])()

md <- read_csv(here(basedir, "data/00_metadata/samples.csv")) |>
  filter(tissue_short %in% include_samples) |>
  mutate(tissue_short = fct(tissue_short, sample_order))

tpm <- read_csv(here(basedir, "data/03_compiled/tpm_no-testis.csv")) |>
  rename_with(
    function(x) {
      ifelse(x %in% md$sample_name, md$tissue_name[match(x, md$sample_name)], x)
    }
  ) |>
  select(all_of(c("gene_id", md$tissue_name)))

cts <- read_csv(here(basedir, "data/03_compiled/counts_no-testis.csv")) |>
  rename_with(
    function(x) {
      ifelse(x %in% md$sample_name, md$tissue_name[match(x, md$sample_name)], x)
    }
  ) |>
  select(all_of(c("gene_id", md$tissue_name)))

genes <- read_tsv(here(genomedir, "results/genes/gene-display-names.tsv"))

mmgenes <- read_tsv(here(
  "~/Dropbox/research/projects/genomes/orthofinder/2025-08-29/results/id-name-mapping/mmus.tsv"
)) |>
  select(gene_id, symbol) |>
  distinct()

crmm <- read_tsv(here(
  "~/Dropbox/research/projects/genomes/orthofinder/2025-08-29/results/pairwise-orthologs/crob_vs_mmus.tsv"
)) |>
  left_join(mmgenes, by = c("mmus" = "gene_id"))


###############################################################################
# subset genes to include in PCA ==============================================
###############################################################################
# make sure that rows of colData correspond to the columns of countData
md <- md[match(md$tissue_name, names(cts)[names(cts) != "gene_id"]), ]
md <- md |> relocate(tissue_name)

# run DESeq
dds <- DESeqDataSetFromMatrix(
  countData = as.data.frame(cts),
  colData = as.data.frame(md),
  design = ~tissue_short,
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
res <- as.data.frame(res_raw) %>%
  filter(!is.na(padj)) %>%
  rownames_to_column("gene_id")

res_de <- res |>
  filter(padj < p_cutoff)

hist(res$pvalue)

include_genes <- res_de$gene_id

###############################################################################
# run PCA =====================================================================
###############################################################################

pca_dat <- tpm |>
  filter(gene_id %in% include_genes) |>
  column_to_rownames("gene_id") |>
  (function(x) x[which(rowSums(x) > 0), ])() |>
  sqrt()

js <- js_run_pca_with_jackstraw(pca_dat, 10, 1000, ndim = 9)

## tidy output ================================================================
ldngs <- js$pca_out$rotation |>
  as.data.frame() |>
  rownames_to_column("gene_id") |>
  pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "loading")

pvals <- js$pval |>
  as.data.frame() |>
  rownames_to_column("gene_id") |>
  pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "pval")

js_tidy <- left_join(ldngs, pvals, by = c("PC", "gene_id")) |>
  group_by(PC) |>
  mutate(padj = p.adjust(pval, method = "BH")) |>
  ungroup() |>
  mutate(PC = fct(PC, paste0("PC", sort(parse_number(unique(PC)))))) |>
  mutate(direction = ifelse(loading < 0, "negative", "positive")) |>
  arrange(PC, direction, -abs(loading)) |>
  left_join(select(genes, gene_id, display_name), by = "gene_id")


add_metadata_to_pc_coords <- function(.pca_out) {
  .pca_out$x |>
    as.data.frame() |>
    rownames_to_column("tissue_name") |>
    left_join(y = md, by = "tissue_name")
}

pc_space <- add_metadata_to_pc_coords(js$pca_out)

## get sig genes ==============================================================
sig_genes <- js_tidy |>
  filter(padj < 0.2) |> # for higher PCs, we need a liberal cutoff
  group_by(PC, direction) |>
  nest() |>
  ungroup()

## plot variance explained ====================================================

var_explained <- get_eig(js$pca_out) |>
  rownames_to_column("PC") |>
  mutate(PC = gsub("Dim\\.", "PC", PC)) |>
  mutate(variance.percent = round(variance.percent, 1))

var_explained |>
  mutate(PC = str_remove(PC, "PC")) |>
  mutate(PC = as.integer(PC)) |>
  ggplot(aes(PC, variance.percent)) +
  geom_point(size = 3) +
  geom_line() +
  scale_x_continuous(breaks = scales::breaks_pretty()) +
  scale_y_continuous(name = "Variance Explained (%)") +
  theme_bw_2()
ggsave(here(resdir, "pca_var-explained.png"), width = 5, height = 3.5)

## plot PC score against segment ==============================================
pc_space |>
  pivot_longer(
    cols = starts_with("PC"),
    names_to = "PC",
    values_to = "PC_score"
  ) |>
  ggplot(aes(sample, PC_score, color = as.factor(rep))) +
  geom_line(aes(group = rep), linewidth = 0.3) +
  geom_point(size = 1.5) +
  facet_wrap(vars(PC), scales = "free") +
  scale_x_discrete(name = "gut segment") +
  scale_y_continuous(name = "PC score") +
  scale_color_discrete(name = "animal") +
  theme_bw_2() +
  theme(
    strip.text = element_text(size = 11)
  )

ggsave(here(resdir, "pca_pc-score-vs-gut-segment.png"), width = 7.5, height = 6)

## plot PC score against animal ===============================================
pc_space |>
  pivot_longer(
    cols = starts_with("PC"),
    names_to = "PC",
    values_to = "PC_score"
  ) |>
  ggplot(aes(as.factor(rep_id), PC_score, color = tissue_short)) +
  geom_line(aes(group = tissue_short), linewidth = 0.3) +
  geom_point(size = 1.5, aes(fill = colour), shape = 21) +
  facet_wrap(vars(PC), scales = "free") +
  scale_x_discrete(name = "replicate id") +
  scale_y_continuous(name = "PC score") +
  scale_color_discrete(name = "gut segment") +
  scale_fill_identity() +
  theme_bw_2() +
  theme(
    strip.text = element_text(size = 11)
  )

ggsave(here(resdir, "pca_pc-score-vs-replicate.pdf"), width = 6.5, height = 6)

## plot ordinations ===========================================================

get_contrasting_text_color <- function(x) {
  lum <- colorspace::coords(as(colorspace::hex2RGB(x), "LAB"))[, "L"]
  ifelse(lum < 50, "white", "black")
}

plot_pc_ordination <- function(dim1, dim2) {
  pc_space |>
    arrange(rep_id, tissue_short) |>
    mutate(rep_id = as.factor(rep_id)) |>
    ggplot(aes(.data[[dim1]], .data[[dim2]], fill = tissue_short)) +
    geom_path(aes(group = animal, linetype = rep_id), color = "grey50") +
    geom_point(size = 9, shape = 21, stroke = 0) +
    geom_text(
      aes(
        label = tissue_name,
        color = after_scale(get_contrasting_text_color(fill))
      ),
      size = 3
    ) +
    scale_fill_viridis_d(name = "gut segment") +
    scale_linetype_manual(
      values = c("solid", "dotted", "dotdash"),
      name = "animal"
    ) +
    coord_equal() +
    labs(
      x = glue(
        "{dim1} ({varexp}%)",
        varexp = var_explained$variance.percent[var_explained$PC == dim1]
      ),
      y = glue(
        "{dim2} ({varexp}%)",
        varexp = var_explained$variance.percent[var_explained$PC == dim2]
      )
    ) +
    theme_bw(base_size = 16, base_family = .use_font) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text = element_text()
    )
}

ordinations <- expand_grid(
  dim1 = paste0("PC", 1:6),
  dim2 = paste0("PC", 1:6)
) |>
  filter(dim1 != dim2) |>
  mutate(xy = map2(dim1, dim2, ~ sort(c(.x, .y)))) |>
  filter(!duplicated(xy)) |>
  select(-xy) |>
  mutate(ord = map2(dim1, dim2, plot_pc_ordination))

fs::dir_create(here(resdir, "pca_ordinations"))

for (i in seq_len(nrow(ordinations))) {
  ggsave(
    here(
      resdir,
      "pca_ordinations",
      glue("{x}-{y}.png", x = ordinations$dim1[i], y = ordinations$dim2[i])
    ),
    plot = ordinations$ord[[i]],
    width = 9,
    height = 7
  )
}


###############################################################################
# GO enrichment ===============================================================
###############################################################################

sig_genes <- sig_genes |>
  mutate(
    go_bg_de = map(
      data,
      ~ run_goseq(.x$gene_id, .bg = res_de$gene_id) |>
        filter(p_adj < 0.05)
    )
  ) |>
  mutate(
    go_bg_all = map(
      data,
      ~ run_goseq(.x$gene_id, .bg = res$gene_id) |>
        filter(p_adj < 0.05)
    )
  )


fs::dir_create(here(resdir, "go/bg_de"))
fs::dir_create(here(resdir, "go/bg_all"))

for (i in seq_len(nrow(sig_genes))) {
  write_tsv(
    sig_genes$go_bg_de[[i]],
    here(
      resdir,
      "go/bg_de",
      paste0(sig_genes$PC[i], "_", sig_genes$direction[i], ".tsv")
    )
  )
  write_tsv(
    sig_genes$go_bg_all[[i]],
    here(
      resdir,
      "go/bg_all",
      paste0(sig_genes$PC[i], "_", sig_genes$direction[i], ".tsv")
    )
  )
}

sig_genes |>
  select(PC, direction, data) |>
  unnest(data) |>
  write_tsv(here(resdir, "jackstraw_sig-genes.tsv"))

###############################################################################
# plot sig genes ==============================================================
###############################################################################

fs::dir_create(here(resdir, "heatmaps"))

plot_heatmap_by_segment <- function(
  .genes,
  .filename,
  .width,
  .height,
  .save = FALSE
) {
  mean_tpm <- tpm |>
    filter(gene_id %in% .genes) |>
    pivot_longer(
      cols = -gene_id,
      names_to = "tissue_name",
      values_to = "tpm"
    ) |>
    left_join(md) |>
    group_by(tissue_short, gene_id) |>
    summarize(tpm = mean(tpm)) |>
    ungroup() |>
    pivot_wider(names_from = tissue_short, values_from = tpm) |>
    column_to_rownames("gene_id")

  rownames(mean_tpm) <- genes$display_name[match(
    rownames(mean_tpm),
    genes$gene_id
  )]

  if (.save) {
    pheatmap(
      mean_tpm,
      cluster_cols = FALSE,
      scale = "row",
      filename = .filename,
      width = .width,
      height = .height
    )
  }

  if (!.save) {
    pheatmap(
      mean_tpm,
      cluster_cols = FALSE,
      scale = "row"
    )
  }
}

get_pc_genes <- function(.pc, .direction) {
  sig_genes$data[[which(
    sig_genes$PC == .pc & sig_genes$direction == .direction
  )]]$gene_id
}


get_pc_genes("PC1", "negative") |>
  plot_heatmap_by_segment(
    .filename = here(resdir, "heatmaps", "PC1_negative.pdf"),
    .width = 11.35,
    .height = 300,
    .save = TRUE
  )

get_pc_genes("PC1", "positive") |>
  plot_heatmap_by_segment(
    .filename = here(resdir, "heatmaps", "PC1_positive.pdf"),
    .width = 11,
    .height = 165,
    .save = TRUE
  )


get_pc_genes("PC2", "negative") |>
  plot_heatmap_by_segment(
    .filename = here(resdir, "heatmaps", "PC2_negative.pdf"),
    .width = 19,
    .height = 55,
    .save = TRUE
  )

get_pc_genes("PC2", "positive") |>
  plot_heatmap_by_segment(
    .filename = here(resdir, "heatmaps", "PC2_positive.pdf"),
    .width = 22,
    .height = 95,
    .save = TRUE
  )


get_pc_genes("PC3", "negative") |>
  plot_heatmap_by_segment(
    .filename = here(resdir, "heatmaps", "PC3_negative.pdf"),
    .width = 9,
    .height = 7,
    .save = TRUE
  )

get_pc_genes("PC3", "positive") |>
  plot_heatmap_by_segment(
    .filename = here(resdir, "heatmaps", "PC3_positive.pdf"),
    .width = 8,
    .height = 1.5,
    .save = TRUE
  )

###############################################################################
# get mouse homologs of PC genes for metascape ================================
###############################################################################

sig_genes <- sig_genes |>
  mutate(
    mm_orth = map(data, function(x) {
      crmm |>
        filter(crob %in% x$gene_id) |>
        pull(symbol) |>
        unique()
    })
  ) |>
  mutate(geneset_name = paste0(PC, "_", direction))

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

meta <- generate_metascape_input(sig_genes$geneset_name, sig_genes$mm_orth)

fs::dir_create(here(resdir, "metascape"))
write_csv(meta, here(resdir, "metascape", "input.csv"))

# end =========================================================================

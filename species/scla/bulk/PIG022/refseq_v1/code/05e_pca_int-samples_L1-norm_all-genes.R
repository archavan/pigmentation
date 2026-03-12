# Arun Chavan
# Started: 2025-09-03

# background ==================================================================

# setup =======================================================================
library(tidyverse)
library(here)
library(glue)
library(factoextra)
library(ggrepel)
library(colorspace)
library(pheatmap)
library(conflicted)

basedir <- here("species/scla/bulk/PIG022/refseq_v1")
genomedir <- here("species/scla/genome/refseq_v1")

source(here(basedir, "code/utilities_ggplot.R"))
source(here(basedir, "code/utilities_jackstraw.R"))
source(here(basedir, "code/utilities_goseq.R"))

resdir <- here(basedir, "results/pca/int-samples_L1-norm_all-genes")
fs::dir_create(resdir)

conflicts_prefer(
  dplyr::filter,
  purrr::map,
  dplyr::select,
  dplyr::count,
  base::as.factor
)

# data ========================================================================
include_samples <- c("I1", "I2", "I3", "I4", "I5", "I6")
sample_order <- c("S", "I1", "I2", "I3", "I4", "I5", "I6") |>
  (function(x) x[x %in% include_samples])()

md <- read_csv(here(basedir, "data/samples_clean.csv")) |>
  filter(sample %in% include_samples) |>
  mutate(sample = fct(sample, sample_order))

tpm <- read_csv(here(basedir, "data/03_compiled/tpm_clean.csv")) |>
  select(all_of(c("gene_id", md$sample_id)))

cts <- read_csv(here(basedir, "data/03_compiled/counts_clean.csv")) |>
  select(all_of(c("gene_id", md$sample_id)))

uss <- read_csv(here(basedir, "data/03_compiled/uss_intestine.csv")) |>
  select(all_of(c("gene_id", md$sample_id)))

genes <- read_tsv(here(genomedir, "results/genes/gene-display-names.tsv"))

mmgenes <- read_tsv(here(
  "~/Dropbox/research/projects/genomes/orthofinder/2025-08-29/results/id-name-mapping/mmus.tsv"
)) |>
  select(gene_id, symbol) |>
  distinct()

scmm <- read_tsv(here(
  "~/Dropbox/research/projects/genomes/orthofinder/2025-08-29/results/pairwise-orthologs/scla_vs_mmus.tsv"
)) |>
  left_join(mmgenes, by = c("mmus" = "gene_id"))

###############################################################################
# run PCA =====================================================================
###############################################################################

pca_dat <- uss |>
  column_to_rownames("gene_id")

js <- js_run_pca_with_jackstraw(pca_dat, 10, 1000, ndim = 9, scale. = FALSE)

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
    rownames_to_column("sample_id") |>
    left_join(y = md, by = "sample_id")
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
  ggplot(aes(as.factor(rep), PC_score, color = sample)) +
  geom_line(aes(group = sample), linewidth = 0.3) +
  geom_point(size = 1.5) +
  facet_wrap(vars(PC), scales = "free") +
  scale_x_discrete(name = "replicate id") +
  scale_y_continuous(name = "PC score") +
  scale_color_discrete(name = "gut segment") +
  theme_bw_2() +
  theme(
    strip.text = element_text(size = 11)
  )

ggsave(here(resdir, "pca_pc-score-vs-replicate.png"), width = 6.5, height = 6)

## plot ordinations ===========================================================

get_contrasting_text_color <- function(x) {
  lum <- colorspace::coords(as(colorspace::hex2RGB(x), "LAB"))[, "L"]
  ifelse(lum < 50, "white", "black")
}

plot_pc_ordination <- function(dim1, dim2) {
  pc_space |>
    arrange(rep, sample) |>
    mutate(rep = as.factor(rep)) |>
    ggplot(aes(.data[[dim1]], .data[[dim2]], fill = sample)) +
    geom_path(aes(group = animal, linetype = rep), color = "grey50") +
    geom_point(size = 9, shape = 21, stroke = 0) +
    geom_text(
      aes(
        label = sample_id,
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
  dim1 = paste0("PC", 1:4),
  dim2 = paste0("PC", 1:4)
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
    go_bg_all = map(
      data,
      ~ run_goseq(.x$gene_id, .bg = rownames(pca_dat)) |>
        filter(p_adj < 0.05)
    )
  )


fs::dir_create(here(resdir, "go/bg_all"))

for (i in seq_len(nrow(sig_genes))) {
  write_tsv(
    sig_genes$go_bg_all[[i]],
    here(
      resdir,
      "go/bg_all",
      paste0(sig_genes$PC[i], "_", sig_genes$direction[i], ".tsv")
    )
  )
}

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
  mean_uss <- uss |>
    filter(gene_id %in% .genes) |>
    pivot_longer(cols = -gene_id, names_to = "sample_id", values_to = "uss") |>
    left_join(md) |>
    group_by(sample, gene_id) |>
    summarize(uss = mean(uss)) |>
    ungroup() |>
    pivot_wider(names_from = sample, values_from = uss) |>
    column_to_rownames("gene_id")

  rownames(mean_uss) <- genes$display_name[match(
    rownames(mean_uss),
    genes$gene_id
  )]

  if (.save) {
    pheatmap(
      mean_uss,
      cluster_cols = FALSE,
      scale = "none",
      filename = .filename,
      width = .width,
      height = .height
    )
  }

  if (!.save) {
    pheatmap(
      mean_uss,
      cluster_cols = FALSE,
      scale = "none"
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
    .width = 10,
    .height = 200,
    .save = TRUE
  )

get_pc_genes("PC1", "positive") |>
  plot_heatmap_by_segment(
    .filename = here(resdir, "heatmaps", "PC1_positive.pdf"),
    .width = 10.5,
    .height = 150,
    .save = TRUE
  )


get_pc_genes("PC2", "negative") |>
  plot_heatmap_by_segment(
    .filename = here(resdir, "heatmaps", "PC2_negative.pdf"),
    .width = 9.8,
    .height = 45,
    .save = TRUE
  )

get_pc_genes("PC2", "positive") |>
  plot_heatmap_by_segment(
    .filename = here(resdir, "heatmaps", "PC2_positive.pdf"),
    .width = 10.25,
    .height = 100,
    .save = TRUE
  )


###############################################################################
# get mouse homologs of PC genes for metascape ================================
###############################################################################
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

meta_genes <- sig_genes |>
  filter(PC %in% paste0("PC", 1:2)) |>
  mutate(
    mm_orth = map(data, function(x) {
      scmm |>
        filter(scla %in% x$gene_id) |>
        pull(symbol) |>
        unique()
    })
  ) |>
  mutate(geneset_name = paste0(PC, "_", direction))


meta <- generate_metascape_input(meta_genes$geneset_name, meta_genes$mm_orth)

fs::dir_create(here(resdir, "metascape"))
write_csv(meta, here(resdir, "metascape", "input.csv"))

# end =========================================================================

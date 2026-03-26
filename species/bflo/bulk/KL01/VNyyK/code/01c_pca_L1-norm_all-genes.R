# Arun Chavan
# Started: 2026-03-26

# setup =======================================================================
library(tidyverse)
library(here)
library(glue)
library(factoextra)
library(ggrepel)
library(colorspace)
library(pheatmap)
library(conflicted)

basedir <- here("species/bflo/bulk/KL01/VNyyK")
genomedir <- here("species/bflo/genome/VNyyK")

source(here(basedir, "code/utilities_ggplot.R"))
source(here(basedir, "code/utilities_jackstraw.R"))
source(here(basedir, "code/utilities_goseq.R"))

resdir <- here(basedir, "results/01c_pca_L1-norm_all-genes")
fs::dir_create(resdir)

conflicts_prefer(
  dplyr::filter,
  purrr::map,
  dplyr::select,
  dplyr::count,
  base::as.factor
)

# data ========================================================================
sample_order <- c("mdiv", "mgut", "hgut")

md <- read_tsv(here(basedir, "data/00_metadata/samples_gut.tsv")) |>
  mutate(tissue = fct(tissue, sample_order))

tpm <- read_tsv(here(basedir, "data/03_compiled/tpm_gut.tsv"))
cts <- read_tsv(here(basedir, "data/03_compiled/counts_gut.tsv"))
uss <- read_tsv(here(basedir, "data/03_compiled/uss.tsv"))

genes <- read_tsv(here(genomedir, "results/genes/gene-display-names.tsv"))

mmgenes <- read_tsv(here(
  "~/Dropbox/research/projects/genomes/orthofinder/2025-08-29/results/id-name-mapping/mmus.tsv"
)) |>
  select(gene_id, symbol) |>
  distinct()

bfmm <- read_tsv(here(
  "~/Dropbox/research/projects/genomes/orthofinder/2025-08-29/results/pairwise-orthologs/bflo_vs_mmus.tsv"
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
    rownames_to_column("sample_name") |>
    left_join(y = md, by = "sample_name")
}

pc_space <- add_metadata_to_pc_coords(js$pca_out)

## get sig genes ==============================================================
sig_genes <- js_tidy |>
  filter(padj < 0.2) |> # for higher PCs, we need a liberal cutoff
  group_by(PC, direction) |>
  nest() |>
  ungroup()

fs::dir_create(here(resdir, "jackstraw-genes"))
for (i in seq_len(nrow(sig_genes))) {
  write_tsv(
    sig_genes$data[[i]],
    here(
      resdir,
      "jackstraw-genes",
      paste0(sig_genes$PC[[i]], "_", sig_genes$direction[[i]], ".tsv")
    )
  )
}

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
  ggplot(aes(tissue, PC_score, color = as.factor(rep_id))) +
  geom_line(aes(group = rep_id), linewidth = 0.3) +
  geom_point(size = 1.5) +
  facet_wrap(vars(PC), scales = "free") +
  scale_x_discrete(name = "gut segment") +
  scale_y_continuous(name = "PC score") +
  scale_color_discrete(name = "animal") +
  theme_bw_2() +
  theme(
    strip.text = element_text(size = 11)
  )

ggsave(here(resdir, "pca_pc-score-vs-gut-segment.pdf"), width = 6, height = 5)

## plot PC score against animal ===============================================
pc_space |>
  pivot_longer(
    cols = starts_with("PC"),
    names_to = "PC",
    values_to = "PC_score"
  ) |>
  ggplot(aes(as.factor(rep_id), PC_score, color = tissue)) +
  geom_line(aes(group = tissue), linewidth = 0.3) +
  geom_point(size = 1.5) +
  facet_wrap(vars(PC), scales = "free") +
  scale_x_discrete(name = "replicate id") +
  scale_y_continuous(name = "PC score") +
  scale_color_discrete(name = "gut segment") +
  theme_bw_2() +
  theme(
    strip.text = element_text(size = 11)
  )

ggsave(here(resdir, "pca_pc-score-vs-replicate.pdf"), width = 6, height = 5)

## plot ordinations ===========================================================
p_ord <- pc_space |>
  ggplot(aes(PC1, PC2, fill = tissue)) +
  geom_polygon(
    data = function(x) {
      x |>
        slice(chull(PC1, PC2), .by = tissue)
    },
    alpha = 0.5
  ) +
  geom_point(shape = 21, size = 4, show.legend = FALSE) +
  geom_text(aes(label = rep_id), size = 2.5) +
  geom_text(
    data = function(x) {
      x |>
        summarize(PC1 = median(PC1), PC2 = median(PC2), .by = tissue)
    },
    mapping = aes(PC1 + 1, PC2 + 1, label = tissue),
    size = 3,
    inherit.aes = FALSE
  ) +
  scale_fill_discrete(name = "gut segment") +
  labs(
    x = glue(
      "PC1 ({varexp}%)",
      varexp = var_explained$variance.percent[var_explained$PC == "PC1"]
    ),
    y = glue(
      "PC2 ({varexp}%)",
      varexp = var_explained$variance.percent[var_explained$PC == "PC2"]
    )
  ) +
  theme_bw_2()

ggsave(here(resdir, "ordination_PC1-PC2.pdf"), p_ord, width = 6, height = 4)


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
      bfmm |>
        filter(bflo %in% x$gene_id) |>
        pull(symbol) |>
        unique()
    })
  ) |>
  mutate(geneset_name = paste0(PC, "_", direction))


meta <- generate_metascape_input(meta_genes$geneset_name, meta_genes$mm_orth)

fs::dir_create(here(resdir, "metascape"))
write_csv(meta, here(resdir, "metascape", "input.csv"))

# end =========================================================================

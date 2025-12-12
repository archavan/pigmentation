# Arun Chavan
# Started: 2025-12-01

# setup =======================================================================
library(tidyverse)
library(here)
library(glue)
library(ggrepel)

basedir <- here("species/scla/bulk/PIG022/refseq_v1")
genomedir <- here("species/scla/genome/refseq_v1")
resdir <- here("species/scla/bulk/PIG022/refseq_v1/results/09a_zwick-genes")
fs::dir_create(resdir)

mmdir <- here("species/mmus/zwick")
orthdir <- here("~/Dropbox/research/projects/genomes/orthofinder/2025-08-29")

source(here(basedir, "code/utilities_ggplot.R"))

###############################################################################
# data ========================================================================
###############################################################################

tpm <- read_csv(here(basedir, "data/03_compiled/tpm.csv"))
samples <- read_csv(here(basedir, "data/samples.csv"))

domain_markers <- read_csv(here(mmdir, "results/markers/top-markers.csv"))

###############################################################################
# genomic data ================================================================
###############################################################################

genes <- read_tsv(here(genomedir, "results/genes/gene-display-names.tsv"))

## orthology data ====
mgenes <- read_tsv(here(orthdir, "results/id-name-mapping/mmus.tsv")) %>%
  select(gene_id, symbol) %>%
  distinct()

orth <- read_tsv(here(orthdir, "results/pairwise-orthologs/mmus_vs_scla.tsv"))

orth <- left_join(orth, mgenes, by = c("mmus" = "gene_id"))
orth <- left_join(
  orth,
  select(genes, gene_id, display_name),
  by = c("scla" = "gene_id")
)

get_orthologs <- function(x) {
  stopifnot(length(x) == 1)
  stopifnot(x %in% orth$symbol)
  orth |>
    filter(symbol == x) |>
    pull(scla) |>
    unique()
}

domain_markers <- domain_markers |>
  mutate(ortholog = map(gene, possibly(get_orthologs))) |>
  unnest(col = ortholog)

###############################################################################
# prepare for plotting ========================================================
###############################################################################

tissue_order <- c("S", "I1", "I2", "I3", "I4", "I5", "I6")

samples <- samples |>
  mutate(tissue_short = fct(sample, tissue_order))

tpm_long <- tpm |>
  pivot_longer(cols = -gene_id, names_to = "sample_id", values_to = "tpm") |>
  left_join(samples, by = "sample_id")

display_names <- deframe(select(genes, gene_id, display_name))

###############################################################################
# faceted gene plots ==========================================================
###############################################################################

fs::dir_create(here(resdir, "faceted-genes"))

order_genes_by_expression <- function(.genes) {
  tpm_long |>
    filter(gene_id %in% .genes) |>
    group_by(gene_id) |>
    summarize(tpm = mean(tpm)) |>
    ungroup() |>
    arrange(-tpm) |>
    pull(gene_id)
}

extract_facet_layout <- function(p) {
  q <- ggplot_build(p)
  n_facet_rows <- max(summarise_layout(q)$row)
  n_facet_cols <- max(summarise_layout(q)$col)
  list(nrows = n_facet_rows, ncols = n_facet_cols)
}

plot_genes_faceted <- function(.genes, .domain = NULL, .save = FALSE) {
  .title <- paste0("Domain ", .domain)
  gene_order <- order_genes_by_expression(.genes)
  pdat <- tpm_long |>
    filter(gene_id %in% .genes) |>
    mutate(gene_id = fct(gene_id, gene_order))

  pp <- ggplot(pdat, aes(tissue_short)) +
    geom_jitter(aes(y = tpm), width = 0.1, color = "grey") +
    stat_summary(aes(y = tpm, group = 1), fun = mean, geom = "line") +
    stat_summary(aes(y = tpm), fun = mean, geom = "point") +
    facet_wrap(
      vars(gene_id),
      scales = "free_y",
      labeller = labeller(gene_id = display_names)
    ) +
    labs(title = .title, x = "segment") +
    theme_bw(base_size = 8) +
    theme(
      strip.text = element_text(size = 5, hjust = 0),
      plot.title = element_text(size = 8, face = 1),
      panel.grid.minor = element_blank()
    )

  if (.save) {
    nfacets <- extract_facet_layout(pp)
    ggsave(
      filename = here(
        resdir,
        "faceted-genes",
        paste0("domain-", .domain, ".png")
      ),
      plot = pp,
      width = nfacets$ncols * 2.2,
      height = nfacets$nrows * 1.6 + 0.6,
      dpi = 300
    )
  } else {
    return(pp)
  }
}


for (i in unique(domain_markers$domain)) {
  domain_markers |>
    filter(domain == i) |>
    pull(ortholog) |>
    unique() |>
    plot_genes_faceted(i, .save = TRUE)
}


###############################################################################
#  heatmaps ===================================================================
###############################################################################
fs::dir_create(here(resdir, "heatmaps"))

find_expressed_genes <- function(.genes) {
  tpm_long |>
    filter(gene_id %in% .genes) |>
    group_by(gene_id) |>
    summarize(max_tpm = max(tpm)) |>
    ungroup() |>
    filter(max_tpm > 2) |>
    pull(gene_id)
}

order_genes_by_expression_pattern <- function(.pdat) {
  .pdat |>
    group_by(gene_id) |>
    slice_max(ztpm, n = 1, with_ties = FALSE) |>
    ungroup() |>
    arrange(tissue_short) |>
    group_by(tissue_short) |>
    arrange(desc(ztpm), .by_group = TRUE) |>
    ungroup() |>
    pull(gene_id) |>
    as.character() |>
    rev()
}

plot_module_heatmap <- function(.genes, .domain = NULL, .save = FALSE) {
  .title <- paste0("Domain ", .domain)
  pdat <- tpm_long |>
    filter(gene_id %in% find_expressed_genes(.genes)) |>
    group_by(gene_id) |>
    mutate(ztpm = scale(tpm)[, 1]) |>
    ungroup()
  pdat <- pdat |>
    mutate(gene_id = fct(gene_id, order_genes_by_expression_pattern(pdat)))

  pp <- ggplot(pdat, aes(x = tissue_short, y = gene_id, fill = ztpm)) +
    geom_tile(color = "black") +
    scale_y_discrete(
      labels = function(x) str_trunc(display_names[x], 90, "right"),
      position = "right"
    ) +
    scale_x_discrete(name = "segment") +
    scale_fill_gradient(low = "white", high = "blue") +
    coord_fixed() +
    labs(title = .title) +
    theme_bw(base_size = 8) +
    theme(
      legend.position = "left",
      axis.title.y = element_blank(),
      axis.text.y = element_text(hjust = 0),
      plot.title = element_text(size = 8, face = 1),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank()
    )

  if (.save) {
    ngenes <- length(unique(pdat$gene_id))
    ggsave(
      filename = here(
        resdir,
        "heatmaps",
        paste0("domain-", .domain, ".png")
      ),
      plot = pp,
      width = 6,
      height = (ngenes * 0.15) + 1,
      dpi = 300
    )
  } else {
    return(pp)
  }
}

for (i in unique(domain_markers$domain)) {
  domain_markers |>
    filter(domain == i) |>
    pull(ortholog) |>
    unique() |>
    plot_module_heatmap(i, .save = TRUE)
}

# end =========================================================================

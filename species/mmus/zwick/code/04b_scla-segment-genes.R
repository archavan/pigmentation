# Arun Chavan
# Started: 2025-12-12

# setup =======================================================================
library(tidyverse)
library(here)
library(glue)
library(Seurat)
library(ggrepel)

basedir <- here("species/mmus/zwick")
genomedir <- here("species/mmus/genome")
orthdir <- here("~/Dropbox/research/projects/genomes/orthofinder/2025-08-29")
scdir <- here("species/scla/bulk/PIG022/refseq_v1")
scgenomedir <- here("species/scla/genome/refseq_v1")

resdir <- here(basedir, "results/scla-segment-genes")
fs::dir_create(resdir)

source(here(basedir, "code/utilities_ggplot.R"))
source(here(basedir, "code/utilities_seurat-plots.R"))

###############################################################################
# data ========================================================================
###############################################################################

ec <- read_rds(here(basedir, "data/enterocytes.rds"))
ec$domain <- factor(ec$domain, LETTERS[1:5])
ec <- SetIdent(ec, value = "domain")

ps <- read_csv(here(basedir, "data/enterocytes_agg-by-domain.csv"))
ps_long <- ps |>
  pivot_longer(
    cols = -gene,
    names_to = "domain",
    values_to = "norm_expr"
  ) |>
  mutate(
    domain = factor(domain, levels = LETTERS[1:5])
  )

sc_kmeans_genes <- tibble(clust = paste0("cluster-", 1:12))
sc_kmeans_genes <- sc_kmeans_genes |>
  mutate(
    genes = map(
      clust,
      ~ read_csv(here(
        scdir,
        "results/cluster-genes/all-samples/k12/genes",
        paste0(.x, ".csv")
      ))
    )
  )


###############################################################################
# genomic data ================================================================
###############################################################################

scgenes <- read_tsv(here(scgenomedir, "results/genes/gene-display-names.tsv"))

## orthology data ====
mgenes <- read_tsv(here(orthdir, "results/id-name-mapping/mmus.tsv")) %>%
  select(gene_id, symbol) %>%
  distinct()

orth <- read_tsv(here(orthdir, "results/pairwise-orthologs/mmus_vs_scla.tsv"))

orth <- left_join(orth, mgenes, by = c("mmus" = "gene_id"))
orth <- left_join(
  orth,
  select(scgenes, gene_id, display_name),
  by = c("scla" = "gene_id")
)

get_orthologs <- function(x) {
  stopifnot(length(x) == 1)
  stopifnot(x %in% orth$scla)
  orth |>
    filter(scla == x) |>
    pull(symbol) |>
    unique()
}


mm_orthologs <- sc_kmeans_genes |>
  unnest(col = genes) |>
  mutate(mm_orth = map(gene_id, possibly(get_orthologs))) |>
  select(clust, mm_orth) |>
  unnest(col = mm_orth) |>
  group_by(clust) |>
  summarize(mm_orth = list(unique(mm_orth)))


###############################################################################
# plotting setup ==============================================================
###############################################################################

plot_dotplot2 <- function(
  .features,
  .title,
  .reorder.genes = TRUE
) {
  dd <- dp_extract_dotplot_data(ec, .features)
  if (.reorder.genes) {
    dd <- dp_define_feature_order(dd)
  }
  dp_split_dotplot_data(dd) |>
    purrr::map(dp_plot) |>
    purrr::map(dp_add_title, .title) |>
    purrr::map(function(p) {
      p +
        theme(plot.title = element_text(size = 8, face = 1))
    })
}

###############################################################################
# dotplots ====================================================================
###############################################################################

normalize_expression <- function(ct, total_ct) {
  log1p((ct / total_ct) * 10000)
}

ec_agg <- AggregateExpression(
  ec,
  group.by = "orig.ident",
  return.seurat = FALSE
)$RNA |>
  as.data.frame() |>
  rownames_to_column("gene") |>
  rename(cts = V1) |>
  mutate(norm_exp = normalize_expression(cts, sum(cts)))

get_top_expressed_genes <- function(.genes, .n = 50) {
  ec_agg |>
    filter(gene %in% .genes) |>
    slice_max(order_by = norm_exp, n = .n) |>
    pull(gene)
}

mm_orthologs <- mm_orthologs |>
  mutate(mm_orth_topexp = map(mm_orth, get_top_expressed_genes)) |>
  mutate(dp = map2(mm_orth_topexp, clust, ~ plot_dotplot2(.x, .y, TRUE)))

fs::dir_create(here(resdir, "dotplots"))

walk2(
  mm_orthologs$dp,
  mm_orthologs$clust,
  ~ save_dotplot(.x, here(resdir, "dotplots", paste0(.y, ".pdf")))
)

# end =========================================================================

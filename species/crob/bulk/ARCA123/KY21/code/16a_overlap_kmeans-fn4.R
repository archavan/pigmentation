# Arun Chavan
# Started: 2026-03-12

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)
library(glue)
library(ggtext)

basedir <- here("species/crob/bulk/ARCA123/KY21")
genomedir <- here("species/crob/genome/KY21")

kmeansdir <- here(basedir, "results/07a_cluster-genes_int-samples/k10")

resdir <- here(basedir, "results/16a_overlap_kmeans-fn4")
fs::dir_create(here(resdir))

source(here(basedir, "code/utilities_fisher-overlap-test.R"))

###############################################################################
# genomic data ================================================================
###############################################################################

genes <- read_tsv(here(genomedir, "results/genes/gene-display-names.tsv"))
fnmod <- read_tsv(here(genomedir, "results/genes/functional-modules_gut.tsv"))

fn <- fnmod |>
  select(module.1, module.2, module.3, module.4, ortholog) |>
  filter(!is.na(ortholog)) |>
  mutate(
    mod = pmap_chr(
      list(module.1, module.2, module.3, module.4),
      ~ c(...) |> na.omit() |> paste0(collapse = " - ")
    )
  ) |>
  summarize(gene_ids = list(ortholog) |> unique(), .by = mod)

###############################################################################
# data ========================================================================
###############################################################################

## expression data ====

tpm <- read_csv(here(basedir, "data/03_compiled/tpm_no-testis.csv"))
uss <- read_csv(here(basedir, "data/03_compiled/uss_intestine.csv"))
samples <- read_csv(here(basedir, "data/00_metadata/samples.csv")) |>
  filter(tissue_short != "S")

tpm <- tpm |>
  select(-starts_with("stm")) |>
  rename_with(
    function(x) samples$tissue_name[match(x, samples$sample_name)],
    .cols = -gene_id
  )

## kmeans clustering results ====
km <- list.files(here(kmeansdir, "genes")) |>
  set_names() |>
  map(~ read_csv(here(kmeansdir, "genes", .x))) |>
  bind_rows(.id = "cluster") |>
  mutate(cluster = str_remove(cluster, ".csv")) |>
  mutate(cluster = str_replace(cluster, "cluster-", "kcr")) |>
  summarize(gene_ids = list(gene_id) |> unique(), .by = cluster)

###############################################################################
# run overlap test ============================================================
###############################################################################

p_thresh <- 0.05

## define universe ============================================================

universe <- tpm |>
  column_to_rownames("gene_id") |>
  (function(x) x[rowSums(x) > 0, ])() |>
  rownames()

fn <- fn |>
  mutate(gene_ids = map(gene_ids, ~ .x[.x %in% universe]))

## run test ===================================================================

olaps <- expand_grid(
  mod = fn$mod,
  kclust = km$cluster
) |>
  mutate(
    genes_mod = fn$gene_ids[match(mod, fn$mod)],
    genes_kclust = km$gene_ids[match(kclust, km$cluster)]
  ) |>
  mutate(overlap = map2(genes_mod, genes_kclust, intersect)) |>
  mutate(
    n_genes_mod = map_int(genes_mod, length),
    n_genes_kclust = map_int(genes_kclust, length),
    n_overlap = map_int(overlap, length)
  ) |>
  filter(n_overlap > 0) |>
  mutate(
    fisher = map2(
      genes_mod,
      genes_kclust,
      \(x, y) {
        run_fisher_overlap_test(x, y, universe = universe)
      },
      .progress = TRUE
    )
  ) |>
  mutate(
    odds_ratio = map_dbl(fisher, \(f) f$estimate),
    pval = map_dbl(fisher, \(f) f$p.value)
  ) |>
  mutate(adj_pval = p.adjust(pval, method = "BH")) |>
  arrange(adj_pval)

## diagnostics ================================================================
ggplot(olaps, aes(pval)) +
  geom_histogram(bins = 20, boundary = 0, color = "white") +
  theme_bw()

ggsave(here(resdir, "overlap-pval-histogram.pdf"), width = 4, height = 3)

## save results ===============================================================
olaps |>
  select(-genes_mod, -genes_kclust, -overlap, -fisher) |>
  write_tsv(here(resdir, "overlap-stats.tsv"))

olaps |>
  filter(adj_pval < p_thresh) |>
  select(mod, kclust, overlap) |>
  unnest(overlap) |>
  rename(gene_id = overlap) |>
  left_join(select(genes, gene_id, display_name)) |>
  write_tsv(here(resdir, "overlap-genes_significant.tsv"))

###############################################################################
# plot ========================================================================
###############################################################################

## ordering & filtering for plot ----
hclust_order <- function(mat) {
  h <- hclust(dist(mat))
  h$labels[h$order]
}

mat_padj <- olaps |>
  pivot_wider(
    id_cols = kclust,
    names_from = mod,
    values_from = adj_pval
  ) |>
  column_to_rownames("kclust") |>
  as.matrix()
mat_padj[is.na(mat_padj)] <- 1

order_kclust <- hclust_order(mat_padj)
order_mod <- hclust_order(t(mat_padj))

sig_mod <- olaps |>
  filter(adj_pval < p_thresh, .by = mod) |>
  pull(mod) |>
  unique()

sig_kclust <- olaps |>
  filter(adj_pval < p_thresh, .by = kclust) |>
  pull(kclust) |>
  unique()


olaps |>
  filter(mod %in% sig_mod, kclust %in% sig_kclust) |>
  mutate(
    # Inf odds ratio don't get plotted. Replace with max finite OR.
    odds_ratio_capped = pmin(
      odds_ratio,
      max(odds_ratio[is.finite(odds_ratio)])
    ),
    kclust = fct(kclust, levels = order_kclust),
    mod = fct(mod, levels = order_mod)
  ) |>
  ggplot(aes(kclust, mod)) +
  geom_point(
    aes(
      fill = -log10(adj_pval),
      size = log(odds_ratio_capped),
      color = adj_pval < p_thresh
    ),
    shape = 21,
    stroke = 0.5
  ) +
  scale_size(range = c(0, 6), limits = c(0, NA), name = "log(Odds Ratio)") +
  scale_fill_viridis_c(
    direction = -1,
    name = "--log<sub>10</sub>(P<sub>adj</sub>)"
  ) +
  scale_color_manual(
    values = c("TRUE" = "#fa0202ff", "FALSE" = "white"),
    name = glue("P<sub>adj</sub> < {p_thresh}")
  ) +
  scale_x_discrete(
    position = "top",
    name = "kmeans cluster",
  ) +
  scale_y_discrete(
    position = "right",
    name = "functional module"
  ) +
  coord_fixed() +
  theme_bw() +
  theme(
    axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0),
    legend.position = "left",
    legend.direction = "vertical",
    legend.title.position = "left",
    legend.title = element_markdown(hjust = 0.5, angle = 90),
    plot.caption = element_markdown(hjust = 1)
  )

quartz(
  type = "pdf",
  file = here(resdir, "overlaps-significant.pdf"),
  width = 6.75,
  height = 6
)
print(last_plot())
dev.off()

# end #########################################################################

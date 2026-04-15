# Arun Chavan
# Started: 2026-03-18

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)
library(glue)
library(pvclust)
library(ggtree)
library(ggtext)
library(conflicted)

conflicts_prefer(
  dplyr::filter,
  dplyr::select,
  dplyr::rename,
  dplyr::slice
)

basedir <- here("figures")
dir_cr <- here("species/crob/bulk/ARCA123/KY21")
dir_cr_genome <- here("species/crob/genome/KY21")

resdir <- here(basedir, "results/01a_crob_fn-mods")
fs::dir_create(resdir)

source(here(basedir, "code/utils_themes.R"), local = TRUE)
source(here(basedir, "code/utils_colors.R"), local = TRUE)

###############################################################################
# data ========================================================================
###############################################################################

## genomic data ===============================================================
fn4 <- read_tsv(here(dir_cr_genome, "results/genes/functional-modules_gut.tsv"))

fn4_idmap <- fn4 |>
  select(id, label) |>
  distinct() |>
  arrange(id)

## bulk RNA-seq data ==========================================================

md_cr <- read_csv(here(dir_cr, "data/00_metadata/samples.csv")) |>
  select(tissue_name, tissue_short, sample_name, rep_id, colour) |>
  rename(sample_id = tissue_name, sample = tissue_short) |>
  filter(sample != "S") |>
  mutate(
    pigment = case_when(
      sample %in% c("I2", "I4") ~ "dark",
      sample %in% c("I1", "I3") ~ "white"
    )
  )

stopifnot(!anyNA(md_cr$pigment))

tpm_cr <- read_csv(here(dir_cr, "data/03_compiled/tpm_no-testis.csv")) |>
  select(gene_id, all_of(md_cr$sample_name)) |>
  rename_with(function(x) {
    ifelse(
      x %in% md_cr$sample_name,
      md_cr$sample_id[match(x, md_cr$sample_name)],
      x
    )
  })

uss_cr <- read_csv(here(dir_cr, "data/03_compiled/uss_intestine.csv"))

###############################################################################
# PCA =========================================================================
###############################################################################

pca_int <- uss_cr |>
  column_to_rownames("gene_id") |>
  t() |>
  prcomp()

pca_int_var <- summary(pca_int)$importance["Proportion of Variance", ] * 100
pca_int_x <- pca_int$x |>
  as.data.frame() |>
  rownames_to_column("sample_id") |>
  inner_join(md_cr, by = "sample_id")

p_pca <- ggplot(pca_int_x, aes(x = PC1, y = PC2)) +
  geom_polygon(
    data = function(x) {
      x |>
        group_by(sample) |>
        slice(chull(PC1, PC2)) |>
        ungroup()
    },
    aes(x = PC1, y = PC2, group = sample, fill = pigment),
    color = "black",
    linewidth = lwd_thin,
  ) +
  geom_point(shape = 21, stroke = lwd_medium, aes(fill = pigment)) +
  geom_text(
    data = tribble(
      ~PC1  , ~PC2 , ~sample ,
      -8    ,    5 , "I3"    ,
      -4    ,   -2 , "I1"    ,
       4.75 ,   -4 , "I2"    ,
       7    ,    4 , "I4"
    ),
    aes(x = PC1, y = PC2, label = sample),
    size = ts_base / .pt
  ) +
  scale_fill_manual(values = clrs_pigment, guide = "none") +
  labs(
    x = glue("PC1 ({round(pca_int_var[['PC1']])}%)"),
    y = glue("PC2 ({round(pca_int_var[['PC2']])}%)")
  ) +
  theme_figure() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(linewidth = lwd_thin)
  )

save_quartz(here(resdir, "pca.pdf"), .width = 1.6, .height = 1.3, p_pca)

###############################################################################
# hierarchical clustering =====================================================
###############################################################################
set.seed(2398853)
pv <- uss_cr |>
  column_to_rownames("gene_id") |>
  pvclust(nboot = 1000, parallel = FALSE)

p_pv <- ggtree(pv$hclust, linewidth = 0.35) %<+%
  md_cr +
  geom_tiplab(offset = 0.05, size = ts_base / .pt) +
  geom_tippoint(
    aes(fill = pigment),
    shape = 21,
    size = 1.5,
    stroke = lwd_medium
  ) +
  scale_fill_manual(values = clrs_pigment) +
  scale_x_continuous(limits = c(NA, 0.125)) +
  theme_figure() +
  theme(
    panel.border = element_blank(),
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )

save_quartz(here(resdir, "hclust.pdf"), .width = 1.35, .height = 1.45, p_pv)

###############################################################################
# module heatmap ==============================================================
###############################################################################

sample_labels <- md_cr |>
  distinct(sample, pigment) |>
  mutate(
    label = glue(
      "<span style='color:{clrs_pigment_small[pigment]}; font-size:{ts_large}pt'>\u25cf</span> {sample}"
    )
  ) |>
  pull(label, name = sample)

# calculate median USS for each module and then scale
uss_long <- uss_cr |>
  pivot_longer(-gene_id, names_to = "sample_id", values_to = "uss") |>
  left_join(md_cr, by = "sample_id")

zuss_fn4 <- fn4 |>
  arrange(id) |>
  filter(!is.na(ortholog)) |>
  rename(gene_id = ortholog) |>
  select(id, label, gene_id) |>
  inner_join(uss_long, by = "gene_id") |>
  summarize(uss = median(uss), .by = c(id, label, sample, gene_id)) |>
  summarize(uss = median(uss), .by = c(id, label, sample)) |>
  mutate(zuss = scale(uss)[, 1], .by = c(id, label))

# for manual inspection
zuss_fn4 |>
  pivot_wider(names_from = sample, values_from = zuss, id_cols = label) |>
  column_to_rownames("label") |>
  pheatmap::pheatmap(cluster_cols = FALSE)

# categorize modules by their scaled median expression
enriched_pigmented <- c(
  "Luminal hydrolysis",
  "BB hydrolysis: carb",
  "BB hydrolysis: peptide",
  "Uptake: carb",
  "Uptake: AA",
  "Uptake: lipid",
  "Endocytosis",
  "Storage: lipid droplet"
)

enriched_unpigmented <- c(
  "Endocrine",
  "Synthesis: AA",
  "Synthesis: cholesterol",
  "Catabolism: carb",
  "Catabolism: FA",
  "Catabolism: AA",
  "Storage: glycogen",
  "Export: AA"
)

other <- fn4_idmap$label[
  !fn4_idmap$label %in% c(enriched_pigmented, enriched_unpigmented)
]

fn4_idmap <- fn4_idmap |>
  mutate(
    enriched_in = case_when(
      label %in% enriched_pigmented ~ "pigmented",
      label %in% enriched_unpigmented ~ "unpigmented",
      .default = "other"
    )
  ) |>
  arrange(enriched_in, id)

write_tsv(fn4_idmap, here(resdir, "module-expression-pattern.tsv"))

# plot

p_fn4_enriched <- zuss_fn4 |>
  inner_join(fn4_idmap, by = "label") |>
  mutate(
    label = fct(
      label,
      c(rev(enriched_pigmented), rev(enriched_unpigmented), other)
    )
  ) |>
  filter(enriched_in != "other") |>
  ggplot(aes(sample, label, fill = zuss)) +
  geom_tile(color = "white", linewidth = lwd_thin) +
  facet_wrap(
    vars(enriched_in),
    ncol = 1,
    space = "free_y",
    scales = "free_y",
    labeller = labeller(enriched_in = function(x) {
      glue("Enriched in<br>**{x}**")
    })
  ) +
  scale_fill_gradient2(
    low = clrs_gradient2$low,
    mid = clrs_gradient2$mid,
    high = clrs_gradient2$high,
    name = "z(median(normalized expression))",
    breaks = c(-1, 0, 1)
  ) +
  scale_x_discrete(labels = sample_labels) +
  scale_y_discrete(position = "right") +
  theme_figure() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(linewidth = lwd_thin, color = "black"),
    axis.title = element_blank(),
    axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom",
    legend.location = "plot",
    legend.title.position = "top",
    legend.title = element_markdown(
      hjust = 0.5,
      size = ts_small,
      angle = 0,
      margin = margin(b = 2)
    ),
    legend.text = element_text(size = ts_small, margin = margin(t = 2)),
    legend.key.width = unit(0.75, "lines"),
    legend.key.height = unit(0.25, "line"),
    legend.box.spacing = unit(0.1, "lines"),
    strip.clip = "off",
    strip.text = element_markdown(
      hjust = 0,
      vjust = 0,
      margin = margin(b = 1),
    )
  )

save_quartz(
  here(resdir, "module-heatmap.pdf"),
  .width = 1.42,
  .height = 2.9,
  p_fn4_enriched
)

###############################################################################
# plot individual genes from select enriched modules ==========================
###############################################################################

display_names <- fn4 |>
  filter(!is.na(ortholog)) |>
  select(gene_id = ortholog, hog_label) |>
  filter(!is.na(hog_label)) |>
  separate_wider_delim(hog_label, "_", names = c("nparalogs", "hog_symbol")) |>
  mutate(
    short_display_name = glue(
      "{hog_symbol}<span style='color:\"{muted_text}\"; font-size:{ts_tiny}pt'> / {gene_id} / {nparalogs}</span>"
    )
  ) |>
  select(gene_id, short_display_name) |>
  deframe()

# keep most informative and key genes for brevity
genes_to_exclude <- c(
  # BB hydrolysis: peptide
  "KY21.Chr9.1143",
  "KY21.Chr9.534",
  "KY21.Chr11.1215",
  # Uptake: carb (only keep active transporters)
  "KY21.Chr7.811",
  "KY21.Chr7.760",
  # Uptake: AA (only keep Slc15)
  "KY21.Chr3.979",
  "KY21.Chr3.660",
  "KY21.Chr3.135",
  "KY21.Chr2.764",
  "KY21.Chr2.461",
  "KY21.Chr2.1323",
  # Uptake: lipid (only keep Cd36)
  "KY21.Chr4.54",
  "KY21.Chr2.1376",
  # Catabolism: FA
  "KY21.Chr8.983"
)

genes_to_plot <- fn4 |>
  select(id, label, gene_id = ortholog) |>
  filter(!is.na(gene_id)) |>
  filter(gene_id %in% names(display_names)) |>
  filter(!gene_id %in% genes_to_exclude) |>
  filter(
    label %in%
      c(
        "BB hydrolysis: carb",
        "BB hydrolysis: peptide",
        "Uptake: carb",
        "Uptake: AA",
        "Uptake: lipid",
        "Endocytosis",
        "Synthesis: AA",
        "Synthesis: cholesterol",
        "Catabolism: carb",
        "Catabolism: FA",
        "Catabolism: AA",
        "Storage: glycogen"
      )
  ) |>
  distinct()

# NOTE: Gls is a part of both Synthesis: AA and Catabolism: AA.

genes_to_plot_data <- uss_long |>
  inner_join(genes_to_plot, by = "gene_id") |>
  summarize(median_uss = median(uss), .by = c(gene_id, sample, id, label)) |>
  mutate(zuss = scale(median_uss)[, 1], .by = c(gene_id, id, label)) |>
  left_join(fn4_idmap, by = c("id", "label")) |>
  mutate(label = fct(label, c(enriched_pigmented, enriched_unpigmented)))

p_genes_select <- genes_to_plot_data |>
  ggplot(aes(sample, gene_id, fill = zuss)) +
  geom_tile(color = "white", linewidth = lwd_thin) +
  facet_wrap(
    vars(label),
    space = "free_y",
    scales = "free_y",
  ) +
  scale_fill_gradient2(
    low = clrs_gradient2$low,
    mid = clrs_gradient2$mid,
    high = clrs_gradient2$high,
    name = "z(normalized expression)",
    breaks = c(-1, 0, 1),
    limits = c(-1.5, 1.5)
  ) +
  scale_y_discrete(position = "right", labels = display_names) +
  scale_x_discrete(labels = sample_labels) +
  theme_figure() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(linewidth = lwd_thin, color = "black"),
    axis.title = element_blank(),
    axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y.right = element_markdown(size = ts_small),
    legend.position = "bottom",
    legend.location = "plot",
    legend.title.position = "top",
    legend.title = element_markdown(
      hjust = 0.5,
      size = ts_small,
      angle = 0,
      margin = margin(b = 2)
    ),
    legend.text = element_text(size = ts_small, margin = margin(t = 2)),
    legend.key.width = unit(0.75, "lines"),
    legend.key.height = unit(0.25, "line"),
    legend.box.spacing = unit(0.1, "lines"),
    strip.clip = "off",
    strip.text = element_markdown(
      hjust = 0,
      vjust = 0,
      margin = margin(b = 1)
    )
  )

save_quartz(
  here(resdir, "gene-heatmap.pdf"),
  .width = 2.15,
  .height = 5.5,
  p_genes_select
)

# end =========================================================================

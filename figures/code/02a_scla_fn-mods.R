# Arun Chavan
# Started: 2026-04-15

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
dir_sc <- here("species/scla/bulk/PIG022/refseq_v1")
dir_sc_genome <- here("species/scla/genome/refseq_v1")

resdir <- here(basedir, "results/02a_scla_fn-mods")
fs::dir_create(resdir)

source(here(basedir, "code/utils_themes.R"), local = TRUE)
source(here(basedir, "code/utils_colors.R"), local = TRUE)

###############################################################################
# data ========================================================================
###############################################################################

## genomic data ===============================================================
fn4 <- read_tsv(here(dir_sc_genome, "results/genes/functional-modules_gut.tsv"))

fn4_idmap <- fn4 |>
  select(id, label) |>
  distinct() |>
  arrange(id)

## bulk RNA-seq data ==========================================================

md_sc <- read_csv(here(dir_sc, "data/samples_clean.csv")) |>
  select(sample_id, sample, rep) |>
  filter(sample != "S") |>
  mutate(
    pigment = case_when(
      sample %in% c("I1") ~ "white",
      sample %in% c("I2", "I3", "I4") ~ "dark",
      sample %in% c("I5", "I6") ~ "light"
    )
  )

stopifnot(!anyNA(md_sc$pigment))

tpm_sc <- read_csv(here(dir_sc, "data/03_compiled/tpm_clean.csv")) |>
  select(gene_id, all_of(md_sc$sample_id))

uss_sc <- read_csv(here(dir_sc, "data/03_compiled/uss_intestine.csv"))

## crob module enrichment =====================================================
# We want to only plot the modules that are enriched in crob pigmented and
# unpigmented segments, and specifically just the ones that highlight the
# division of labor between uptake and metabolism.
crob_enrichment <- read_tsv(here(
  basedir,
  "results/01a_crob_fn-mods/module-expression-pattern.tsv"
))

###############################################################################
# PCA =========================================================================
###############################################################################

pca_int <- uss_sc |>
  column_to_rownames("gene_id") |>
  t() |>
  prcomp()

pca_int_var <- summary(pca_int)$importance["Proportion of Variance", ] * 100
pca_int_x <- pca_int$x |>
  as.data.frame() |>
  rownames_to_column("sample_id") |>
  inner_join(md_sc, by = "sample_id")


label_nudge <- tribble(
  ~sample , ~nudge_x , ~nudge_y ,
  "I1"    , -0.6     ,  1.4     ,
  "I2"    ,  0.0     ,  1.0     ,
  "I3"    , -1.0     , -0.4     ,
  "I4"    ,  0.3     ,  0.7     ,
  "I5"    , -1.2     ,  0.2     ,
  "I6"    ,  0.6     ,  0.8     ,
)

label_pos <- pca_int_x |>
  summarize(PC1 = median(PC1), PC2 = median(PC2), .by = sample) |>
  left_join(label_nudge, by = "sample")

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
  geom_point(shape = 21, stroke = 0.35, aes(fill = pigment)) +
  geom_text(
    data = label_pos,
    aes(PC1 + nudge_x, PC2 + nudge_y, label = sample),
    size = ts_small / .pt
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
# module heatmap ==============================================================
###############################################################################

sample_labels <- md_sc |>
  distinct(sample, pigment) |>
  mutate(
    label = glue(
      "<span style='color:{clrs_pigment_small[pigment]}; font-size:{ts_large}pt'>\u25cf</span> {sample}"
    )
  ) |>
  pull(label, name = sample)

# calculate median USS for each module and then scale
uss_long <- uss_sc |>
  pivot_longer(-gene_id, names_to = "sample_id", values_to = "uss") |>
  inner_join(md_sc, by = "sample_id")

zuss_fn4 <- fn4 |>
  arrange(id) |>
  filter(!is.na(ortholog)) |>
  rename(gene_id = ortholog) |>
  select(id, label, gene_id) |>
  inner_join(uss_long, by = "gene_id") |>
  summarize(uss = median(uss), .by = c(id, label, sample, gene_id)) |> # over replicates
  summarize(uss = median(uss), .by = c(id, label, sample)) |> # over modules
  mutate(zuss = scale(uss)[, 1], .by = c(id, label))

keep_mods <- c(
  "Uptake: carb",
  "BB hydrolysis: carb",
  "Endocytosis",
  "BB hydrolysis: peptide",
  "Uptake: AA",
  "Uptake: lipid",
  "Catabolism: FA",
  "Synthesis: cholesterol",
  "Catabolism: AA",
  "Synthesis: AA",
  "Catabolism: carb",
  "Storage: glycogen"
)

fn4_categories <- crob_enrichment |>
  filter(label %in% keep_mods) |>
  mutate(
    category = case_when(
      enriched_in == "pigmented" ~ "Uptake",
      enriched_in == "unpigmented" ~ "Metabolism",
      .default = NA_character_
    )
  ) |>
  mutate(category = fct(category, c("Uptake", "Metabolism")))

# plot

p_fn4_enriched <- zuss_fn4 |>
  inner_join(fn4_categories, by = c("id", "label")) |>
  mutate(label = fct(label, keep_mods)) |>
  ggplot(aes(sample, label, fill = zuss)) +
  geom_tile(color = "white", linewidth = lwd_thin) +
  facet_wrap(
    vars(category),
    ncol = 1,
    space = "free_y",
    scales = "free_y",
    labeller = labeller(category = function(x) {
      glue("**{x}**")
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
  .width = 1.6,
  .height = 2.2,
  p_fn4_enriched
)

# end =========================================================================

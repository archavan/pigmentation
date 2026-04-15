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

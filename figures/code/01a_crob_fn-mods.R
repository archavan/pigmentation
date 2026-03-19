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

basedir <- here("figures")
dir_cr <- here("species/crob/bulk/ARCA123/KY21")
dir_cr_10x <- here(
  "~/Dropbox/research/projects/endoderm/species/crob/10x/run-02/KY21"
)
dir_cr_genome <- here("species/crob/genome/KY21")

resdir <- here(basedir, "results/01a_crob_fn-mods")
fs::dir_create(resdir)

source(here(basedir, "code/utils_themes.R"))
source(here(basedir, "code/utils_colors.R"))

###############################################################################
# data ========================================================================
###############################################################################

## genomic data ===============================================================
genes <- read_tsv(here(dir_cr_genome, "results/genes/gene-display-names.tsv"))

## bulk RNA-seq data ==========================================================

md_cr <- read_csv(here(dir_cr, "data/00_metadata/samples.csv")) |>
  select(tissue_name, tissue_short, sample_name, rep_id, colour) |>
  rename(sample_id = tissue_name, sample = tissue_short) |>
  filter(sample != "S") |>
  mutate(
    pigment = case_when(
      sample %in% c("I2", "I4") ~ "light",
      sample %in% c("I1", "I3") ~ "white"
    )
  )

tpm_cr <- read_csv(here(dir_cr, "data/03_compiled/tpm_no-testis.csv")) |>
  select(gene_id, all_of(md_cr$sample_name)) |>
  rename_with(\(x) {
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
  column_to_rownames("gene_id") %>%
  t() |>
  prcomp()

pca_int_var <- summary(pca_int)[[6]][2, ] * 100
pca_int_x <- pca_int$x |>
  as.data.frame() |>
  rownames_to_column("sample_id") |>
  (function(.y) right_join(x = md_cr, y = .y, by = "sample_id"))()

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
      ~PC1  , ~PC2  , ~sample ,
      -4    ,  5.2  , "I3"    ,
      -5    , -4    , "I1"    ,
       4.75 , -2.75 , "I2"    ,
       7.5  ,  3.5  , "I4"
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

quartz(type = "pdf", file = here(resdir, "pca.pdf"), width = 1.5, height = 1.2)
print(p_pca)
dev.off()

###############################################################################
# hierarchical clustering =====================================================
###############################################################################

pv <- uss_cr |>
  column_to_rownames("gene_id") |>
  pvclust(parallel = TRUE)

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

quartz(
  type = "pdf",
  file = here(resdir, "hclust.pdf"),
  width = 1.3,
  height = 1.35
)
print(last_plot())
dev.off()

# end #########################################################################

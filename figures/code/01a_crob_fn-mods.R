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
genes <- read_tsv(here(dir_cr_genome, "results/genes/gene-display-names.tsv"))
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

quartz(type = "pdf", file = here(resdir, "pca.pdf"), width = 1.5, height = 1.2)
print(p_pca)
dev.off()

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

quartz(
  type = "pdf",
  file = here(resdir, "hclust.pdf"),
  width = 1.3,
  height = 1.35
)
print(last_plot())
dev.off()

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
  select(id, label, gene_id, display_name) |>
  inner_join(uss_long, by = "gene_id") |>
  summarize(uss = median(uss), .by = c(id, label, sample)) |>
  mutate(zuss = scale(uss)[, 1], .by = c(id, label))

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
  "Export: AA",
  "Bile acid metabolism"
)

other <- fn4_idmap$label[
  !fn4_idmap$label %in% c(enriched_pigmented, enriched_unpigmented)
]

fn4_idmap <- fn4_idmap |>
  mutate(
    group = case_when(
      label %in% enriched_pigmented ~ "Enriched in<br>**pigmented**",
      label %in% enriched_unpigmented ~ "Enriched in<br>**unpigmented**",
      .default = "Other"
    )
  ) |>
  arrange(group, id)

write_tsv(fn4_idmap, here(resdir, "module-expression-pattern.tsv"))

# plot
zuss_fn4 |>
  inner_join(fn4_idmap, by = "label") |>
  mutate(
    label = fct(
      label,
      c(rev(enriched_pigmented), rev(enriched_unpigmented), other)
    )
  ) |>
  filter(group != "Other") |>
  ggplot(aes(sample, label, fill = zuss)) +
  geom_tile(color = "white", linewidth = lwd_thin) +
  facet_wrap(
    vars(group),
    ncol = 1,
    space = "free_y",
    scales = "free_y"
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
      margin = margin(b = 1)
    ),
    legend.text = element_text(size = ts_small),
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

quartz(
  type = "pdf",
  file = here(resdir, "module-heatmap.pdf"),
  width = 1.42,
  height = 3
)
print(last_plot())
dev.off()

# end #########################################################################

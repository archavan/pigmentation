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
fn4 <- read_tsv(here(dir_cr_genome, "results/genes/functional-modules_gut.tsv"))

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


uss_long <- uss_cr |>
  pivot_longer(-gene_id, names_to = "sample_id", values_to = "uss") |>
  left_join(md_cr, by = "sample_id")

zuss_fn4 <- fn4 |>
  filter(!is.na(ortholog)) |>
  rename(gene_id = ortholog) |>
  select(id, label, gene_id, display_name) |>
  inner_join(uss_long, by = "gene_id") |>
  summarize(uss = median(uss), .by = c(id, label, sample)) |>
  mutate(zuss = scale(uss)[, 1], .by = c(id, label))

keep_mods <- zuss_fn4 |>
  select(id, label) |>
  distinct() |>
  arrange(id) |>
  filter(
    !label %in%
      c(
        "Tissue barrier",
        "Mucus barrier",
        "Ciliary motility",
        "Export: AA",
        "Acidic digestion"
      )
  )


zuss_fn4_wide <- zuss_fn4 |>
  filter(label %in% keep_mods$label) |>
  pivot_wider(id_cols = "label", names_from = "sample", values_from = "zuss") |>
  column_to_rownames("label")

fn4_order <- zuss_fn4_wide |>
  dist() |>
  hclust() |>
  (function(x) x$label[x$order])()

fn4_order_manual <- tribble(
  ~label                   , ~group                    ,
  "BB hydrolysis: peptide" , "Enriched in pigmented"   ,
  "BB hydrolysis: carb"    , "Enriched in pigmented"   ,
  "Uptake: AA"             , "Enriched in pigmented"   ,
  "Uptake: carb"           , "Enriched in pigmented"   ,
  "Uptake: lipid"          , "Enriched in pigmented"   ,
  "Endocytosis"            , "Enriched in pigmented"   ,
  "Storage: lipid droplet" , "Enriched in pigmented"   ,
  "Storage: glycogen"      , "Enriched in unpigmented" ,
  "Bile acid metabolism"   , "Enriched in unpigmented" ,
  "Catabolism: AA"         , "Enriched in unpigmented" ,
  "Catabolism: FA"         , "Enriched in unpigmented" ,
  "Catabolism: carb"       , "Enriched in unpigmented" ,
  "Synthesis: AA"          , "Enriched in unpigmented" ,
  "Synthesis: cholesterol" , "Enriched in unpigmented" ,
  "Endocrine"              , "Enriched in unpigmented" ,
  "Urea cycle"             , "Other"                   ,
  "Phagocytosis"           , "Other"                   ,
  "Reesterification"       , "Other"                   ,
  "Detoxification"         , "Other"                   ,
  "Synthesis: carb"        , "Other"                   ,
  "Luminal hydrolysis"     , "Other"                   ,
  "Antimicrobial"          , "Other"                   ,
  "Lysosomal digestion"    , "Other"                   ,
  "Uptake: micronutrient"  , "Other"                   ,
  "Export: carb"           , "Other"                   ,
  "Efflux: cholesterol"    , "Other"                   ,
  "Transport: lipid"       , "Other"                   ,
  "Chylomicron assembly"   , "Other"
)

stopifnot(length(fn4_order[!fn4_order %in% fn4_order_manual$label]) == 0)

zuss_fn4 |>
  filter(label %in% fn4_order_manual$label) |>
  left_join(fn4_order_manual, by = "label") |>
  mutate(label = fct(label, fn4_order_manual$label)) |>
  ggplot(aes(sample, label, fill = zuss)) +
  geom_tile(color = "white") +
  facet_wrap(
    ncol = 1,
    vars(group),
    space = "free_y",
    scales = "free_y",
    labeller = label_wrap_gen(50),
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
    strip.text = element_text(
      hjust = 0,
      vjust = 0,
      margin = margin(b = 3),
      face = 2
    )
  )

quartz(
  type = "pdf",
  file = here(resdir, "module-heatmap.pdf"),
  width = 1.4,
  height = 4
)
print(last_plot())
dev.off()

# end #########################################################################

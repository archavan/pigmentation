# Arun Chavan
# Started: 2026-05-10

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(readxl)
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
resdir <- here(basedir, "results/07a_lipidomics")
fs::dir_create(here(resdir))

datadir <- here("expts/PIG005_lipidomics/Expt5_Parastou")
xlpath <- here(datadir, "20250716Lipid6500_Exp5_CionaLipid.xlsx")

source(here(basedir, "code/utils_themes.R"), local = TRUE)
source(here(basedir, "code/utils_colors.R"), local = TRUE)

clrs_group <- c(
  S = clrs_pigment_small$light,
  W = clrs_pigment_small$white,
  Y = clrs_pigment_small$dark
)

###############################################################################
# data ========================================================================
###############################################################################

add_labels <- function(df) {
  df |>
    mutate(
      group = case_when(
        GroupName == "Stomach" ~ "S",
        GroupName == "IntW" ~ "W",
        GroupName == "IntY" ~ "Y"
      ),
      rep = row_number(),
      .by = GroupName
    ) |>
    mutate(sample = paste0(group, "_", rep)) |>
    select(
      -c(`...1`, ExpNum, GroupNum, GroupName, SampleID, SampleNorm, NormType)
    )
}

# Use Composition sheets (not weight-normalized; see email from Kevin -- weight
# estimation is unreliable due to water content): class % of total lipid,
# species % within class.

## raw composition data ========================================================

lpd_class <- read_excel(xlpath, sheet = "Class Composit") |> add_labels()

row_totals <- lpd_class |>
  select(-c(group, rep, sample)) |>
  rowSums(na.rm = TRUE)
stopifnot(all(near(row_totals, 100)))

lpd_species <- read_excel(xlpath, sheet = "Species Composit") |> add_labels()

known_classes <- lpd_class |> select(-c(group, rep, sample)) |> names()
spc_names <- lpd_species |> select(-c(group, rep, sample)) |> names()

## species-to-class mapping ====================================================

# Species Composit uses "FA" prefix; the corresponding class in Class Composit is "FFA".
class_prefixes <- setNames(known_classes, known_classes)
class_prefixes["FFA"] <- "FA"

pfx_pattern <- paste0(
  "^(",
  paste(
    class_prefixes[order(nchar(class_prefixes), decreasing = TRUE)],
    collapse = "|"
  ),
  ")(?=[ /]|$)"
)
species_to_class <- tibble(
  species = spc_names,
  class = names(class_prefixes)[match(
    str_extract(spc_names, pfx_pattern, group = 1),
    class_prefixes
  )]
)
stopifnot(!anyNA(species_to_class$class))

spc_mat <- lpd_species |> select(-c(group, rep, sample))
class_sums <- split.default(spc_mat, species_to_class$class) |>
  vapply(function(cls) rowSums(cls, na.rm = TRUE), numeric(nrow(spc_mat)))
stopifnot(all(near(class_sums, 100)))

## species % of total lipid ====================================================

# species % of total = (species % within class) × (class % of total) / 100

lpd_species_pct <- lpd_species |>
  pivot_longer(
    -c(group, rep, sample),
    names_to = "species",
    values_to = "pct_within"
  ) |>
  left_join(species_to_class, by = "species") |>
  left_join(
    lpd_class |>
      pivot_longer(
        -c(group, rep, sample),
        names_to = "class",
        values_to = "pct_class"
      ),
    by = c("group", "rep", "sample", "class")
  ) |>
  mutate(pct_total = pct_within * pct_class / 100) |>
  select(group, rep, sample, species, pct_total) |>
  pivot_wider(names_from = species, values_from = pct_total)

row_totals_spc <- lpd_species_pct |>
  select(-c(group, rep, sample)) |>
  rowSums(na.rm = TRUE)
stopifnot(all(near(row_totals_spc, 100)))

md_lpd <- lpd_class |> select(sample, group)

###############################################################################
# PCA =========================================================================
###############################################################################

lpd_mat <- lpd_species_pct |>
  select(-c(group, rep)) |>
  column_to_rownames("sample") |>
  mutate(across(everything(), function(x) replace_na(x, 0))) |>
  as.matrix()

lpd_mat <- lpd_mat[, apply(lpd_mat, 2, var) > 0]

pca_res <- prcomp(lpd_mat, scale. = TRUE)
pca_var <- summary(pca_res)$importance["Proportion of Variance", ] * 100

pca_x <- pca_res$x |>
  as.data.frame() |>
  rownames_to_column("sample") |>
  left_join(md_lpd, by = "sample")

label_nudge <- tribble(
  ~group , ~nudge_x , ~nudge_y ,
  "S"    ,       -7 ,       -3 ,
  "W"    ,        4 ,       -2 ,
  "Y"    ,       -9 ,        3
)

label_pos <- pca_x |>
  summarize(PC1 = median(PC1), PC2 = median(PC2), .by = group) |>
  left_join(label_nudge, by = "group")

p_pca <- ggplot(pca_x, aes(x = PC1, y = PC2)) +
  geom_polygon(
    data = function(x) x |> slice(chull(PC1, PC2), .by = group),
    aes(group = group, fill = group),
    color = "black",
    linewidth = lwd_thin,
    alpha = 0.3
  ) +
  geom_point(shape = 21, stroke = lwd_medium, aes(fill = group)) +
  geom_text(
    data = label_pos,
    mapping = aes(x = PC1 + nudge_x, y = PC2 + nudge_y, label = group),
    size = ts_small,
    size.unit = "pt"
  ) +
  scale_fill_manual(values = clrs_group, guide = "none") +
  labs(
    x = glue("PC1 ({round(pca_var[['PC1']])}%)"),
    y = glue("PC2 ({round(pca_var[['PC2']])}%)")
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

hclust_xlim <- 0.00015 # right margin for tip labels

set.seed(4719)
pv <- pvclust(t(lpd_mat), nboot = 1000, parallel = FALSE)

p_hclust <- ggtree(pv$hclust, linewidth = lwd_thin) %<+%
  md_lpd +
  geom_tiplab(offset = 0.001, size = ts_base / .pt) +
  geom_tippoint(
    aes(fill = group),
    shape = 21,
    size = 1.5,
    stroke = lwd_medium
  ) +
  scale_fill_manual(values = clrs_group) +
  scale_x_continuous(limits = c(NA, hclust_xlim)) +
  theme_figure() +
  theme(
    panel.border = element_blank(),
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )

save_quartz(here(resdir, "hclust.pdf"), .width = 1.2, .height = 1.1, p_hclust)

# end #########################################################################

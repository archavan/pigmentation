# Arun Chavan
# Started: 2025-08-19

# background ==================================================================

# setup =======================================================================
library(tidyverse)
library(here)
library(ggrepel)
library(pvclust)
library(factoextra)
library(glue)

basedir <- here("species/scla/bulk/PIG022/refseq_v1")

resdir <- here(basedir, "results/03_clustering_pca")
fs::dir_create(resdir)

source(here(basedir, "code/utilities_ggplot.R"), local = TRUE)

# data ========================================================================

cts <- read_csv(here(basedir, "data/03_compiled/counts.csv"))
tpm <- read_csv(here(basedir, "data/03_compiled/tpm.csv"))
metadata <- read_csv(here(basedir, "data/samples.csv"))

sample_order <- c("S", "I1", "I2", "I3", "I4", "I5", "I6")

###############################################################################
# distributions ===============================================================
###############################################################################

(cts |>
  pivot_longer(
    cols = where(is.numeric),
    names_to = "sample_id",
    values_to = "counts"
  ) |>
  left_join(metadata, by = "sample_id") |>
  mutate(sample = fct(sample, sample_order)) |>
  ggplot(aes(log10(counts + 1))) +
  geom_histogram(binwidth = 0.2, boundary = 1) +
  facet_grid(cols = vars(sample), rows = vars(animal)) +
  theme_bw_2()) |>
  (function(x) {
    ggsave(
      here(resdir, "distribution_counts_log10.png"),
      x,
      width = 9,
      height = 5
    )
  })()

(tpm |>
  pivot_longer(
    cols = where(is.numeric),
    names_to = "sample_id",
    values_to = "tpm"
  ) |>
  left_join(metadata, by = "sample_id") |>
  mutate(sample = fct(sample, sample_order)) |>
  ggplot(aes(sqrt(tpm))) +
  geom_histogram(binwidth = 0.25, boundary = 0) +
  facet_grid(cols = vars(sample), rows = vars(animal)) +
  coord_cartesian(xlim = c(0, 10)) +
  theme_bw_2()) |>
  (function(x) {
    ggsave(
      here(resdir, "distribution_tpm_sqrt.png"),
      x,
      width = 9,
      height = 5
    )
  })()

(tpm |>
  pivot_longer(
    cols = where(is.numeric),
    names_to = "sample_id",
    values_to = "tpm"
  ) |>
  left_join(metadata, by = "sample_id") |>
  mutate(sample = fct(sample, sample_order)) |>
  ggplot(aes(log10(tpm + 1))) +
  geom_histogram(binwidth = 0.25, boundary = 0) +
  facet_grid(cols = vars(sample), rows = vars(animal)) +
  theme_bw_2()) |>
  (function(x) {
    ggsave(
      here(resdir, "distribution_tpm_log10.png"),
      x,
      width = 9,
      height = 5
    )
  })()

(tpm |>
  pivot_longer(
    cols = where(is.numeric),
    names_to = "sample_id",
    values_to = "tpm"
  ) |>
  left_join(metadata, by = "sample_id") |>
  mutate(sample = fct(sample, sample_order)) |>
  ggplot(aes(tpm)) +
  geom_histogram(binwidth = 2.5, boundary = 0) +
  facet_grid(cols = vars(sample), rows = vars(animal)) +
  coord_cartesian(xlim = c(0, 100)) +
  theme_bw_2()) |>
  (function(x) {
    ggsave(
      here(resdir, "distribution_tpm.png"),
      x,
      width = 9,
      height = 5
    )
  })()

###############################################################################
# clustering ==================================================================
###############################################################################

pv <- tpm |>
  select(where(is.numeric)) |>
  mutate(across(everything(), sqrt)) |>
  pvclust(
    method.hclust = "complete",
    method.dist = "corr",
    nboot = 1000,
    parallel = TRUE
  )

png(here(resdir, "pvclust.png"), width = 8, height = 6, units = "in", res = 300)
plot(pv, print.num = FALSE, print.pv = "bp", cex.pv = 0.7)
dev.off()

###############################################################################
# PCA =========================================================================
###############################################################################

## run pca ====================================================================

run_pca <- function(.tpm) {
  .tpm |>
    column_to_rownames("gene_id") |>
    sqrt() |>
    (function(x) x[which(rowSums(x) > 0), ])() |>
    t() |>
    prcomp(x = _, center = TRUE, scale. = TRUE)
}

add_metadata_to_pc_coords <- function(.pca_out) {
  .pca_out$x |>
    as.data.frame() |>
    rownames_to_column("sample_id") |>
    left_join(y = metadata, by = "sample_id")
}

pca_out <- run_pca(tpm)
pc_space <- add_metadata_to_pc_coords(pca_out)
pc_space <- mutate(pc_space, sample = fct(sample, sample_order))

## plot variance explained ====================================================

var_explained <- get_eig(pca_out) |>
  rownames_to_column("PC") |>
  mutate(PC = gsub("Dim\\.", "PC", PC)) |>
  mutate(variance.percent = round(variance.percent, 1))

var_explained |>
  mutate(PC = str_remove(PC, "PC")) |>
  mutate(PC = as.integer(PC)) |>
  ggplot(aes(PC, variance.percent)) +
  geom_point(size = 3) +
  geom_line() +
  scale_x_continuous(breaks = scales::breaks_pretty()) +
  scale_y_continuous(name = "Variance Explained (%)") +
  theme_bw_2()
ggsave(here(resdir, "pca_var-explained.png"), width = 6, height = 4)

## plot PC score against segment ==============================================
pc_space |>
  pivot_longer(
    cols = starts_with("PC"),
    names_to = "PC",
    values_to = "PC_score"
  ) |>
  filter(PC %in% paste0("PC", 1:12)) |>
  mutate(PC = fct(PC, paste0("PC", 1:12))) |>
  ggplot(aes(sample, PC_score, color = as.factor(rep))) +
  geom_line(aes(group = rep), linewidth = 0.3) +
  geom_point(size = 1.5) +
  facet_wrap(vars(PC), scales = "free") +
  scale_x_discrete(name = "gut segment") +
  scale_y_continuous(name = "PC score") +
  scale_color_discrete(name = "animal") +
  theme_bw_2() +
  theme(
    strip.text = element_text(size = 11)
  )

ggsave(here(resdir, "pca_pc-score-vs-gut-segment.png"), width = 9.5, height = 6)


## plot ordinations ===========================================================

plot_pc_ordination <- function(dim1, dim2) {
  pc_space |>
    arrange(rep, sample) |>
    mutate(rep = as.factor(rep)) |>
    mutate(
      textcolor = case_when(
        sample %in% c("S", "I1", "I2") ~ "white",
        sample %in% c("I3", "I4", "I5", "I6", "I7") ~ "black"
      )
    ) |>
    ggplot(aes(.data[[dim1]], .data[[dim2]])) +
    geom_path(aes(group = animal, linetype = rep), color = "grey50") +
    geom_point(aes(fill = sample), size = 9, shape = 21, stroke = 0) +
    geom_text(aes(label = sample_id, color = textcolor), size = 3) +
    scale_fill_viridis_d(name = "gut segment") +
    scale_color_identity() +
    scale_linetype_manual(
      values = c("solid", "dotted", "dotdash", "longdash"),
      name = "animal"
    ) +
    coord_equal() +
    labs(
      x = glue(
        "{dim1} ({varexp}%)",
        varexp = var_explained$variance.percent[var_explained$PC == dim1]
      ),
      y = glue(
        "{dim2} ({varexp}%)",
        varexp = var_explained$variance.percent[var_explained$PC == dim2]
      )
    ) +
    theme_bw(base_size = 16, base_family = .use_font) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text = element_text()
    )
}

ordinations <- expand_grid(
  dim1 = paste0("PC", 1:8),
  dim2 = paste0("PC", 1:8)
) |>
  filter(dim1 != dim2) |>
  mutate(xy = map2(dim1, dim2, ~ sort(c(.x, .y)))) |>
  filter(!duplicated(xy)) |>
  select(-xy) |>
  mutate(ord = map2(dim1, dim2, plot_pc_ordination))

fs::dir_create(here(resdir, "pca_ordinations"))

for (i in seq_len(nrow(ordinations))) {
  ggsave(
    here(
      resdir,
      "pca_ordinations",
      glue("{x}-{y}.png", x = ordinations$dim1[i], y = ordinations$dim2[i])
    ),
    plot = ordinations$ord[[i]],
    width = 9,
    height = 7
  )
}

# end =========================================================================

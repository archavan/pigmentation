# Arun Chavan
# Started: 2026-04-23

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)
library(glue)
library(ggtext)

basedir <- here("figures")

resdir <- here(basedir, "results/06a_brdu-plots")
fs::dir_create(resdir)

source(here(basedir, "code/utils_themes.R"), local = TRUE)
source(here(basedir, "code/utils_colors.R"), local = TRUE)

###############################################################################
# data ========================================================================
###############################################################################

cv <- read_csv(here(basedir, "data/BrdU-count_CV.csv"))
int <- read_csv(here(basedir, "data/BrdU-count_I1-I4.csv"))

cv <- cv |>
  rowid_to_column("cv_unit") |>
  pivot_longer(
    cols = -cv_unit,
    names_to = c("tissue", "rep"),
    names_sep = "\\.",
    values_to = "count"
  ) |>
  mutate(
    tissue = case_when(
      tissue == "crypt" ~ "Crypt",
      tissue == "v" ~ "Villus",
      .default = NA_character_
    )
  ) |>
  mutate(
    color = case_when(
      tissue == "Crypt" ~ "white",
      tissue == "Villus" ~ "dark",
      .default = NA_character_
    )
  )

int <- int |>
  rename(rep = Animal) |>
  pivot_longer(
    cols = -rep,
    names_to = "tissue",
    values_to = "count"
  ) |>
  mutate(
    color = case_when(
      tissue %in% c("I1", "I3") ~ "white",
      tissue %in% c("I2", "I4") ~ "dark"
    )
  )


###############################################################################
# plots =======================================================================
###############################################################################

cv |>
  summarize(count = mean(count), .by = c(tissue, rep, color)) |>
  ggplot(aes(count, tissue)) +
  geom_jitter(
    aes(fill = color),
    shape = 21,
    size = 1.5,
    position = position_jitter(width = 0, height = 0.3, seed = 2)
  ) +
  stat_summary(fun = mean, geom = "point", shape = "|", size = 1.5) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    linewidth = lwd_thin
  ) +
  scale_fill_manual(
    breaks = names(clrs_pigment_small),
    values = clrs_pigment_small,
    guide = "none"
  ) +
  scale_x_continuous(name = "BrdU<sup>+</sup> nuclei") +
  annotate(geom = "segment", x = 15.2, xend = 15.2, y = 1, yend = 2) +
  annotate(
    geom = "text",
    x = 15.6,
    y = 1.5,
    label = "P = 0.0087",
    angle = 270,
    vjust = 0,
    size = ts_tiny / .pt
  ) +
  coord_cartesian(clip = "off") +
  theme_figure() +
  theme(
    panel.border = element_blank(),
    panel.widths = unit(1.35, "in"),
    panel.heights = unit(0.6, "in"),
    axis.line = element_line(color = "black", linewidth = lwd_thin),
    axis.title.y = element_blank(),
    axis.title.x = element_markdown()
  )

save_quartz(here(resdir, "brdu_cv.pdf"), .width = 1.65, .height = 1)


int |>
  ggplot(aes(count, tissue)) +
  geom_jitter(
    aes(fill = color),
    shape = 21,
    size = 1.5,
    position = position_jitter(width = 0, height = 0.2, seed = 1)
  ) +
  stat_summary(fun = mean, geom = "point", shape = "|", size = 1.5) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    linewidth = lwd_thin
  ) +
  scale_fill_manual(
    breaks = names(clrs_pigment_small),
    values = clrs_pigment_small,
    guide = "none"
  ) +
  scale_x_continuous(
    name = "BrdU<sup>+</sup> nuclei",
    breaks = seq(0, 80, length.out = 5),
    limits = c(0, 80)
  ) +
  coord_cartesian(clip = "off") +
  theme_figure() +
  theme(
    panel.border = element_blank(),
    panel.widths = unit(1.35, "in"),
    panel.heights = unit(1, "in"),
    axis.line = element_line(color = "black", linewidth = lwd_thin),
    axis.title.y = element_blank(),
    axis.title.x = element_markdown()
  )

save_quartz(here(resdir, "brdu_int.pdf"), .width = 1.6, .height = 1.5)

# end #########################################################################

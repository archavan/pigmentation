# Arun Chavan
# Started: 2026-10-05

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)
library(glue)
library(ggtext)
library(Seurat)
library(patchwork)
library(grid)
library(conflicted)

conflicts_prefer(
  dplyr::filter,
  dplyr::select,
  dplyr::rename
)

basedir <- here("figures")

dir_cr <- here("species/crob/bulk/ARCA123/KY21")
dir_cr_genome <- here("species/crob/genome/KY21")

dir_sc <- here("species/scla/bulk/PIG022/refseq_v1")
dir_sc_genome <- here("species/scla/genome/refseq_v1")

dir_mm <- here("species/mmus/zwick")
dir_mm_genome <- here("species/mmus/genome")

resdir <- here(basedir, "results/08a_carb-transporters")
fs::dir_create(resdir)

source(here(basedir, "code/utils_themes.R"), local = TRUE)
source(here(basedir, "code/utils_colors.R"), local = TRUE)
source(here(basedir, "code/utils_seurat-plots.R"), local = TRUE)

###############################################################################
# data ========================================================================
###############################################################################
uss_sc <- read_csv(here(dir_sc, "data/03_compiled/uss_intestine.csv"))

tpm_cr <- read_csv(here(dir_cr, "data/03_compiled/tpm_no-testis.csv"))
tpm_sc <- read_csv(here(dir_sc, "data/03_compiled/tpm_clean.csv")) |>
  filter(gene_id %in% uss_sc$gene_id)


ec_mm <- read_rds(here(dir_mm, "data/enterocytes.rds"))
ec_mm$domain <- factor(ec_mm$domain, LETTERS[1:5])
ec_mm <- SetIdent(ec_mm, value = "domain")

md_cr <- read_csv(here(dir_cr, "data/00_metadata/samples.csv")) |>
  select(tissue_name, tissue_short, sample_name, rep_id, colour) |>
  rename(sample_id = tissue_name, sample = tissue_short) |>
  filter(sample != "S") |>
  mutate(
    pigment = case_when(
      sample %in% c("I2", "I4") ~ "dark",
      sample %in% c("I1", "I3") ~ "white",
      .default = NA_character_
    )
  )

md_sc <- read_csv(here(dir_sc, "data/samples_clean.csv")) |>
  select(sample_id, sample, rep) |>
  filter(sample != "S") |>
  mutate(
    pigment = case_when(
      sample %in% c("I1") ~ "white",
      sample %in% c("I2", "I3", "I4") ~ "dark",
      sample %in% c("I5", "I6") ~ "light",
      .default = NA_character_
    )
  )

fn4_cr <- read_tsv(here(
  dir_cr_genome,
  "results/genes/functional-modules_gut.tsv"
))

fn4_sc <- read_tsv(here(
  dir_sc_genome,
  "results/genes/functional-modules_gut.tsv"
))

fn4_mm <- read_tsv(here(
  dir_mm_genome,
  "gene-modules/functional-modules_gut.tsv"
))


###############################################################################
# prepare and plot: crob and scla =============================================
###############################################################################

## prepare data ===============================================================
dd_cr <- tpm_cr |>
  pivot_longer(cols = -gene_id, names_to = "sample_name", values_to = "tpm") |>
  inner_join(md_cr, by = "sample_name") |>
  inner_join(
    fn4_cr |>
      select(id, label, ortholog, hog_label) |>
      filter(str_detect(label, "Uptake\\: carb")),
    by = c("gene_id" = "ortholog")
  ) |>
  mutate(
    trans_type = case_when(
      str_detect(hog_label, "Slc5") ~ "Active",
      str_detect(hog_label, "Slc2") ~ "Passive",
      .default = NA_character_
    )
  ) |>
  mutate(
    lab = map2_chr(
      hog_label,
      gene_id,
      function(x, y) {
        hog <- gsub("(.+)_(.+)", "\\2", x)
        paralog_number <- gsub("(.+)_(.+)", "\\1", x)
        paste0(hog, "<br>", y, " / ", paralog_number)
      }
    )
  )

dd_sc <- tpm_sc |>
  pivot_longer(cols = -gene_id, names_to = "sample_name", values_to = "tpm") |>
  inner_join(md_sc, by = c("sample_name" = "sample_id")) |>
  inner_join(
    fn4_sc |>
      select(id, label, ortholog, hog_label) |>
      filter(str_detect(label, "Uptake\\: carb")),
    by = c("gene_id" = "ortholog")
  ) |>
  mutate(
    trans_type = case_when(
      str_detect(hog_label, "Slc5") ~ "Active",
      str_detect(hog_label, "Slc2") ~ "Passive",
      .default = NA_character_
    )
  ) |>
  mutate(
    lab = map2_chr(
      hog_label,
      gene_id,
      function(x, y) {
        hog <- gsub("(.+)_(.+)", "\\2", x)
        paralog_number <- gsub("(.+)_(.+)", "\\1", x)
        paste0(hog, "<br>", y, " / ", paralog_number)
      }
    )
  )

## functions ==================================================================

plot_gene <- function(.gene_id, .df) {
  .df |>
    filter(gene_id == .gene_id) |>
    ggplot(aes(sample, tpm)) +
    geom_line(
      data = function(x) {
        x |> summarize(tpm = mean(tpm), .by = c(sample, gene_id))
      },
      mapping = aes(sample, tpm, group = gene_id),
      linewidth = 0.5,
      color = "grey"
    ) +
    geom_point(
      aes(fill = pigment),
      shape = 21,
      size = 1.5,
      stroke = 0.4
    ) +
    scale_fill_manual(values = clrs_pigment_small, guide = "none") +
    scale_y_continuous(limits = range(.df$tpm), name = "TPM") +
    scale_x_discrete(name = "Segment") +
    ggtitle(unique(.df$lab[.df$gene_id == .gene_id])) +
    theme_classic() +
    theme(
      plot.margin = margin(),
      plot.title = element_markdown(lineheight = 1.2, size = 5),
      axis.text = element_text(size = 6),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.title = element_text(size = 6),
      axis.line = element_line(linewidth = 0.25),
      axis.ticks = element_line(linewidth = 0.25)
    )
}

label_grob <- function(.trans_type) {
  grobTree(
    linesGrob(
      x = unit(c(0, 1), "npc"),
      y = unit(c(0, 0), "npc"),
      gp = gpar(col = "grey", lwd = 2)
    ),
    textGrob(
      .trans_type,
      x = unit(0.5, "npc"),
      y = unit(0.6, "npc"),
      gp = gpar(fontsize = 6)
    )
  )
}

title_row <- function(txt, face = 1) {
  ggplot() +
    labs(title = txt) +
    theme_void() +
    theme(
      plot.title = element_text(
        size = 7,
        face = face,
        margin = margin(t = 6, unit = "pt"),
        color = "black"
      ),
      plot.title.position = "panel"
    )
}


## plot crob ==================================================================
pp_active_cr <- unique(dd_cr$gene_id[dd_cr$trans_type == "Active"]) |>
  map(plot_gene, .df = dd_cr) |>
  wrap_plots(nrow = 1, axis_titles = "collect")

pp_passive_cr <- unique(dd_cr$gene_id[dd_cr$trans_type == "Passive"]) |>
  map(plot_gene, .df = dd_cr) |>
  wrap_plots(nrow = 1, axis_titles = "collect")

active_cr <- (wrap_elements(label_grob("Active transporters")) +
  labs(tag = "a") +
  theme(plot.tag = element_text(size = 8, face = 2))) /
  pp_active_cr +
  plot_layout(heights = c(0.1, 0.9))

passive_cr <- (wrap_elements(label_grob("Passive transporters")) +
  labs(tag = "b") +
  theme(plot.tag = element_text(size = 8, face = 2))) /
  pp_passive_cr +
  plot_layout(heights = c(0.1, 0.9))

pp_cr <- wrap_plots(
  title_row("Ciona", face = 4),
  wrap_plots(active_cr, passive_cr, widths = c(3, 2)),
  nrow = 2,
  heights = c(0.01, 1)
)

save_quartz(
  here(resdir, "crob.pdf"),
  .width = 5.5,
  .height = 2.4,
  .plot = pp_cr
)

## plot scla ==================================================================
pp_active_sc <- unique(dd_sc$gene_id[dd_sc$trans_type == "Active"]) |>
  map(plot_gene, .df = dd_sc) |>
  wrap_plots(nrow = 1, axis_titles = "collect")

pp_passive_sc <- unique(dd_sc$gene_id[dd_sc$trans_type == "Passive"]) |>
  map(plot_gene, .df = dd_sc) |>
  wrap_plots(nrow = 1, axis_titles = "collect")

active_sc <- (wrap_elements(label_grob("Active transporters")) +
  labs(tag = "c") +
  theme(plot.tag = element_text(size = 8, face = 2))) /
  pp_active_sc +
  plot_layout(heights = c(0.1, 0.9))

passive_sc <- (wrap_elements(label_grob("Passive transporters")) +
  labs(tag = "d") +
  theme(plot.tag = element_text(size = 8, face = 2))) /
  pp_passive_sc +
  plot_layout(heights = c(0.1, 0.9))

pp_sc <- wrap_plots(
  title_row("Styela", face = 4),
  wrap_plots(active_sc, passive_sc, widths = c(4, 2)),
  nrow = 2,
  heights = c(0.01, 1)
)

save_quartz(
  here(resdir, "scla.pdf"),
  .width = 6.5,
  .height = 2.4,
  .plot = pp_sc
)

###############################################################################
# prepare and plot mouse ======================================================
###############################################################################

carb_trans_genes_mm <- fn4_mm |>
  filter(str_detect(label, "Uptake\\: carb")) |>
  pull(genes) |>
  str_split_1("; ")

## A-P domain =================================================================

dp_mm <- dp_extract_dotplot_data(ec_mm, carb_trans_genes_mm) |>
  mutate(
    trans_type = case_when(
      str_detect(features.plot, "Slc5") ~ "Active transporter",
      str_detect(features.plot, "Slc2") ~ "Passive transporter",
      .default = NA_character_
    )
  ) |>
  mutate(lab = glue("*{features.plot}* ({trans_type})")) |>
  ggplot() +
  geom_point(
    aes(id, lab, size = pct.exp, color = avg.exp.scaled),
    shape = 16,
    stroke = 0
  ) +
  scale_size(
    range = c(0, 4),
    name = "% Cells with Expression",
    guide = guide_legend(
      theme = theme(
        legend.key.height = unit(0.6, "lines"),
        legend.key.width = unit(0.5, "lines")
      )
    )
  ) +
  scale_color_gradient(
    low = "white",
    high = clrs_gradient2$high,
    limits = c(-1.5, 1.5),
    breaks = c(-1, 0, 1),
    oob = scales::squish,
    name = "z(Average Expression)",
    guide = guide_colorbar(
      theme = theme(
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(3, "lines")
      )
    )
  ) +
  scale_y_discrete(labels = function(y) str_trunc(y, 50), position = "right") +
  labs(x = "Domain") +
  coord_fixed() +
  theme_figure() +
  theme(
    plot.tag = element_text(size = 8, face = 2),
    panel.border = element_blank(),
    panel.grid.major = element_line(linewidth = 0.1, color = "grey"),
    axis.line = element_line(linewidth = lwd_thin, color = "black"),
    axis.text.y.right = element_markdown(),
    axis.title.y = element_blank(),
    legend.position = "bottom",
    legend.title.position = "top",
    legend.location = "plot",
    legend.title = element_markdown(
      hjust = 0.5,
      angle = 0,
      size = ts_small
    ),
    strip.clip = "off",
    strip.text = element_markdown(
      hjust = 0,
      vjust = 0,
      margin = margin(b = 1),
    )
  )

## crypt - villus ==============================================================
dd_mm <- FetchData(ec_mm, c("PC_1", carb_trans_genes_mm)) |>
  as.data.frame() |>
  rownames_to_column("barcode") |>
  pivot_longer(
    cols = starts_with("Slc"),
    names_to = "gene",
    values_to = "norm_exp"
  ) |>
  mutate(
    trans_type = case_when(
      str_detect(gene, "Slc5") ~ "Active transporter",
      str_detect(gene, "Slc2") ~ "Passive transporter",
      .default = NA_character_
    )
  ) |>
  mutate(lab = glue("*{gene}* ({trans_type})")) |>
  mutate(lab = factor(lab, levels = unique(lab[order(trans_type, gene)])))

pc_mm <- dd_mm |>
  filter(norm_exp > 0) |>
  ggplot(aes(PC_1, norm_exp)) +
  geom_point(shape = 16, stroke = 0, size = 0.6, color = "grey", alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 0.5) +
  facet_grid(cols = vars(lab), axes = "all") +
  labs(x = "PC1 (Crypt to Villus)", y = "Normalized Expression", tag = "f") +
  theme_classic() +
  theme(
    plot.tag = element_text(size = 8, face = 2),
    strip.background = element_blank(),
    strip.text = element_markdown(hjust = 0, size = 6),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 6),
    axis.line = element_line(linewidth = 0.25),
    axis.ticks = element_line(linewidth = 0.25)
  )

pp_mm <- wrap_plots(
  title_row("Mouse", face = 2),
  wrap_plots(
    wrap_elements(panel = dp_mm) +
      labs(tag = "e") +
      theme(plot.tag = element_text(size = 8, face = 2)),
    pc_mm,
    nrow = 1,
    widths = c(1, 1.5)
  ),
  nrow = 2,
  heights = c(0.01, 1)
)

save_quartz(
  here(resdir, "mmus.pdf"),
  .width = 5.6,
  .height = 2.25,
  .plot = pp_mm
)

###############################################################################
# combine =====================================================================
###############################################################################

composed <- wrap_plots(
  wrap_plots(pp_cr, plot_spacer(), nrow = 1, widths = c(5.5, 1)),
  wrap_plots(pp_sc),
  wrap_plots(pp_mm, plot_spacer(), nrow = 1, widths = c(5.8, 0.7)),
  nrow = 3,
  heights = c(2.4, 2.4, 2.2)
)

save_quartz(
  here(resdir, "composed.pdf"),
  .width = 6.5,
  .height = 6.6,
  .plot = composed
)

# end #########################################################################

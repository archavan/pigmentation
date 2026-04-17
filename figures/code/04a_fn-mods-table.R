# Arun Chavan
# Started: 2026-04-16

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)
library(glue)
library(flextable)
library(officer)

basedir <- here("figures")

resdir <- here(basedir, "results/04a_fn-mods-table")
fs::dir_create(resdir)

###############################################################################
# data ========================================================================
###############################################################################

fn4 <- read_tsv(here(basedir, "data/fn4-modules-for-table.tsv"))

###############################################################################
# make table ==================================================================
###############################################################################

set_flextable_defaults(
  font.size = 7,
  font.family = "Helvetica",
  font.color = "#000000",
  table.layout = "fixed",
  padding.top = 2,
  padding.bottom = 2,
  padding.left = 2,
  padding.right = 2
)

alphabetize_genes <- function(x, sep = ", ") {
  str_split_1(x, sep) |> sort() |> str_flatten(collapse = sep, na.rm = TRUE)
}

fn4_reorg <- fn4 |>
  select(module.1, module.2, module.3, module.4, genes) |>
  rename_with(~ str_replace(.x, "module", "Hierarchy")) |>
  mutate(across(starts_with("Hierarchy"), function(x) {
    str_replace_all(x, "\\.", " ") |> str_to_sentence()
  })) |>
  mutate(genes = str_replace_all(genes, ";\\s*", ", ")) |>
  mutate(genes = map_chr(genes, alphabetize_genes))

ft <- flextable(fn4_reorg)

ft <- ft |>
  merge_v(j = ~ Hierarchy.1 + Hierarchy.2 + Hierarchy.3 + Hierarchy.4) |>
  set_table_properties(layout = "autofit") |>
  set_header_labels(
    Hierarchy.1 = "First",
    Hierarchy.2 = "Second",
    Hierarchy.3 = "Third",
    Hierarchy.4 = "Fourth",
    genes = "Genes"
  ) |>
  add_header_row(
    values = c("Functional Hierarchy Level", "Genes"),
    colwidths = c(4, 1)
  ) |>
  merge_v(j = 5, part = "header") |>
  style(
    j = 5,
    pr_t = fp_text_default(italic = TRUE)
  ) |>
  style(
    j = 1,
    pr_t = fp_text_default(bold = TRUE)
  ) |>
  bg(
    i = seq(2, nrow(fn4_reorg), 2),
    j = ~genes,
    bg = "#EEEEEE"
  ) |>
  theme_box() |>
  align(i = 1, part = "header", align = "center")


save_as_docx(
  ft,
  path = here(resdir, "fn4-modules-table.docx"),
  pr_section = prop_section(page_size = page_size(width = 8.5, height = 11))
)

system2("open", here(resdir, "fn4-modules-table.docx"))

# end #########################################################################

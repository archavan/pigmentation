# Arun Chavan
# Started: 2025-11-29

# setup =======================================================================
library(tidyverse)
library(here)
library(glue)
library(Seurat)
library(ggrepel)

basedir <- here("species/mmus/zwick")
genomedir <- here("species/mmus/genome")

resdir <- here(basedir, "results/functional-modules")
fs::dir_create(resdir)

source(here(basedir, "code/utilities_ggplot.R"))
source(here(basedir, "code/utilities_seurat-plots.R"))

# data ========================================================================
ec <- read_rds(here(basedir, "data/enterocytes.rds"))
ec$domain <- factor(ec$domain, LETTERS[1:5])
ec <- SetIdent(ec, value = "domain")


ps <- read_csv(here(basedir, "data/enterocytes_agg-by-domain.csv"))
ps_long <- ps |>
  pivot_longer(
    cols = -gene,
    names_to = "domain",
    values_to = "norm_expr"
  ) |>
  mutate(
    domain = factor(domain, levels = LETTERS[1:5])
  )

# genomic data ================================================================
modules <- read_tsv(here(genomedir, "gene-modules/functional-modules_gut.tsv"))

modules$module <- NA_character_
for (i in seq_len(nrow(modules))) {
  modules$module[i] <- modules[i, 2:5] |>
    unlist() |>
    na.omit() |>
    paste0(collapse = " - ")
}

modules <- modules |>
  separate_longer_delim(genes, "; ")

# plotting setup ==============================================================
plot_dotplot2 <- function(
  .features,
  .title,
  .reorder.genes = TRUE,
  .wrap.title = TRUE
) {
  if (.wrap.title) {
    .title <- str_replace_all(.title, " - ", " → \n  ")
  }

  dd <- dp_extract_dotplot_data(ec, .features)
  if (.reorder.genes) {
    dd <- dp_define_feature_order(dd)
  }
  dp_split_dotplot_data(dd) |>
    purrr::map(dp_plot) |>
    purrr::map(dp_add_title, .title) |>
    purrr::map(function(p) {
      p +
        theme(plot.title = element_text(size = 8, face = 1))
    })
}

###############################################################################
# plot functional modules =====================================================
###############################################################################

mod_nest <- modules |>
  select(module, genes) |>
  filter(!is.na(genes)) |>
  distinct() |>
  group_by(module) |>
  summarize(genes = list(genes)) |>
  ungroup()

## plot individual genes per module ===========================================
### dotplots

mod_nest <- mod_nest |>
  mutate(dp = map2(genes, module, plot_dotplot2))

fs::dir_create(here(resdir, "module-dotplots"))
for (i in seq_len(nrow(mod_nest))) {
  pdf(
    file = here(resdir, "module-dotplots", paste0(mod_nest$module[i], ".pdf")),
    width = 4,
    height = 10,
    onefile = TRUE
  )
  print(mod_nest$dp[[i]])
  dev.off()
}

### pseudobulk

fs::dir_create(here(resdir, "module-pseudobulk"))

plot_and_save_psedobulk_plot <- function(.module.idx, .save = FALSE) {
  genes_in_module <- mod_nest$genes[[.module.idx]]

  p <- ps_long |>
    filter(gene %in% genes_in_module) |>
    ggplot(aes(x = domain, y = norm_expr)) +
    geom_line(aes(group = gene)) +
    geom_point() +
    facet_wrap(vars(gene), scales = "free_y") +
    labs(
      title = str_replace_all(mod_nest$module[[.module.idx]], " - ", " → \n  "),
      x = "Enterocyte Domain",
      y = "Normalized Expression"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 8, face = 1)
    )

  q <- ggplot_build(p)
  n_facet_rows <- max(summarise_layout(q)$row)
  n_facet_cols <- max(summarise_layout(q)$col)

  if (.save) {
    ggsave(
      filename = here(
        resdir,
        "module-pseudobulk",
        paste0(mod_nest$module[[.module.idx]], ".pdf")
      ),
      plot = p,
      width = n_facet_cols * 2.5,
      height = n_facet_rows * 2 + 1
    )
  } else {
    return(p)
  }
}


for (i in seq_len(nrow(mod_nest))) {
  plot_and_save_psedobulk_plot(i, .save = TRUE)
}

# aggregated (mean for all genes in module) ===================================

modules_grp <- modules |>
  select(module, genes) |>
  filter(!is.na(genes)) |>
  distinct() |>
  group_by(module) |>
  summarize(genes = list(genes)) |>
  ungroup() |>
  mutate(
    normd = map(genes, function(x) {
      FetchData(ec, x, layer = "data") |>
        rowMeans()
    })
  )

clust_ids <- FetchData(ec, "domain") |>
  as.data.frame() |>
  rownames_to_column("barcode")

modules_grp |>
  select(module, normd) |>
  unnest_wider(col = normd) |>
  pivot_longer(cols = -module, names_to = "barcode", values_to = "normd") |>
  left_join(clust_ids, by = "barcode") |>
  group_by(module, domain) |>
  summarize(normd = mean(normd)) |>
  ungroup() |>
  group_by(module) |>
  mutate(z = scale(normd)[, 1]) |>
  ungroup() |>
  mutate(domain = factor(domain, levels = LETTERS[1:5])) |>
  ggplot(aes(domain, module, fill = z)) +
  geom_tile() +
  scale_fill_viridis_c(direction = -1, name = "Z-scaled\naggr. expr.") +
  scale_y_discrete(position = "right") +
  coord_fixed() +
  theme_tileplot() +
  theme(
    legend.position = "left"
  )

ggsave(
  here(resdir, "all-aggregated.pdf"),
  width = 10,
  height = 7
)

# end =========================================================================

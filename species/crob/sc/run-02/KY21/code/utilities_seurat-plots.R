# Helper functions for making plots from a Seurat object with custom setup. 

# Dependencies:
#   - Packages
#     - Seurat
#     - tidyverse
#     - glue
#   - Other scripts:
#     - ./utilities_ggplot.R
#   - data
#     - here(genomedir, "genes/gene_display-names.tsv")

###############################################################################
# Split vector into chunks ####
###############################################################################

#' Split vector into multiple small vectors
#'
#' This is used internally in plot_dotplot() to split long gene lists into chunks
#'
#' @param vec vector
#' @param chunk_size length of chunks
#'
#' @return a list

split_vector <- function(vec, chunk_size) {
  split(vec, ceiling(seq_along(vec)/chunk_size))
}

###############################################################################
# extract dotplot data ####
###############################################################################

#' extract dotplot data
#'
#' @param .obj Seurat object
#' @param .features vector of gene ids
#' @param assay assay of seurat object to use
#'
#' @return data frame

dp_extract_dotplot_data <- function(.obj, .features, assay = DefaultAssay(.obj)) {
  p <- DotPlot(object = .obj, assay = assay, features = .features)
  dat <- p$data
  ## remove genes with zero expression ===
  # These genes have NaN in avg.exp.scaled since we can't calculate their
  # variance. It's genrally fine to keep them around, but they cause at least two
  # problems: (1) they cause hclust to fail, so at least we have to remove them
  # from hclust call and add their levels back; (2) as long as they are with
  # other non-NaN features, they plot fine, but when all genes in a plot end up
  # being NaN (let's say due to splitting of feature list below), ggplot fails.
  # We can technically take care of this by dynamically determining chunk size in
  # split_vector() so that the last chunk doesn't end up with all NaN genes, but
  # at this point it's not worth the effort. For now we will just get rid of
  # these genes
  dat <- filter(dat, !is.nan(avg.exp.scaled))
  return(dat)
}

###############################################################################
# define cluster order ####
###############################################################################

#' Change factor levels for clusters to reorder them 
#'
#' @param p.dat dotplot data extracted from Seurat::DotPlot() output
#' @param .cluster_order vector specifying cluster order
#'
#' @return dotplot data with releveled cluster order

dp_define_cluster_order <- function(p.dat, .cluster_order) {
  p.dat %>% 
    mutate(id = fct(x = as.character(id), 
                    levels = as.character(.cluster_order)))
} 

###############################################################################
# define feature order ####
###############################################################################

#' Change factor level for genes to reorder them 
#'
#' @param p.dat dotplot data extracted from Seurat::DotPlot() output
#'
#' @return dotplot data with releveled gene order

dp_define_feature_order <- function(p.dat) {
  .feature_order <- p.dat %>% 
    mutate(combined.exp = pct.exp * avg.exp.scaled) %>% 
    group_by(features.plot) %>% 
    slice_max(combined.exp, n = 1, with_ties = FALSE) %>% 
    ungroup() %>% 
    arrange(id) %>% 
    group_by(id) %>% 
    arrange(desc(combined.exp), .by_group = TRUE) %>% 
    ungroup() %>% 
    pull(features.plot) %>% 
    as.character() %>% 
    rev()
  p.dat %>% 
    mutate(features.plot = fct(x = as.character(features.plot),
                               levels = .feature_order))
}

###############################################################################
# split dotplot data ####
###############################################################################

#' split dotplot data by number of genes
#'
#' @param p.dat dotplot data 
#' @param .n number of genes per plot
#'
#' @return list of data dfs for dotplots

dp_split_dotplot_data <- function(p.dat, .n = 50) {
  stopifnot(is.factor(p.dat$features.plot))
  stopifnot(is.factor(p.dat$id))
  gene_levels <- levels(p.dat$features.plot)
  ident_levels <- levels(p.dat$id)
  p.dat <- p.dat %>% 
    mutate(features.plot = as.character(features.plot),
           id = as.character(id))
  genes_split <- split_vector(rev(gene_levels), .n)
  p.dat_split <- purrr::map(genes_split, 
                            ~ p.dat %>% 
                              filter(features.plot %in% .x) %>% 
                              mutate(features.plot = fct(features.plot, 
                                                         levels = rev(.x)),
                                     id = fct(id, levels = ident_levels)))
  return(p.dat_split)
  # the two rev() calls above are for keeping later gene levels in early chunk
  # after splitting the vector
}

###############################################################################
# create dotplot ####
###############################################################################

#' create dotplot
#'
#' @param .df dotplot data 
#'
#' @return gg

dp_plot <- function(.df) {
  ggplot(.df, aes(id, features.plot)) +
    geom_point(aes(color = avg.exp.scaled, size = pct.exp), 
               shape = 16, 
               stroke = 0) +
    scale_size(range = c(0, 6)) +
    scale_color_gradient(low = "white", high = pal_okabe[5]) +
    coord_fixed() +
    theme_half_open_grid(font_size = 10) +
    theme(axis.text.x = element_text(angle = 90, 
                                     hjust = 1, 
                                     vjust = 0.5),
          axis.title = element_blank(),
          panel.grid.major = element_line(linewidth = 0.1))
}

###############################################################################
# add title ####
###############################################################################

#' add title to dotplot
#'
#' @param .gg gg object
#' @param .title character string
#'
#' @return gg object
#' 
dp_add_title <- function(.gg, .title) {
  .gg + 
    ggtitle(.title)
}

###############################################################################
# add gene display name ####
###############################################################################

#' add gene display names to dotplot data
#'
#' @param .gg ggplot object
#' @param .display_names named (gene_id as name)character vector of display names
#'
#' @return data frame

dp_replace_gene_id_with_display_name <- function(.gg, 
                                                 .display_names,
                                                 .trunc_width = 90) {
  .gg +
    scale_y_discrete(labels = function(x) str_trunc(.display_names[x], .trunc_width, "right"),
                     position = "right") +
    theme(legend.position = "left")
}

###############################################################################
# add tissue of origin info ####
###############################################################################

#' add tissue origin info
#'
#' @param .gg ggplot object
#' @param .tissue_of_origin named vector where names are clusters and values are tissue of origin
#'
#' @return ggplot object

dp_add_tissue_of_origin <- function(.gg, .tissue_of_origin) {
  .gg +
    scale_x_discrete(label = function(x) {
      x <- as.character(x)
      glue("({.tissue_of_origin[x]}) {x}")
  })
}

###############################################################################
# wrapper to make dotplot ####
###############################################################################

#' wrapper to make dotplots
#'
#' @param .obj seurat object
#' @param .title plot title
#' @param .features vector of features, which will be split into chunks of 50
#'
#' @return a list of ggplot objects

plot_dotplot <- function(.obj, 
                         .features, 
                         .title = NULL,
                         .cluster_order,
                         .display_names = gene_display_names,
                         .tissue_of_origin = tissue_origin) {
  dp_extract_dotplot_data(.obj = .obj, .features = .features) %>% 
    dp_define_cluster_order(.cluster_order = .cluster_order) %>%
    dp_define_feature_order() %>% 
    dp_split_dotplot_data() %>% 
    purrr::map(dp_plot) %>% 
    purrr::map(~ dp_replace_gene_id_with_display_name(.gg = .x, .display_names = .display_names)) %>% 
    purrr::map(~ dp_add_tissue_of_origin(.gg = .x, .tissue_of_origin = .tissue_of_origin)) %>% 
    purrr::map(dp_add_title, .title = .title)
}

###############################################################################
# save list of dotplots to multipage pdf ####
###############################################################################

#' save dotplots to pdf using cairo_pdf
#'
#' @param p plot (can be a list of plots, which will be saved to multipage pdf)
#' @param filename filename
#' @param width width
#' @param height height
#'
#' @return pdf written to disk

save_dotplot <- function(p, 
                         filename, 
                         width = 13, 
                         height = 11) {
  
  cairo_pdf(filename, 
            width = width, 
            height = height, 
            onefile = TRUE)
  print(p)
  dev.off()
  
}

# end =========================================================================

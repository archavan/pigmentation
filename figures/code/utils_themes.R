# text sizes (pts) ============================================================
ts_tiny <- 5
ts_small <- 5.5
ts_base <- 6
ts_medium <- 6.5
ts_large <- 7
ts_huge <- 8

# -- line widths ==============================================================
lwd_thin <- 0.25
lwd_medium <- 0.5
lwd_thick <- 1

# font ========================================================================
font_family <- "Helvetica"

# theme =======================================================================
theme_figure <- function(
  family = font_family,
  base_size = ts_base,
  ...
) {
  theme_bw(
    base_size = base_size,
    base_family = family,
    base_line_size = lwd_thin,
    base_rect_size = lwd_thin
  ) %+replace%
    theme(
      # blank backgrounds (no white boxes in Illustrator)
      plot.background = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank(),
      strip.background = element_blank(),

      # text sizes
      text = element_text(family = family, size = ts_base),
      plot.title = element_text(size = ts_large, hjust = 0),
      plot.subtitle = element_text(size = ts_medium),
      axis.title = element_text(size = ts_base),
      axis.text = element_text(size = ts_base, color = "black"),
      legend.text = element_text(size = ts_base),
      legend.title = element_text(size = ts_base),
      strip.text = element_text(size = ts_base),

      # lines & ticks
      axis.ticks = element_line(linewidth = lwd_thin),

      complete = TRUE
    )
}

# save ========================================================================
save_quartz <- function(
  .filename,
  .width,
  .height,
  .plot = get_last_plot(),
  ...
) {
  quartz(
    type = "pdf",
    file = .filename,
    width = .width,
    height = .height,
    ...
  )
  print(.plot)
  dev.off()
}

###############################################################################

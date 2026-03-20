# pigment colors ==============================================================
clrs_pigment <- list(
  dark = "#fac597",
  light = "#faeba7",
  white = "#ffffff"
)

clrs_pigment_small <- lapply(clrs_pigment, function(x) {
  ifelse(x == "#ffffff", x, colorspace::darken(x, amount = 0.2))
})

clrs_gradient2 <- list(low = "#578dd3", mid = "white", high = "#dd6e68")
# #4575b4 #d73027

###############################################################################

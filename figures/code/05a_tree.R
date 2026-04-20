# Arun Chavan
# Started: 2026-04-18

###############################################################################
# setup =======================================================================
###############################################################################

library(ape)
library(here)
library(glue)
library(tidyverse)
library(ggtree)
library(conflicted)

conflicts_prefer(dplyr::filter, dplyr::select)

basedir <- here("figures")
resdir <- here(basedir, "results/05a_tree")
fs::dir_create(here(resdir))

source(here(basedir, "code/utils_themes.R"))

###############################################################################
# data ========================================================================
###############################################################################

# Delsuc, F., Philippe, H., Tsagkogeorga, G. et al. A phylogenomic framework and
# timescale for comparative studies of tunicates. BMC Biol 16, 39 (2018).
# https://doi.org/10.1186/s12915-018-0499-2

tr <- read.tree(
  "https://raw.githubusercontent.com/psimion/SuppData_Delsuc_BMCBiol_2018_Dating_Tunicata/refs/heads/master/Delsuc_BMCBiol_2018_Fig1_PhyloBayes_CAT-GTR-G.tree"
)

# drop non-chordates
tr <- drop.tip(
  tr,
  c(
    "Xenoturbella_bocki",
    "Saccoglossus_kowalevskii",
    "Ptychodera_flava",
    "Patiria_miniata",
    "Ophionotus_victoriae",
    "Holothuria_glaberrima",
    "Apostichopus_japonicus",
    "Eucidaris_tribuloides",
    "Strongylocentrotus_purpuratus",
    "Lytechinus_variegatus"
  )
)

# prune vertebrates
tr <- drop.tip(
  tr,
  c(
    "Anas_platyrhynchos",
    "Taeniopygia_guttata",
    "Crocodylus_niloticus",
    "Sphenodon_punctatus",
    "Sarcophilus_harrisii",
    "Cyclorana_alboguttata",
    "Eptatretus_burgeri",
    "Scyliorhinus_canicula",
    "Phrynops_hilarii",
    "Pelodiscus_chinensis",
    "Pelodiscus_sinensis",
    "Chrysemys_picta_bellii",
    "Chelonoidis_nigra",
    "Notophthalmus_viridescens",
    "Ambystoma_mexicanum",
    "Rana_chensinensis",
    "Silurana_tropicalis",
    "Eublepharis_macularius",
    "Python_regius",
    "Protopterus_annectens",
    "Callorhinchus_milii",
    "Oreochromis_niloticus"
  )
)

# replace Styela plicata with clava
tr$tip.label[tr$tip.label == "Styela_plicata"] <- "Styela_clava"


###############################################################################
# plot ========================================================================
###############################################################################

nodes_to_label <- c(
  tunicates = getMRCA(tr, c("Ciona_robusta", "Oikopleura_dioica")),
  vertebrates = getMRCA(tr, c("Petromyzon_marinus", "Mus_musculus")),
  chordates = getMRCA(tr, c("Branchiostoma_floridae", "Mus_musculus"))
)

ggtree(tr, layout = "rectangular", size = 0, ladderize = TRUE) +
  geom_highlight(
    node = nodes_to_label[["tunicates"]],
    fill = "#FFFACD",
    alpha = 0.5
  ) +
  geom_tree(size = lwd_thin) + # redraws edges on top of the highlight
  geom_tiplab(
    aes(
      label = gsub("_", " ", label),
      fontface = ifelse(
        label %in%
          c(
            "Ciona_robusta",
            "Styela_clava",
            "Mus_musculus",
            "Halocynthia_roretzi"
          ),
        "bold.italic",
        "italic"
      )
    ),
    geom = "text",
    size = ts_small / .pt,
    offset = 0
  ) +
  geom_nodelab(
    aes(
      subset = (node == !!nodes_to_label[["chordates"]]),
      label = "Chordates"
    ),
    angle = 90,
    hjust = 0.5,
    vjust = 0,
    nudge_x = -0.1,
    nudge_y = 00,
    size = ts_small / .pt
  ) +
  geom_nodelab(
    aes(
      subset = (node == !!nodes_to_label[["tunicates"]]),
      label = "Tunicates"
    ),
    angle = 90,
    hjust = 0,
    vjust = 0,
    nudge_x = -0.1,
    nudge_y = 0.2,
    size = ts_small / .pt
  ) +
  geom_nodelab(
    aes(
      subset = (node == !!nodes_to_label[["vertebrates"]]),
      label = "Vertebrates"
    ),
    angle = 90,
    hjust = 0,
    vjust = 0,
    nudge_x = -0.1,
    nudge_y = 0.2,
    size = ts_small / .pt
  ) +
  scale_x_continuous(limits = c(NA, 7)) +
  coord_cartesian(clip = "off")

save_quartz(here(resdir, "tree.pdf"), .width = 2.7, .height = 3.1)

# end #########################################################################

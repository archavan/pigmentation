# Arun Chavan
# Started: 2025-03-22

# background ==================================================================

# Here we generate plots for various categories of candidate genes.

#'#############################################################################
# setup =======================================================================
#'#############################################################################

library(GO.db)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(tidyverse)
library(Seurat)
library(here)
library(ggrepel)
library(glue)

conflicted::conflicts_prefer(
  dplyr::filter(),
  dplyr::select(),
  purrr::map(),
  purrr::reduce(),
  base::union(),
  base::intersect()
)

## define paths ===============================================================

basedir <- here("species/crob/bulk/ARCA123/KY21")
genomedir <- here("species/crob/genome/KY21")
resdir <- here(basedir, "results/13a_candidate-genes")
fs::dir_create(resdir)

source(here(basedir, "code/utilities_ggplot.R"))

#'#############################################################################
# data ========================================================================
#'#############################################################################

## rna-seq data ===============================================================

tpm <- read_csv(here(basedir, "data/03_compiled/tpm_no-testis.csv"))
samples <- read_csv(here(basedir, "data/00_metadata/samples.csv"))

highly_expressed_genes <- tpm %>%
  column_to_rownames("gene_id") %>%
  rowSums() %>%
  .[. > 5] %>%
  names()

## genomic data ===============================================================

genes <- read_tsv(here(genomedir, "results/genes/gene-display-names.tsv"))

go_emapper <- read_tsv(here(
  genomedir,
  "results/emapper/processed/emapper_go.tsv"
))

go_prost <- read_tsv(here(
  genomedir,
  "results/prost/processed/go.tsv.gz"
))

mgenes <- read_tsv(here(
  "orthofinder/results/2025-03-06/id-name-mapping/mmus.tsv"
)) %>%
  select(gene_id, symbol) %>%
  distinct()

mmcr <- read_tsv(here(
  "orthofinder/results/2025-03-06/pairwise-orthologs/mmus_vs_crob.tsv"
))

mmcr <- left_join(mmcr, mgenes, by = c("mmus" = "gene_id"))
mmcr <- left_join(
  mmcr,
  select(genes, gene_id, display_name),
  by = c("crob" = "gene_id")
)

#'#############################################################################
# functions for getting gene lists ============================================
#'#############################################################################

get_mm_genes_for_GO_and_offspring <- function(.go_id) {
  offspring_terms <- unique(c(
    GOBPOFFSPRING[[.go_id]],
    GOMFOFFSPRING[[.go_id]],
    GOCCOFFSPRING[[.go_id]]
  ))

  all_go_terms <- unique(c(.go_id, offspring_terms))

  entrez_ids <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = all_go_terms,
    keytype = "GO",
    columns = "ENTREZID"
  )$ENTREZID %>%
    unique() %>%
    na.omit()

  gene_symbols <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = entrez_ids,
    keytype = "ENTREZID",
    columns = "SYMBOL"
  )$SYMBOL

  return(unique(gene_symbols))
}

get_cr_genes_for_GO_and_offspring <- function(.go_id, .go_db = go_emapper) {
  offspring_terms <- unique(c(
    GOBPOFFSPRING[[.go_id]],
    GOMFOFFSPRING[[.go_id]],
    GOCCOFFSPRING[[.go_id]]
  ))
  all_go_terms <- unique(c(.go_id, offspring_terms))

  .go_db %>%
    filter(go_id %in% all_go_terms) %>%
    pull(gene_id) %>%
    unique()
}

get_cr_orthologs <- function(symbols) {
  mmcr %>%
    filter(symbol %in% symbols) %>%
    pull(crob) %>%
    unique()
}

#'#############################################################################
# plotting setup ==============================================================
#'#############################################################################

tissue_order <- c("stm.1.s", "int.1.w", "int.2.y", "int.3.w", "int.4.y")

tpm_long <- tpm %>%
  pivot_longer(cols = -gene_id, names_to = "sample_name", values_to = "tpm") %>%
  left_join(samples) %>%
  mutate(tissue = factor(tissue, tissue_order))

plot_genes_faceted <- function(
  .gene_ids,
  .title,
  .order_by_expression = FALSE,
  .nrow = 8,
  .ncol = 5
) {
  # create plot data
  plot_data <- tpm_long %>%
    filter(gene_id %in% .gene_ids) %>%
    left_join(genes) %>%
    mutate(
      display_name = paste0(
        "[",
        gene_id,
        "] ",
        symbol,
        "\n",
        str_trunc(description, 72)
      )
    )

  # order genes
  if (.order_by_expression) {
    gene_order <- plot_data %>%
      summarise(tpm = mean(tpm), .by = display_name) %>%
      arrange(-tpm) %>%
      pull(display_name)
  } else {
    gene_order <- plot_data %>%
      select(gene_id, display_name) %>%
      distinct() %>%
      mutate(gene_id = factor(gene_id, levels = .gene_ids)) %>%
      arrange(gene_id) %>%
      pull(display_name)
  }

  plot_data <- mutate(
    plot_data,
    display_name = factor(display_name, gene_order)
  )
  plot_data <- arrange(plot_data, display_name)

  # plot
  ggplot(plot_data, aes(x = tissue)) +
    geom_jitter(aes(y = tpm), width = 0.1, color = "grey") +
    stat_summary(aes(y = tpm, group = 1), fun = mean, geom = "line") +
    stat_summary(aes(y = tpm), fun = mean, geom = "point") +
    facet_wrap(
      vars(display_name),
      scales = "free_y",
      nrow = .nrow,
      ncol = .ncol
    ) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    theme_bw_2() +
    theme(
      strip.text = element_text(size = 6, hjust = 0),
      plot.caption = element_text(hjust = 0)
    ) +
    ggtitle(.title)
}

save_plot <- function(p, name, h = 11, w = 9) {
  ggsave(
    here(resdir, paste0(name, ".pdf")),
    p,
    width = w,
    height = h,
    device = cairo_pdf
  )
}

ncol_used <- function(.gene_ids, .ncol_max) {
  ngenes <- length(.gene_ids)
  ifelse(ngenes < .ncol_max, ngenes, .ncol_max)
}

nrow_used <- function(.gene_ids, .ncol) {
  ngenes <- length(.gene_ids)
  ceiling(ngenes / .ncol)
}

calc_plot_dims <- function(.gene_ids, .ncol_max = 5) {
  ncol_used <- min(length(.gene_ids), .ncol_max)
  nrow_used <- ceiling(length(.gene_ids) / ncol_used)
  h <- (nrow_used * 1.5) + 0.5
  w <- (ncol_used * 1.75) + 0.5
  list(h = h, w = w)
}

plot_genes_and_save <- function(
  .gene_ids,
  .title,
  .nrow = 8,
  .ncol = 5,
  ...
) {
  pdims <- calc_plot_dims(.gene_ids, .ncol)
  p <- plot_genes_faceted(
    .gene_ids,
    .title,
    .nrow = .nrow,
    .ncol = .ncol,
    ...
  )
  save_plot(p, .title, pdims$h, pdims$w)
}

#'#############################################################################
# plots =======================================================================
#'#############################################################################

## genes from Zwick et al ======================================================

# fmt:skip
zwick_genes <- list(
  zoneA = c("Meis2", "Otop3", "Sult6b2", "Foxq1", "Fth1"),
  zoneB = c("Fabp1", "Ces1f", "Papss2", "Adh1", "Aadac", "Aldh1a1"),
  zoneC = c("Fam213a", "Sis", "Slc2a2", "Apoc3", "Glud1", ""),
  zoneD = c("Fam151a", "Duoxa2", "Anxa4", "Slc15a1", "Pla2g4c", "Plb1"),
  zoneE = c("Slc10a2", "Suox", "Cd74", "Slc25a48", "Hmgcs1", "Fabp6")
)

zwick_genes_orth <- map(zwick_genes, get_cr_orthologs)

zwick_genes_plots <- vector("list", length(zwick_genes_orth))
names(zwick_genes_plots) <- names(zwick_genes_orth)

for (i in names(zwick_genes_orth)) {
  zwick_genes_plots[[i]] <- plot_genes_faceted(
    zwick_genes_orth[[i]],
    .title = i,
    TRUE
  )
  save_plot(
    zwick_genes_plots[[i]],
    name = paste0("zwick_", i),
    h = calc_plot_dims(zwick_genes_orth[[i]], 5)$h,
    w = calc_plot_dims(zwick_genes_orth[[i]], 5)$w
  )
}

## Fabp genes =================================================================
fabp_genes <- genes %>%
  filter(
    gene_id %in%
      union(
        mmcr %>%
          filter(str_detect(symbol, "Fabp")) %>%
          pull(crob),
        genes %>%
          filter(str_detect(
            display_name,
            regex("fatty acid binding|fabp", TRUE)
          )) %>%
          pull(gene_id)
      )
  ) %>%
  pull(gene_id) %>%
  unique()

plot_genes_and_save(fabp_genes, "fabp", .order_by_expression = TRUE)

## hormones ===================================================================
go_hormones <- "GO:0005179"
hormones_mm <- get_mm_genes_for_GO_and_offspring(go_hormones)

list(
  get_cr_orthologs(hormones_mm),
  get_cr_genes_for_GO_and_offspring(go_hormones, go_emapper)
) %>%
  unlist() %>%
  unique() %>%
  plot_genes_and_save("hormones")

## neuropeptide hormone activity ==============================================
go_neuropeptides <- "GO:0005184"
neuropeptides_mm <- get_mm_genes_for_GO_and_offspring(go_neuropeptides)

neuropeptides_cr <- list(
  get_cr_orthologs(neuropeptides_mm),
  get_cr_genes_for_GO_and_offspring(go_neuropeptides, go_prost)
) %>%
  unlist() %>%
  unique() %>%
  .[. %in% expressed_genes]

plot_genes_and_save(
  neuropeptides_cr[1:100],
  "neuropeptides_1",
  .order_by_expression = TRUE,
  .ncol = 5,
  .nrow = 20
)

plot_genes_and_save(
  neuropeptides_cr[101:200],
  "neuropeptides_2",
  .order_by_expression = TRUE,
  .ncol = 5,
  .nrow = 20
)

plot_genes_and_save(
  neuropeptides_cr[201:300],
  "neuropeptides_3",
  .order_by_expression = TRUE,
  .ncol = 5,
  .nrow = 20
)

## endocrine genes ============================================================
go_endocrine_processing <- "GO:0016486"
peptide_hormone_processing_mm <- get_mm_genes_for_GO_and_offspring(
  go_endocrine_processing
)

endocrine_mm <- c(
  peptide_hormone_processing_mm,
  "Tph1",
  "Chga",
  "Chgb",
  "Neurog3"
)

endocrine_cr <- list(
  get_cr_orthologs(endocrine_mm),
  get_cr_genes_for_GO_and_offspring(go_endocrine_processing)
) %>%
  unlist() %>%
  unique()

plot_genes_and_save(
  endocrine_cr,
  "endocrine",
  .order_by_expression = TRUE,
  .ncol = 5,
  .nrow = 8
)

endocrine_mm_select <- c("Cpe", "Pcsk2", "Tph1", "Neurog3", "Pdx1")

mmcr %>%
  filter(symbol %in% endocrine_mm_select) %>%
  pull(crob) %>%
  append("KY21.Chr12.180") %>% # CPE
  plot_genes_and_save("endocrine_select")


## Nuclear receptors ==========================================================
go_nuclear_receptors <- "GO:0004879"
nuclear_receptors_mm <- get_mm_genes_for_GO_and_offspring(go_nuclear_receptors)

list(
  get_cr_orthologs(nuclear_receptors_mm),
  get_cr_genes_for_GO_and_offspring(go_nuclear_receptors)
) %>%
  unlist() %>%
  unique() %>%
  plot_genes_and_save("nuclear-receptors")


## specific transcription factors =============================================
specific_tfs_mm <- mgenes %>%
  filter(
    symbol == "Nr1h3" | # LXR alpha
      symbol == "Nr1h2" | # LXR beta
      symbol == "Id3" |
      str_detect(symbol, "^Hox") |
      symbol == "Spdef" |
      str_detect(symbol, "Hnf") |
      str_detect(symbol, "Foxa") |
      str_detect(symbol, "Nkx") |
      str_detect(symbol, "Ppar") |
      str_detect(symbol, "Pdx1") |
      str_detect(symbol, "Ptf1") |
      str_detect(symbol, "Neuro") |
      str_detect(symbol, "Mnx") |
      str_detect(symbol, "Spi1|Spib|Spic$") |
      str_detect(symbol, "Cebp") |
      str_detect(symbol, "Gata") |
      str_detect(symbol, "Pou") |
      str_detect(symbol, "Rar[a-z]$") |
      str_detect(symbol, "Rxr[a-z]$")
  ) %>%
  arrange(symbol) %>%
  pull(symbol)

specific_tfs_mm %>%
  get_cr_orthologs() %>%
  unique() %>%
  plot_genes_and_save("select-transcription-factors")

## xenobiotic genes ===========================================================
go_xenobiotic <- "GO:0006805"
xenobiotic_mm <- get_mm_genes_for_GO_and_offspring(go_xenobiotic)

list(
  get_cr_orthologs(xenobiotic_mm),
  get_cr_genes_for_GO_and_offspring(go_xenobiotic)
) %>%
  unlist() %>%
  unique() %>%
  plot_genes_and_save("xenobiotic", .nrow = 27, .order_by_expression = TRUE)

## bile acid synthesis ========================================================
go_bile_synthesis <- "GO:0006699"
bile_synthesis_mm <- get_mm_genes_for_GO_and_offspring(go_bile_synthesis)

list(
  get_cr_orthologs(bile_synthesis_mm),
  get_cr_genes_for_GO_and_offspring(go_bile_synthesis)
) %>%
  unlist() %>%
  unique() %>%
  plot_genes_and_save(
    "bile-acid-synthesis",
    .nrow = 9,
    .order_by_expression = TRUE
  )

## peroxisome =================================================================
go_peroxisome <- "GO:0005777"
peroxisome_mm <- get_mm_genes_for_GO_and_offspring(go_peroxisome)

peroxisome_cr <- list(
  get_cr_orthologs(peroxisome_mm),
  get_cr_genes_for_GO_and_offspring(go_peroxisome)
) %>%
  unlist() %>%
  unique()

plot_genes_and_save(
  peroxisome_cr[1:104],
  "peroxisome_1",
  .order_by_expression = TRUE,
  .ncol = 5,
  .nrow = 21
)

plot_genes_and_save(
  peroxisome_cr[105:205],
  "peroxisome_2",
  .order_by_expression = TRUE,
  .ncol = 5,
  .nrow = 20
)

## carotenoid metabolism ======================================================
go_carotenoid_metabolism <- "GO:0016116"
carotenoid_metabolism_mm <- get_mm_genes_for_GO_and_offspring(
  go_carotenoid_metabolism
)

list(
  get_cr_orthologs(carotenoid_metabolism_mm),
  get_cr_genes_for_GO_and_offspring(go_carotenoid_metabolism)
) %>%
  unlist() %>%
  unique() %>%
  plot_genes_and_save("carotenoid-metabolism")

## cell cycle genes ===========================================================
mmhs <- read_tsv(
  "https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt"
) %>%
  janitor::clean_names()

get_mm_orths <- function(x) {
  dbid <- mmhs$db_class_key[which(mmhs$symbol == x)]
  mmhs %>%
    filter(db_class_key %in% dbid) %>%
    filter(ncbi_taxon_id == 10090) %>%
    pull(symbol)
}

s_genes <- Seurat::cc.genes$s.genes %>%
  map(get_mm_orths) %>%
  unlist() %>%
  unique() %>%
  get_cr_orthologs() %>%
  unique()

g2m_genes <- Seurat::cc.genes$g2m.genes %>%
  map(get_mm_orths) %>%
  unlist() %>%
  unique() %>%
  get_cr_orthologs() %>%
  unique()

plot_genes_and_save(s_genes, "cell-cycle_S", .order_by_expression = TRUE)
plot_genes_and_save(
  g2m_genes,
  "cell-cycle_G2M",
  .nrow = 9,
  .order_by_expression = TRUE
)

## lipid uptake and storage genes =============================================
mmcr %>%
  filter(symbol %in% c("Cd36", "Pnpla2", "Lipe", "Plin2", "Nr1h3")) %>%
  pull(crob) %>%
  unique() %>%
  plot_genes_and_save(
    .title = "lipid-uptake-and-storage",
    .order_by_expression = FALSE,
    .ncol = 3,
    .nrow = 2
  )

# EXPERIMENTAL BELOW ==========================================================
stop()

mmcr %>%
  filter(str_detect(symbol, regex("Top\\d"))) %>%
  pull(crob) %>%
  unique() %>%
  plot_genes_faceted("top")


mmcr %>%
  filter(str_detect(symbol, regex("Top1"))) %>%
  pull(crob) %>%
  unique() %>%
  plot_genes_faceted("Top1")


genes %>%
  filter(str_detect(symbol, regex("muc", TRUE))) %>%
  pull(gene_id) %>%
  unique() %>%
  plot_genes_faceted("Muc", .order_by_expression = TRUE)


# catalase and peroxosomal genes ==============================================
genes %>%
  filter(str_detect(display_name, regex("catalase", TRUE))) %>%
  pull(gene_id) %>%
  plot_genes_faceted(.title = "Catalase")

genes %>%
  filter(str_detect(display_name, regex("superoxide", TRUE))) %>%
  pull(gene_id) %>%
  plot_genes_faceted(.title = NULL)

orth %>%
  filter(str_detect(symbol, "Pex")) %>%
  pull(crob) %>%
  unique() %>%
  plot_genes_faceted(.title = "Pex genes")


orth %>%
  filter(str_detect(symbol, "Prdx")) %>%
  pull(crob) %>%
  unique() %>%
  plot_genes_faceted(.title = "Prdx genes")


orth %>%
  filter(str_detect(symbol, "Sod")) %>%
  pull(crob) %>%
  unique() %>%
  plot_genes_faceted(.title = "Crob genes")

stop()

mmcr %>%
  filter(str_detect(symbol, "Abcg1")) %>%
  pull(crob) %>%
  unique() %>%
  plot_genes_faceted(.title = "Abcg1")

orth %>%
  filter(str_detect(symbol, "Npc")) %>%
  pull(crob) %>%
  unique() %>%
  plot_genes_faceted(.title = "")

get_cr_orthologs(nuclear_receptors_mm) %>%
  unique() %>%
  plot_genes_faceted(.title = "nuclear receptors")

unique(unlist(nuclear_receptors_cr)) %>%
  plot_genes_faceted(.title = "nuclear receptors")

c("Mnx1", "Cyp26a1") %>%
  get_cr_orthologs() %>%
  unique() %>%
  plot_genes_faceted(.title = "Rar signaling")

c("Nr1h4") %>%
  get_cr_orthologs() %>%
  plot_genes_faceted(.title = "Nr1h4")

plot_genes_faceted(specific_tfs_cr, "tf")

unique(unlist(xenobiotic_cr))[1:40] %>%
  plot_genes_faceted(.title = "xenobiotic")
unique(unlist(xenobiotic_cr))[41:80] %>%
  plot_genes_faceted(.title = "xenobiotic")
unique(unlist(xenobiotic_cr))[81:120] %>%
  plot_genes_faceted(.title = "xenobiotic")


unique(unlist(endocrine_cr)) %>%
  plot_genes_faceted("endocrine")

unique(unlist(bile_synthesis_cr))[1:40] %>%
  plot_genes_faceted("bile")

unique(unlist(peroxisome_cr))[41:80] %>%
  plot_genes_faceted("peroxi")


unique(unlist(carotenoid_metabolism_cr))[1:40] %>%
  plot_genes_faceted("carotenoid")


plot_genes_faceted(stem_cr, "stem")
plot_genes_faceted(cellcycle_cr, "cc")


# Hox genes ===================================================================
orth %>%
  filter(str_detect(symbol, "Hox|Cdx|Pdx\\d|Gsx")) %>%
  pull(crob) %>%
  unique() %>%
  plot_genes_faceted("Hox and Parahox")

hox <- cand %>%
  filter(group %in% c("hox", "parahox")) %>%
  filter(!is.na(gene_id)) %>%
  pull(gene_id) %>%
  plot_genes_faceted("Hox and ParaHox genes") +
  labs(caption = paste0("ParaHox genes: Gsx, IPF1/Pdx1, Cdx"))

save_plot(hox, "hox-and-parahox", 5.5)

# enteroendocrine prohormone processing enzymes ===============================

c(
  orth %>%
    filter(
      str_detect(symbol, "Pdx\\d") |
        str_detect(symbol, "Pcsk2") |
        str_detect(symbol, "Cpe") |
        str_detect(symbol, "Ins\\d")
    ) %>%
    pull(crob),
  genes %>%
    filter(str_detect(symbol, regex("Ins", TRUE))) %>%
    pull(gene_id)
) %>%
  unique() %>%
  plot_genes_faceted("Enteroendocrin / Prohormone processing")

eec_genes <- c(
  "ncbi:445683", # ipf1/pdx1,
  "ncbi:100181292", # one2one ortholog of Mouse Pcsk2
  "ncbi:100184198", # carboxypeptidase E
  "ncbi:100170003", # insulin-like 1
  "ncbi:100170005", # insulin-like 2
  "ncbi:100170004" # insulin-like 3
)

eec_genes_plot <- plot_genes_faceted(
  eec_genes,
  "Enteroendocrine-related genes"
) +
  labs(
    caption = paste0(
      "ncbi:445683/Ipf1 = Pdx1",
      "\n",
      "ncbi:100181292/cipc2 = 1:1 ortholog of mouse Pcsk2",
      "\n",
      "    (Ortholog for Pcsk1 not found)",
      "\n",
      "    (Ortholog for Pcsk2 cofactor Scg5 not found)"
    )
  )

save_plot(eec_genes_plot, "eec-genes", 5)

# brush border membrane enzymes ===============================================
bbme_plot <- cand %>%
  filter(group == "brush border membrane enzyme") %>%
  pull(gene_id) %>%
  plot_genes_faceted("Brush border membrane enzymes") +
  labs(
    caption = paste0(
      c(
        "Lactase-like:",
        "    ncbi:104265662 / GBA3",
        "    ncbi:108950792 / LCTL",
        "    ncbi:104266271 / GBA3",
        "Maltase-glucoamylase:",
        "    ncbi:100175461 / SI",
        "Aminopeptidase A:",
        "    ncbi:100186666 / ENPEP",
        "Aminopeptidase N:",
        "    ncbi:100180384 / ANPEP",
        "    ncbi:100176712 / ANPEP",
        "    ncbi:100187078 / ANPEP"
      ),
      collapse = "\n"
    )
  )

save_plot(bbme_plot, "brush-border-enzymes", 5.8)

# membrane transporters =======================================================
memb_transporteres_plot <- cand %>%
  filter(group == "membrane transporter") %>%
  pull(gene_id) %>%
  plot_genes_faceted("Membrane transporters") +
  labs(
    caption = paste0(
      c(
        "SGLT1 (SLC5A1):",
        "    ncbi:100184349 / SLC5A9",
        "    ncbi:100177241 / SLC5A9",
        "GLUT5 (SLC2A5):",
        "    ncbi:100181984 / SLC2A5",
        "    ncbi:104265822 / SLC2A5",
        "PEPT1 (SLC15A1):",
        "    ncbi:100178449 / SLC15A1",
        "    ncbi:100184982 / SLC15A1",
        "    ncbi:100186649 / SLC15A1",
        "B0AT1 (SLC6A19):",
        "    ncbi:100186059 / SLC6A19",
        "    ncbi:100175124 / SLC6A19"
      ),
      collapse = "\n"
    )
  )

save_plot(memb_transporteres_plot, "membrane-transporters", 5.8)

# pinocytosis =================================================================
pinocytosis_plot <- cand %>%
  filter(group == "pinocytosis") %>%
  pull(gene_id) %>%
  plot_genes_faceted("Pinocytosis")

save_plot(pinocytosis_plot, "pinocytosis", 2.75, 6)

# phagocytosis ================================================================
phagocytosis_plot <- cand %>%
  filter(group == "phagocytosis") %>%
  pull(gene_id) %>%
  plot_genes_faceted("Phagocytosis")

save_plot(phagocytosis_plot, "phagocytosis", 2.75, 6)

# endocytosis =================================================================
endocytosis_plot <- cand %>%
  filter(group == "endocytosis") %>%
  pull(gene_id) %>%
  plot_genes_faceted("Endocytosis") +
  labs(
    caption = paste0(c(
      "Cathepsins = lysosomal proteases"
    ))
  )

save_plot(endocytosis_plot, "endocytosis", 5.75)

# digestive enzymes ===========================================================
digestive_enzyme_genes <- c(
  (cand %>%
    filter(group == "digestive enzyme") %>%
    pull(gene_id)),
  (genes %>%
    filter(str_detect(display_name, "PTF1A")) %>%
    pull(gene_id)),
  (genes %>%
    filter(str_detect(display_name, "chymotrypsinogen A")) %>%
    pull(gene_id))
)

digestive_enzyme_plot <- plot_genes_faceted(
  digestive_enzyme_genes,
  "Digestive enzymes and PTF1A"
)

save_plot(digestive_enzyme_plot, "digestive-enzymes", 2.5, 9)

# lumen protection ============================================================
mucin_genes <- c(
  (cand %>%
    filter(symbol_or_name == "MUC5AC") %>%
    pull(gene_id)),
  (genes %>%
    filter(str_detect(display_name, "mucin")) %>%
    filter(
      gene_id %in%
        paste0(
          "ncbi:", # genes with reasonable expression
          c(
            "100182274",
            "100187099",
            "113474321",
            "100177485",
            "100178531",
            "100186539",
            "104265530",
            "100176669",
            "113474448",
            "104266756",
            "100181221",
            "104266831",
            "104266522",
            "100185544",
            "108949449",
            "108950998"
          )
        )
    ) %>%
    pull(gene_id))
)

mucin_genes_plot <- plot_genes_faceted(
  mucin_genes,
  "Mucin genes",
  .order_by_expression = TRUE
) +
  labs(
    caption = paste0(
      c(
        "ncbi:104266451 = MUC5AC (Iguchi et al 2023)",
        "\n",
        "secreted: MUC2, MUC5, MUC6, MUC19, MUC7, MUC8, MUC9, MUC20",
        "membrane-bound: MUC1, MUC3, MUC4, MUC12, MUC13, MUC15, MUC16, MUC17, MUC21, MUC22"
      ),
      collapse = "\n"
    )
  )

save_plot(mucin_genes_plot, "mucin", 8)

# adipocyte signatures ========================================================
# Markers of adipocytes or genes likely associated with cells that contain lipid
# droplets.

## CD36 / fatty acid translocase
adipocyte_genes <- append(
  adipocyte_genes,
  orth %>%
    filter(str_detect(mm_symbol, "Cd36")) %>%
    pull(cr_gene_id) %>%
    unique()
)

# hormone-sensitive lipase / lipase E
adipocyte_genes <- append(
  adipocyte_genes,
  orth %>%
    filter(str_detect(mm_symbol, "Lipe")) %>%
    pull(cr_gene_id) %>%
    unique()
)

# ATGL / adipose tridlyceride lipase
adipocyte_genes <- append(
  adipocyte_genes,
  orth %>%
    filter(str_detect(mm_symbol, "Pnpla2")) %>% # ATGL
    pull(cr_gene_id) %>%
    unique()
)

# Perilipin
adipocyte_genes <- append(
  adipocyte_genes,
  orth %>%
    filter(str_detect(mm_symbol, "Plin")) %>%
    pull(cr_gene_id) %>%
    unique()
)

p_adipocyte <- plot_genes_faceted(adipocyte_genes, "Adipocyte-related genes") +
  labs(
    caption = paste0(
      c(
        "CD36 = Fatty acid translocase",
        "Adipocyte lipases:",
        "    LIPE = Hormone sensitive lipase / Lipase E",
        "    PNPLA2 = ATGL (Adipose triglyceride lipase)"
      ),
      collapse = "\n"
    )
  )

save_plot(p_adipocyte, "adipocyte-related", 3.25)


## FABPs from Orito et al 2015 ================================================
p_fabp <- xdb %>%
  filter(
    str_detect(ghost_gene_id, "KH.C6.206") | # intestinal
      str_detect(ghost_gene_id, "KH.C7.153") | # liver
      str_detect(ghost_gene_id, "KH.C12.472") | # brain / adipocyte
      str_detect(ghost_gene_id, "KH.L96.77") | # liver
      str_detect(ghost_gene_id, "KH.L119.12")
  ) %>%
  pull(gene_id) %>%
  unique() %>%
  plot_genes_faceted("FABP")

save_plot(p_fabp, "fabp", 2.75)

# A-P genes from Zwick et al ==================================================
orth %>%
  filter(str_detect(mm_symbol, "Meis2|Fabp1|Fabp6|Plb|Slc10a2")) %>%
  pull(cr_gene_id) %>%
  unique() %>%
  plot_genes_faceted("test")

# beta-galactosidase ==========================================================

genes %>%
  filter(str_detect(
    display_name,
    regex("beta-galactosidase", ignore_case = TRUE)
  )) %>%
  pull(gene_id) %>%
  plot_genes_faceted("beta-galactosidase") %>%
  save_plot("beta-galactosidase", 2.75, 3)

# end =========================================================================

# Arun Chavan
# Started: 2024-01-01

# background ==================================================================

# Here we generate plots for various categories of candidate genes. 

# setup =======================================================================
library(tidyverse)
library(here)
library(ggrepel)

basedir <- here("species/bflo/bulk/KL01")
genomedir <- here("species/bflo/genome/VNyyK")
resdir <- here(basedir, "results/05a_candidate-genes")
fs::dir_create(resdir)

source(here(basedir, "code/utilities_ggplot.R"))

# data ========================================================================
tpm <- read_tsv(here(basedir, "data/03_compiled/tpm_gut.tsv"))
samples <- read_tsv(here(basedir, "data/00_metadata/samples_gut.tsv"))

genes <- read_tsv(here(genomedir, "results/genes/gene-display-names.tsv"))
orth <- read_csv(here(genomedir, 
                      "results/emapper/processed/emapper_orthologs_mouse.csv"))

# for plotting ================================================================
tissue_order <- c("mdiv", "mgut", "hgut")

tpm_long <- tpm %>% 
  pivot_longer(cols = -gene_id, names_to = "sample_name", values_to = "tpm") %>% 
  left_join(samples) %>% 
  mutate(tissue = factor(tissue, tissue_order))

plot_genes_faceted <- function(.gene_ids, 
                               .title, 
                               .order_by_expression = FALSE) {
  
  # create plot data
  plot_data <- tpm_long %>% 
    filter(gene_id %in% .gene_ids) %>% 
    left_join(genes) %>% 
    mutate(display_name = paste0("[", gene_id, "] ", symbol, "\n", 
                                 str_trunc(description, 72)))
  
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
  
  plot_data <- mutate(plot_data, display_name = factor(display_name, gene_order))
  plot_data <- arrange(plot_data, display_name)
  
  # plot
  ggplot(plot_data, aes(x = tissue)) +
    geom_jitter(aes(y = tpm), width = 0.1, color = "grey") +
    stat_summary(aes(y = tpm, group = 1), fun = mean, geom = "line") +
    stat_summary(aes(y = tpm), fun = mean, geom = "point") +
    facet_wrap(vars(display_name), 
               scales = "free_y",
               nrow = 8, ncol = 6) +
    scale_x_discrete(guide = guide_axis(n.dodge = 1)) +
    theme_bw_2() +
    theme(
      strip.text = element_text(size = 6, hjust = 0),
      plot.caption = element_text(hjust = 0)
      ) +
    ggtitle(.title)
  
}

save_plot <- function(p, name, h = 11, w = 10) {
  ggsave(here(resdir, paste0(name, ".pdf")), 
         p, width = w, height = h, device = cairo_pdf)
}

# Hox genes ===================================================================
hox <- genes %>% 
  filter(str_detect(display_name, regex("gsx|pdx1|cdx|hox", TRUE))) %>% 
  pull(gene_id) %>% 
  plot_genes_faceted("Hox and ParaHox genes")

hox <- cand %>% 
  filter(group %in% c("hox", "parahox")) %>% 
  filter(!is.na(gene_id)) %>% 
  pull(gene_id) %>% 
  plot_genes_faceted("Hox and ParaHox genes") +
  labs(caption = paste0("ParaHox genes: Gsx, IPF1/Pdx1, Cdx"))

save_plot(hox, "hox-and-parahox", 8)

# enteroendocrine prohormone processing enzymes ===============================

eec_genes <- c(
  "ncbi:118419641", # ipf1/pdx1,
  "ncbi:118409753", # one2one ortholog of Mouse Pcsk1
  "ncbi:118418937", # one2one ortholog of Mouse Pcsk2
  "ncbi:118423235", # one2one ortholog of mouse Scg5 (Pcsk2 cofactor) 
  "ncbi:118411334", # carboxypeptidase E-like 
  "ncbi:118414102", # carboxypeptidase E-like  
  "ncbi:118407887"  # insulin-like peptide
)

eec_genes_plot <- plot_genes_faceted(eec_genes, "Enteroendocrine-related genes") +
  labs(caption = paste0(
    "ncbi:118419641 = Pdx1", "\n",
    "ncbi:118409753 = 1:1 ortholog of mouse Pcsk1", "\n",
    "ncbi:118418937 = 1:1 ortholog of mouse Pcsk2", "\n"
  )) 

save_plot(eec_genes_plot, "eec-genes", 2.8)

# fluorescent proteins ========================================================
fluor_plot <- genes %>% 
  filter(str_detect(description, regex("fluorescent", TRUE))) %>% 
  pull(gene_id) %>% 
  plot_genes_faceted(., "Fluorescent proteins")

save_plot(fluor_plot, "fluorescent-proteins", w = 10.5, h = 6)

# exocrine pancreatic genes ===================================================

###############################################################################
# Below this is the template from crob candidate genes. 
# TODO: Need to make a list of candidate genes in bflo and make plots
###############################################################################

# brush border membrane enzymes ===============================================
bbme_plot <- cand %>% 
  filter(group == "brush border membrane enzyme") %>% 
  pull(gene_id) %>% 
  plot_genes_faceted("Brush border membrane enzymes") +
  labs(caption = paste0(c(
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
  ), collapse = "\n"))

save_plot(bbme_plot, "brush-border-enzymes", 5.8)

# membrane transporters =======================================================
memb_transporteres_plot <- cand %>% 
  filter(group == "membrane transporter") %>% 
  pull(gene_id) %>% 
  plot_genes_faceted("Membrane transporters") +
  labs(caption = paste0(c(
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
  ), collapse = "\n"))
  
save_plot(memb_transporteres_plot, "membrane-transporters", 5.8)

# pinocytosis =================================================================
pinocytosis_plot <-  cand %>% 
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
  labs(caption = paste0(c(
    "Cathepsins = lysosomal proteases"
  )))

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

digestive_enzyme_plot <- plot_genes_faceted(digestive_enzyme_genes, 
                                            "Digestive enzymes and PTF1A")

save_plot(digestive_enzyme_plot, "digestive-enzymes", 2.5, 9)

# lumen protection ============================================================
mucin_genes <- c(
  (cand %>% 
     filter(symbol_or_name == "MUC5AC") %>% 
     pull(gene_id)),
  (genes %>% 
     filter(str_detect(display_name, "mucin")) %>% 
     filter(gene_id %in% paste0("ncbi:", # genes with reasonable expression
                                c("100182274", 
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
                                  "108950998"))) %>% 
     pull(gene_id))
)

mucin_genes_plot <- plot_genes_faceted(mucin_genes, 
                                       "Mucin genes", 
                                       .order_by_expression = TRUE) +
  labs(caption = paste0(c(
    "ncbi:104266451 = MUC5AC (Iguchi et al 2023)", "\n",
    "secreted: MUC2, MUC5, MUC6, MUC19, MUC7, MUC8, MUC9, MUC20",
    "membrane-bound: MUC1, MUC3, MUC4, MUC12, MUC13, MUC15, MUC16, MUC17, MUC21, MUC22"
  ), collapse = "\n"))

save_plot(mucin_genes_plot, "mucin", 8)

# FoxQ2c (Giacomo Gattoni) ====================================================
plot_genes_faceted("ncbi:118416625", "FoxQ2c") %>% 
  save_plot("FoxQ2c", w = 4, h = 3)

# end =========================================================================


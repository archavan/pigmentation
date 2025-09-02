# Arun Chavan
# 2024-11-13

# setup =======================================================================
library(Seurat)
library(here)
library(tidyverse)

scdir <- here("species/crob/10x/run-02/KY21")
basedir <- here("species/crob/bulk/ARCA123/KY21")
genomedir <- here("species/crob/genome/KY21")

# non single-cell data ========================================================
tpm <- read_csv(here(basedir, "data/03_compiled/tpm.csv"))
counts <- read_csv(here(basedir, "data/03_compiled/counts.csv"))

genes <- read_tsv(here(genomedir, "results/genes/gene-display-names.tsv"))

# sc data =====================================================================

merge_samples <- function(.obj_list) {
  merge(x = .obj_list[[1]],
        y = .obj_list[2:length(.obj_list)],
        add.cell.ids = names(.obj_list))
}

cr_samples_with_sperm <- tribble(
  ~animal, ~tissue, ~fileprefix,
  "cr41",  "int",   "crob41_int",
  "cr41",  "stm",   "crob41_stm",
  "cr42",  "int",   "crob42_int",
  "cr42",  "stm",   "crob42_stm"
)

cr_samples_with_sperm <- cr_samples_with_sperm %>% 
  mutate(seur = map(fileprefix, 
                    ~ read_rds(here(scdir, 
                                    "data/processed/01_qc-filtered/",
                                    paste0(.x, "_seurat.rds")))))

cr_samples_with_sperm <- cr_samples_with_sperm %>% 
  mutate(barcode_prefix = paste0(animal, "_", tissue))


sc_sperm <- merge_samples(
  set_names(cr_samples_with_sperm$seur, cr_samples_with_sperm$barcode_prefix)
)

cr_samples_no_sperm <- tribble(
  ~animal, ~fileprefix,
  "cr41",  "crob41",
  "cr42",  "crob42"
) %>% 
  mutate(seur = map(
    fileprefix,
    ~ read_rds(here(scdir,
                    "data/processed/02_cleaned/",
                    paste0(.x, "_clean_v1_no-sperm.rds")))
    ))


sc_clean <- merge_samples(
  set_names(cr_samples_no_sperm$seur, cr_samples_no_sperm$animal)
)

testis_barcodes <- colnames(sc_sperm)[which(
  !colnames(sc_sperm) %in% colnames(sc_clean)
)]

sc_sperm$tissue2 <- ifelse(rownames(sc_sperm[[]]) %in% testis_barcodes,
                           "testis",
                           "gut")

# find testis specific genes ==================================================
sc_process_data <- function(.obj) {
  
  npcs <- 30
  dims <- 1:npcs
  kparam <- 20
  
  .obj <- .obj %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
    ScaleData() %>%
    RunPCA(npcs = npcs)
  .obj <- RunUMAP(.obj, reduction = "pca", dims = dims, n.neighbors = kparam)
  return(.obj)  
}

sc_sperm <- sc_process_data(sc_sperm)
sc_sperm <- SetIdent(sc_sperm, value = "tissue2")

DimPlot(sc_sperm, label = TRUE, repel = TRUE, shuffle = TRUE)

aggregate_then_normalize <- function(x, clust_prefix) {
  aggcts <- GetAssayData(x, slot = "counts") %*% model.matrix(~ 0 + Idents(x))
  colnames(aggcts) <- gsub("Idents\\(x\\)", clust_prefix, colnames(aggcts))
  aggcts <- as.matrix(aggcts)
  normdata <- apply(aggcts, 2, function(y) (y*10000)/sum(y))
  normdata <- log1p(normdata)
  
  normdata %>% 
    as.data.frame() %>% 
    rownames_to_column("gene_id") %>% 
    as_tibble()
}

avg <- aggregate_then_normalize(sc_sperm, "") %>% 
  left_join(select(genes, gene_id, symbol, description))

# after trying various criteria, I am going with the following:
# 1. top 500 genes that have more than 2 fold expression in testis than in gut
# 2. genes with sperm/testis in their description

higher_in_testis <- avg %>% 
  filter(testis > 2*gut) %>% 
  slice_max(order_by = testis, n = 500)

by_description <- avg %>% 
  filter(str_detect(description, regex("sperm|testis", TRUE))) 

testis_genes <- bind_rows(higher_in_testis, by_description) %>% 
  distinct()
  
# write testis genes
fs::dir_create(here(basedir, "results/04b_testis-genes"))
write_csv(testis_genes, here(basedir, "results/04b_testis-genes/testis-enriched-genes.csv"))

# filter tpm and counts =======================================================
tpm %>% 
  filter(!gene_id %in% testis_genes$gene_id) %>% 
  write_csv(here(basedir, "data/03_compiled/tpm_no-testis.csv"))

counts %>% 
  filter(!gene_id %in% testis_genes$gene_id) %>% 
  write_csv(here(basedir, "data/03_compiled/counts_no-testis.csv"))

# plot testis specific genes ==================================================
z_testis <- tpm %>% 
  pivot_longer(cols = -gene_id, names_to = "sample_name", values_to = "tpm") %>% 
  separate(sample_name, into = c("tissue", "rep_id"), sep = "_") %>% 
  filter(gene_id %in% testis_genes$gene_id) %>% 
  group_by(tissue, gene_id) %>% 
  summarize(tpm = mean(tpm)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "tissue", values_from = "tpm") %>% 
  column_to_rownames("gene_id") %>% 
  t() %>% 
  scale() %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("gene_id") %>% 
  pivot_longer(cols = -gene_id, values_to = "z_score", names_to = "tissue") %>% 
  mutate(tissue = factor(tissue, levels = c("stm.1.s", 
                                            "int.1.w", "int.2.y",
                                            "int.3.w", "int.4.y")))

testis_genes_expression_plot <- ggplot(z_testis, aes(tissue, z_score)) +
  geom_line(aes(group = gene_id), alpha = 0.2, color = "grey40") +
  theme_bw()

ggsave(here(basedir, "results/04b_testis-genes/testis-genes-expression.png"), 
       plot = testis_genes_expression_plot, width = 6, height = 4, dpi = 300)

# end ==========================================================================


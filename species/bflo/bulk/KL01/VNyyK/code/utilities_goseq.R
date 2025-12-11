# setup =======================================================================
library(goseq)
library(tidyverse)

# genomic data ================================================================
glen <- read_tsv("species/bflo/genome/VNyyK/results/genes/gene-length-median.tsv")
go <- read_tsv("species/bflo/genome/VNyyK/results/emapper/processed/emapper_go.tsv.gz") 

# function for running go enrichment ==========================================
#' run geneset enrichment with goseq
#'
#' @param .target character vector of target gene ids
#' @param .bg character vector of background gene ids
#' @param .db dataframe specifying relationship between gene_id, cat_id, cat_name
#' @param .cat_id_col name of column with category ids in .db
#' @param .cat_name_col name of column with category names in .db
#' @param .gene_id_col name of column with gene ids in .db
#' @param .gene_lengths named numeric vector of gene lengths (name = gene_id)
#'
#' @return tibble with enrichment results

run_goseq <- function(.target, 
                      .bg,
                      .gene_lengths = deframe(glen)) {
  
  # {goseq} expects a vector of target genes as a named integer vector where,
  # the length of the vector is the length of background gene vector, names of
  # the vector are the names of the background genes, and the values are 1 or 0
  # for whether a gene is in the target set or not, respectively.
  
  .target <- as.integer(.bg %in% .target) %>% 
    set_names(.bg)
  
  # only retain genes that are in the backgroud set
  filtered_gene2cat <- go[which(go[["gene_id"]] %in% .bg), ]
  filtered_gene2cat <- as.data.frame(filtered_gene2cat) # goseq can't handle tibbles
  
  # category id and category name mapping
  cat2name <- go %>% 
    select(all_of(c("go_id", "go_term"))) %>% 
    distinct()
  
  # start running goseq
  pwf <- nullp(DEgenes = .target, 
               bias.data = .gene_lengths[.bg], 
               plot.fit = FALSE)
  
  res <- goseq(pwf = pwf, 
               gene2cat = filtered_gene2cat, 
               method = "Hypergeometric") # ignore gene length 
  
  # augment goseq output
  res_aug <- res %>% 
    as_tibble() %>% 
    select(-any_of(c("under_represented_pvalue", "term"))) %>% 
    dplyr::rename(go_id = category) %>% 
    mutate(enrichment = (numDEInCat / sum(.target)) / 
             (numInCat / length(unique(filtered_gene2cat[["gene_id"]])))) %>% 
    mutate(p_adj = p.adjust(p = over_represented_pvalue, method = "BH")) %>% 
    left_join(cat2name, by = "go_id") # add category name
  
  return(res_aug)
  
}

# end =========================================================================



# Helper functions for clustering of data in Seurat objects. 

# Dependencies:
#     - Seurat
#     - tidyverse
#     - glue

###############################################################################
# rename clusters by size ====
###############################################################################

#' Rename clusters by size
#' 
#' Rename clusters by size and set new idents as default.
#'
#' @param .idents.old idents to rename
#' @param .idents.new name for new idents
#' @param .obj Seurat object
#' @param .remove.old whether to remove old ident column
#'
#' @return Seurat object

rename_clusters_by_size <- function(.obj, 
                                    .idents.old, 
                                    .idents.new,
                                    .remove.old = TRUE) {
  
  clusters_by_freq <- names(sort(table(.obj[[.idents.old, drop = TRUE]]), 
                                 decreasing = TRUE))
  .obj[[.idents.new]] <- factor(.obj[[.idents.old, drop = TRUE]], 
                                levels = clusters_by_freq) %>% 
    as.numeric() %>% 
    as.factor()
  
  .obj <- SetIdent(.obj, value = .idents.new)
  
  if (.remove.old) .obj[[.idents.old]] <- NULL
  
  return(.obj)
  
}

###############################################################################
# calculate pairwise correlations ====
###############################################################################

#' Calculate pairwise correlations of clusters
#' 
#' Calculate pairwise pearson correlation of clusters using cluster-wise average
#'  expression (log1pCPT; `data` slot from the default assay) of highly variable
#'  genes. 
#'
#' @param .obj Seurat object
#'
#' @return data.frame of pairwise pearson correlation coefficients

calculate_pairwise_correlations <- function(.obj) {
  
  avg_exp <- AverageExpression(.obj, 
                               assays = DefaultAssay(.obj), 
                               slot = "data",
                               features = VariableFeatures(.obj), 
                               return.seurat = FALSE, 
                               verbose = FALSE) %>% 
    .[[DefaultAssay(.obj)]] %>% 
    log1p() 
  
  corrs <- cor(avg_exp, avg_exp)
  corrs[lower.tri(corrs, diag = TRUE)] <- NA
  corrs <- as.data.frame(as.table(corrs)) # melt
  corrs <- na.omit(corrs)
  corrs <- dplyr::rename(corrs, corr = Freq, clust1 = Var1, clust2 = Var2)
  
  return(corrs)
  
}

###############################################################################
# get number of DE genes ====
###############################################################################

#' Get number of DE genes
#'
#' @param .obj Seurat object
#' @param .clust1 ident 1
#' @param .clust2 ident 2
#' @param .thresh_log2FC log2 fold change threshold
#' @param .thresh_padj p.adj threshold
#'
#' @return integer
#'
#' @examples get_number_of_de_genes(pbmc_small, 0, 1, 1, 0.05)
get_number_of_de_genes <- function(.obj, 
                                   .clust1, 
                                   .clust2, 
                                   .thresh_log2FC, 
                                   .thresh_padj) {
  
  de_genes <- FindMarkers(.obj,
                          ident.1 = .clust1,
                          ident.2 = .clust2,
                          assay = "RNA", 
                          features = VariableFeatures(.obj), 
                          base = 2,
                          verbose = FALSE) %>% 
    filter(p_val_adj < .thresh_padj, 
           abs(avg_log2FC) > .thresh_log2FC)
  
  n_de <- nrow(de_genes)
  
  return(n_de)
  
}

###############################################################################
# combine two clusters ====
###############################################################################

#' Combine two clusters
#' 
#' Combine two clusters into one. This creates a new metadata column with the 
#'  new cluster names (= oldname1.oldname2) and sets the new clusters as the 
#'  default Ident. 
#'
#' @param .obj Seurat object
#' @param .clust1 first cluster to combine (quoted or unquoted)
#' @param .clust2 second cluster to combine (quoted or unquoted)
#' @param .new.ident name for the new metadata column (character string)
#'
#' @return Seurat object
#'
#' @examples combine_two_clusters(pbmc_small, 2, 1, "temp")

combine_two_clusters <- function(.obj, .clust1, .clust2, .new.ident) {
  
  .clust1 <- as.character(.clust1)
  .clust2 <- as.character(.clust2)
  
  new_cluster_name <- paste0(.clust1, ".", .clust2)
  
  clusters_merged <- as.character(Idents(.obj))
  clusters_merged[clusters_merged == .clust1] <- new_cluster_name
  clusters_merged[clusters_merged == .clust2] <- new_cluster_name
  
  .obj[[.new.ident]] <- factor(clusters_merged, 
                               levels = sort(unique(clusters_merged)))
  
  .obj <- SetIdent(.obj, value = .new.ident)
  
  return(.obj)
  
}

###############################################################################
# merge closest pair of clusters ====
###############################################################################

#' Merge closes pair of clusters
#' 
#' Merge the two most similar clusters and stop. 
#'
#' @param .obj Seurat object
#' @param .merged.idents.col name for the merged idents column
#' @param .thresh_corr threshold corr. coeff. below which not to merge clusters
#' @param .thresh_n_de threshold n of DE genes below which not to merge clusters
#' @param .thresh_log2FC threshold log2FC to consider genes DE
#'
#' @return Seurat object

merge_closest_pair_of_clusters <- function(.obj, 
                                           .merged.idents.col,
                                           .thresh_corr,
                                           .thresh_n_de,
                                           .thresh_log2FC) {
  
  .obj[[.merged.idents.col]] <- Idents(.obj)
  .obj <- SetIdent(.obj, value = .merged.idents.col)
  
  corrs <- calculate_pairwise_correlations(.obj) %>% 
    filter(corr > .thresh_corr) %>% 
    arrange(desc(corr))
  
  if (nrow(corrs) == 0) {
    .obj <- .obj # if corr threshold not met, don't merge anything
  } else {
    for (i in seq_len(nrow(corrs))) {
      
      clust1 <- as.character(corrs$clust1[i])
      clust2 <- as.character(corrs$clust2[i])
      
      n_de <- get_number_of_de_genes(.obj = .obj,
                                     .clust1 = clust1,
                                     .clust2 = clust2,
                                     .thresh_log2FC = .thresh_log2FC,
                                     .thresh_padj = 0.05)
      
      if (n_de < .thresh_n_de) {
        .obj <- combine_two_clusters(.obj, clust1, clust2, .merged.idents.col)
        message(glue("  Merged clusters {clust1} and {clust2} ", 
                     "(corr = {round(corrs$corr[i], 2)}, DE genes = {n_de})"))
        break # right after merging the first pair
      } else next
    }
  }
  
  return(.obj)
  
}

###############################################################################
# merge similar clusters ====
###############################################################################

#' Merge similar clusters
#' 
#' Merge similar clusters that meet certain thresholds: fewer than `.thresh_n_de` 
#'  number of differentially expressed genes and higher correlation coefficient
#'  than `.thresh_corr`. 
#'
#' @param .obj Seurat object
#' @param .idents.new name for the merged clusters column, passed to `rename_clusters_by_size()`
#' @param .thresh_corr passed to `merge_closest_pair_of_clusters()`
#' @param .thresh_n_de passed to `merge_closest_pair_of_clusters()`
#' @param .thresh_log2FC passed to `merge_closest_pair_of_clusters()`
#'
#' @return Seurat object. Merged clusters are in `.idents.new` column which has 
#'   been set as the default ident. 
#' @export
#'
#' @examples
merge_similar_clusters <- function(.obj, 
                                   .idents.new,
                                   .thresh_corr = 0.8,
                                   .thresh_n_de = 20,
                                   .thresh_log2FC = 1) {
  
  # set up
  time_start <- Sys.time()
  counter <- 1
  
  # merge clusters
  repeat {
    message(glue("iteration {counter}"))
    
    n_clusters_old <- length(unique(Idents(.obj)))
    
    .obj <- merge_closest_pair_of_clusters(.obj = .obj,
                                           .merged.idents.col = "merged_temp",
                                           .thresh_corr = .thresh_corr, 
                                           .thresh_n_de = .thresh_n_de, 
                                           .thresh_log2FC = .thresh_log2FC)
    
    # check if any clusters got merged
    n_clusters_new <- length(unique(Idents(.obj)))
    no_clusters_merged <- isTRUE(all.equal(n_clusters_old, n_clusters_new))
    
    if (no_clusters_merged) {
      message("  No clusters merged")
      break
    } else counter <- counter + 1
  }
  
  message("Finished merging clusters.")
  
  # rename newly merged clusters by size
  .obj <- rename_clusters_by_size(.obj = .obj, 
                                  .idents.old = "merged_temp", 
                                  .idents.new = .idents.new)
  message(glue("Renamed {n_clusters} clusters by size.",
               n_clusters = length(unique(.obj[[.idents.new, drop = TRUE]]))))
  
  # wrap up
  time_end <- Sys.time()
  time_elapsed <- time_end - time_start
  message(glue("Elapsed time: {round(time_elapsed[[1]], 2)} ", 
               "{attr(time_elapsed, \"units\")}."))
  
  return(.obj)
  
}

# end =========================================================================



#' Run Fisher's exact test to test significance of overlap between two gene sets.
#'
#' @param geneset_1 vector of gene identifiers representing the first gene set.
#' @param geneset_2 vector of gene identifiers representing the second gene set.
#' @param universe vector of gene identifiers representing the background
#'   universe of genes.
#' @param alternative A character string specifying the alternative hypothesis,
#'   must be one of "two.sided" (default), "greater" or "less".
#'
#' @returns A list containing the results of Fisher's exact test, including
#'   the p-value and odds ratio.
#'
#' @export
#' @examples
run_fisher_overlap_test <- function(
  geneset_1,
  geneset_2,
  universe,
  alternative = "greater"
) {
  geneset_1 <- unique(geneset_1)
  geneset_2 <- unique(geneset_2)
  universe <- unique(universe)

  stopifnot(all(c(geneset_1, geneset_2) %in% universe))

  overlap <- length(intersect(geneset_1, geneset_2))
  only_geneset_1 <- length(setdiff(geneset_1, geneset_2))
  only_geneset_2 <- length(setdiff(geneset_2, geneset_1))
  neither <- length(setdiff(universe, union(geneset_1, geneset_2)))

  contingency_table <- matrix(
    c(overlap, only_geneset_1, only_geneset_2, neither),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(geneset_1 = c("yes", "no"), geneset_2 = c("yes", "no"))
  )

  fisher.test(contingency_table, alternative = alternative)
}

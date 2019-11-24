#' Greedy mapping
#'
#' Map trait name to gene if the trait name shared at least one trait_code with the gene, where
#' trait name to trait_code is specified in map_table and
#' gene to trait_code is specified in silver_table
#'
#' @param map_table a data.frame containing 'trait' as column along with 'phecode' and/or 'HPO' as columns.
#' This table serves as the map between trait name in 'trait' column and trait codes ('phecode', 'HPO', etc).
#' So that it is used to map trait names (in score_table) to trait code (in silver_standard).
#' @param silver_table a data.frame where each row pairs trait in phecode, HPO to gene.
#' @param trait_code trait code to use for mapping: 'phecode' or 'HPO'.
#'
#' @return a data.frame pairing trait name to gene
#'
#'
#' @export
#' @import dplyr

greedy_map = function(map_table, silver_table, trait_code) {
  silver_table = silver_table[!is.na(silver_table[, trait_code]), ]
  map_table = map_table[!is.na(map_table[, trait_code]), ]
  trait_gene_silver = map_table %>% inner_join(silver_table, by = trait_code) %>% select(.data$trait, .data$gene) %>% distinct()
}

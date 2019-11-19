#' Plotting PR/ROC curves
#'
#' Given a silver standard as 'truth', plot PR/ROC curves for the trait/gene pairs and scores listed in score_table, e.g.
#'
#' > score_table = data.frame(trait = ..., gene = ..., score1 = ..., score2 = ..., ...)
#'
#' Silver standard supports HPO, phecode.
#' User should provide the corresponding HPO or phecode for their GWAS traits in map_table, e.g.
#'
#' > map_table1 = data.frame(trait = ..., phecode = ...)
#'
#' > map_table2 = data.frame(trait = ..., HPO = ...)
#'
#'
#' @param score_table a data.frame containing 'trait' and 'gene' as columns along with
#' multiple columns for scores calculated by different method. Note that we interpret HIGHER score as
#' more likely to be causal. If you have the opposite, please flip it.
#' @param map_table a data.frame containing 'trait' as column along with 'phecode' and/or 'HPO' as columns.
#' This table serves as the map between trait name in 'trait' column and trait codes ('phecode', 'HPO', etc).
#' So that it is used to map trait names (in score_table) to trait code (in silver_standard).
#' @param silver_standard shared along with the package. It is a list with meta information and
#' a table (a data.frame) where each row pairs trait in phecode, HPO to gene.
#' New silver_standard data set is awaiting for contribution from the community.
#' @param score_cols a vector of column names in score_table to be used as score columns. If it is NULL, all
#' but 'trait' and 'gene' columns will be used.
#' @param trait_code trait code to use for mapping: 'phecode' or 'HPO'.
#' @param mapping_method the name of the mapping function used to match trait name with gene
#' on the basis of map_table and silver_standard$table
#' New mapping_method is awaiting for contribution from the community.
#'
#'
#' @return list containing PR and ROC plots which are ggplot2 objects
#'
#' @examples
#' silver_standard_proto(
#'   score_table = data.frame(
#'     trait = c(rep('t1', 20), rep('t2', 30), rep('t3', 50)),
#'     gene = paste0('g', sample(1:100, size = 100, replace = TRUE)),
#'     score1 = runif(100),
#'     score2 = runif(100)
#'   ),
#'   map_table = data.frame(
#'     trait = c('t1', 't2', 't3'),
#'     phecode = c('1.1', '23', '12')
#'   ),
#'   silver_standard = list(
#'     table = data.frame(
#'       phecode = c(rep('1.1', 10), rep('23', 10), rep('12', 10), rep('2', 10)),
#'       gene = paste0('g', sample(1:100, size = 40, replace = TRUE))
#'     ),
#'     script_info = 'toy_example'
#'   )
#' )
#'
#' @export
#' @import dplyr
#' @import ggplot2
silver_standard_proto = function(score_table, map_table, silver_standard, score_cols = NULL, trait_code = 'phecode', mapping_method = 'greedy_map') {
  message('Run with silver standard from: ', silver_standard$script_info)
  message('Map trait by: ', trait_code)
  message('Mapper chosen: ', mapping_method)
  mapper = match.fun(mapping_method)
  if(! trait_code %in% colnames(map_table)) {
    message('ERROR: The input map_table does not have the trait_code you select')
  }
  if(is.null(score_cols)) {
    message('Extracting all columns of score_table other than "trait" and "gene" as scores')
    score_cols = colnames(score_table)
    score_cols = score_cols[!score_cols %in% c('gene', 'trait')]
  }
  message(length(score_cols), ' score columns are used')

  # step 1
  # join trait - code (map_table) and code - gene (silver_standard$table) tables
  # to get trait - gene table for silver standard pairs
  trait_gene_silver = mapper(map_table, silver_standard$table, trait_code)

  # step 2
  # subset score_table to keep traits successfully matched to silver standard
  score_table_subset = score_table %>% filter(.data$trait %in% trait_gene_silver$trait)

  # step 3
  # create candidate trait/gene pair pool
  # and silver standard trait/gene pair pool
  trait_gene_candidate_pool = pair_up_trait_and_gene(score_table_subset$trait, score_table_subset$gene)
  trait_gene_silver_pool = pair_up_trait_and_gene(trait_gene_silver$trait, trait_gene_silver$gene)

  # step 4
  # subset silver standard pool to keep trait/gene pairs occurring candidate pool
  trait_gene_silver_pool = trait_gene_silver_pool[trait_gene_silver_pool %in% trait_gene_candidate_pool]

  pr_curves = list()
  for(i in score_cols) {
    pr_curve = gen_fdr_power_curve(true_genes = trait_gene_silver_pool, gene = trait_gene_candidate_pool, score = score_table_subset[, i])
    pr_curves[[length(pr_curves) + 1]] = pr_curve[-nrow(pr_curve), ] %>% mutate(score = i)

  }
  pr_curves = do.call(rbind, pr_curves)
  pr = pr_curves %>% ggplot() + geom_path(aes_string(x = 'power', y = 'precision', color = 'score'))


  roc_curves = list()
  for(i in score_cols) {
    roc_curve = gen_roc_curve(true_genes = trait_gene_silver_pool, gene = trait_gene_candidate_pool, score = score_table_subset[, i])
    roc_curves[[length(roc_curves) + 1]] = roc_curve %>% mutate(score = i)

  }
  roc_curves = do.call(rbind, roc_curves)
  auc = compute_auc(roc_curves)
  roc = roc_curves %>% ggplot() + geom_path(aes_string(x = 'fpr', y = 'tpr', color = 'score')) + coord_equal()
  return(list(pr = pr, roc = roc, roc_auc = auc$roc_auc))
}

pair_up_trait_and_gene = function(trait, gene) {
  paste(trait, gene, sep = '-x-')
}

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
#' @param gwas_loci a data.frame with columns: chromosome (e.g. chr1), start, end,
#' trait (with trait name matching the one in score_table) indicating which GWAS loci you want to focus on.
#' If it is specified as non-NULL, you should also specify gene_annotation as the genomic position for genes.
#' By specifying gwas_loci and gene_annotation, it triggers a more stringent
#' (recommended due to the sparse nature of silver standard) pre-processing step
#' to limit analysis on a subset of trait-gene pairs in score_table, by the following 2 steps:
#' 1. Select GWAS loci overlapping with silver standard gene
#' 2. Select genes that overlap with selected GWAS loci as candidate genes.
#' @param gene_annotation a data.frame with columns:
#' chromosome (e.g. chr1), start, end, gene (Ensembl ID (no dot), e.g. ENSG00000084754)
#' @param score_cols a vector of column names in score_table to be used as score columns. If it is NULL, all
#' but 'trait' and 'gene' columns will be used.
#' @param trait_codes trait codes to use for mapping: 'phecode' and/or 'HPO'.
#' @param mapping_method the name of the mapping function used to match trait name with gene
#' on the basis of map_table and silver_standard$table
#' New mapping_method is awaiting for contribution from the community.
#'
#'
#' @return list containing PR and ROC plots which are ggplot2 objects
#'
#' @examples
#' score_table = data.frame(
#'   trait = c(rep('t1', 20), rep('t2', 30), rep('t3', 50)),
#'   gene = paste0('g', 1:100),
#'   score1 = runif(100),
#'   score2 = runif(100),
#'   stringsAsFactors = FALSE
#' )
#' map_table = data.frame(
#'   trait = c('t1', 't2', 't3'),
#'   phecode = c('1.1', '23', '12'),
#'   stringsAsFactors = FALSE
#' )
#' silver_standard = list(
#'   table = data.frame(
#'     phecode = c(rep('1.1', 10), rep('23', 10), rep('12', 10), rep('2', 10)),
#'     gene = paste0('g', 1:40),
#'     stringsAsFactors = FALSE
#'   ),
#'   script_info = 'toy_example'
#' )
#' gwas_loci = data.frame(
#'   chromosome = paste0('chr', sample(1:1, size = 10, replace = TRUE)),
#'   start = (1:10) * 1e3 + 300, end = (1:10) * 1e3 + 500,
#'   trait = sample(c('t1', 't2'), size = 10, replace = TRUE)
#' )
#' gene_annot = data.frame(
#'   chromosome = paste0('chr', sample(1:1, size = 100, replace = TRUE)),
#'   gene_id = paste0('g', 1:100),
#'   start = 1:10 * 1e3 + 350,
#'   end = 1:10 * 1e3 + 450,
#'   gene_type = c(rep('protein_coding', 60), rep('psuedogene', 40))
#' )
#'
#' @export
#' @import dplyr
#' @import ggplot2
silver_standard_proto = function(score_table, map_table, silver_standard, gwas_loci = NULL, gene_annotation = NULL, score_cols = NULL, trait_codes = 'phecode', mapping_method = 'greedy_map') {
  message('Run with silver standard from: ', silver_standard$script_info)
  message('Map trait by: ', paste(trait_codes, collapse = ' ,'))
  message('Mapper chosen: ', mapping_method)
  mapper = match.fun(mapping_method)
  if(sum(trait_codes %in% colnames(map_table)) != length(trait_codes)) {
    message('ERROR: The input map_table does not have the trait_codes you select')
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
  trait_gene_silver_list = list()
  for(trait_code in trait_codes) {
    trait_gene_silver_list[[length(trait_gene_silver_list) + 1]] = mapper(map_table, silver_standard$table, trait_code) %>% mutate(trait_gene = pair_up_trait_and_gene(.data$trait, .data$gene))
  }
  trait_gene_silver = do.call(rbind, trait_gene_silver_list) %>% distinct()

  # step 2
  # if no pre-processing
  # subset score_table to keep traits successfully matched to silver standard
  # if pre-processing
  # perform pre-processing to limit to GWAS loci having silver standard genes and candidate genes overlapping selected GWAS loci
  message('# trait-gene pairs in score table before step 2: ', nrow(score_table))
  if(is.null(gwas_loci)) {
    message('gwas_loci is set to NULL, skip pre-processing step')
    score_table_subset = score_table %>% filter(.data$trait %in% trait_gene_silver$trait)
  } else {
    message('Start pre-processing step: limit to (1) GWAS loci with silver standard genes and (2) genes overlapping in selected GWAS loci')
    if(is.null(gene_annotation)) {
      message('gene_annotation is NULL, please specify the gene_annotation is you want to perform pre-processing')
      return(NULL)
    }
    message('Pre-processing step 1: select GWAS loci with silver standard genes')
    # gwas_loci = gwas_loci %>% mutate(identifier = merge_up(chromosome, start, end, trait))
    gwas_loci_with_gene = annotate_gwas_loci_with_gene(gwas_loci, gene_annotation) %>% mutate(trait_gene = pair_up_trait_and_gene(.data$trait, .data$gene)) %>% mutate(identifier = merge_up(.data$chromosome, .data$start, .data$end, .data$trait))
    gwas_loci_with_silver = gwas_loci_with_gene %>% filter(.data$trait_gene %in% trait_gene_silver$trait_gene) %>% select(-.data$gene) %>% distinct()
    message('Pre-processing step 2: select candidate gens overlapping selected GWAS loci')
    trait_gene_pairs_overlapping_gwas_loci_with_silver = gwas_loci_with_gene %>% filter(.data$identifier %in% gwas_loci_with_silver$identifier)
    score_table_subset = score_table %>% filter(pair_up_trait_and_gene(.data$trait, .data$gene) %in% trait_gene_pairs_overlapping_gwas_loci_with_silver$trait_gene)
    trait_gene_silver = trait_gene_silver %>% filter(.data$trait_gene %in% gwas_loci_with_silver$trait_gene)
  }
  message('# trait-gene pairs in score table after step 2: ', nrow(score_table_subset))
  if(nrow(score_table_subset) < 3) {
    message('Too few trait-gene pairs to work with. Exit')
    return()
  }

  # step 3
  # create candidate trait/gene pair pool
  # and silver standard trait/gene pair pool
  trait_gene_candidate_pool = pair_up_trait_and_gene(score_table_subset$trait, score_table_subset$gene)
  trait_gene_silver_pool = trait_gene_silver$trait_gene

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
  roc_aucs = list()
  for(i in score_cols) {
    roc_curve = gen_roc_curve(true_genes = trait_gene_silver_pool, gene = trait_gene_candidate_pool, score = score_table_subset[, i])
    roc_curves[[length(roc_curves) + 1]] = roc_curve %>% mutate(score = i)
    roc_aucs[[length(roc_aucs) + 1]] = compute_auc(roc_curves[[length(roc_curves)]]) %>% mutate(score = i)
  }
  roc_curves = do.call(rbind, roc_curves)
  roc_aucs = do.call(rbind, roc_aucs)

  roc = roc_curves %>% ggplot() + geom_path(aes_string(x = 'fpr', y = 'tpr', color = 'score')) + coord_equal()
  return(list(pr = pr, roc = roc, roc_auc = roc_aucs))
}

pair_up_trait_and_gene = function(trait, gene) {
  paste(trait, gene, sep = '-x-')
}

merge_up = function(...) {
  paste(..., sep = '-x-')
}

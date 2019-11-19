#' Silver standard based on OMIM <-> HPo/pheocode mapping
#'
#' The MIM phenotype to MIM gene map is obtained from genemap2.txt provided by OMIM database (accessed at 07/2019).
#'
#' @docType data
#'
#' @usage data(omim_based_silver_standard)
#'
#' @format a list of:
#'   1. table: a data.frame with
#'      column gene (with Ensembl ID, e.g. ENSG00000133019) ,
#'      column HPO (HPO term in the format HPO:0000126),
#'      column phecode (e.g. 595);
#'   2. script_info: meta information on how data set is generated.
#'
#' @keywords datasets
#'
#'
#'
#' @source \href{https://github.com/hakyimlab/gold-standard/blob/master/scripts/phecode_to_omim_to_gene.sh}{Source code generating this silver standard dataset}
#'
#' @examples
#' data(omim_based_silver_standard)
#' head(omim_based_silver_standard$table)
#'
"omim_based_silver_standard"

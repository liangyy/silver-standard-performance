#' Silver standard based on OMIM <-> HPO/pheocode mapping
#'
#' The MIM phenotype to MIM gene map is obtained from genemap2.txt
#' provided by OMIM database (accessed at 07/2019).
#' The MIM phenotypes to HPO/phecode map is provided by Lisa Bastarache.
#'
#' @docType data
#'
#' @usage data(omim_based_silver_standard)
#'
#' @format a list of:
#'   1. table: a data.frame with
#'      column gene (with Ensembl ID, e.g. ENSG00000133019) ,
#'      column HPO (HPO term in the format HP:0000126),
#'      column phecode (e.g. 595),
#'      column EFO (e.g. EFO:0004238);
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


#' Silver standard based on Orphanet <-> HPO/pheocode mapping
#'
#' The Orphanet phenotype to Orphanet gene map and
#' the Orphanet phenotype to phecode/HPO map
#' are provided by Lisa Bastarache.
#'
#' @docType data
#'
#' @usage data(orphanet_based_silver_standard)
#'
#' @format a list of:
#'   1. table: a data.frame with
#'      column gene (with Ensembl ID, e.g. ENSG00000084754) ,
#'      column HPO (HPO term in the format HP:0000512),
#'      column phecode (e.g. 346),
#'      column EFO (e.g. EFO:0004207);
#'   2. script_info: meta information on how data set is generated.
#'
#' @keywords datasets
#'
#'
#'
#' @source \href{https://github.com/hakyimlab/gold-standard/blob/master/scripts/phecode_to_orphanet_to_gene.sh}{Source code generating this silver standard dataset}
#'
#' @examples
#' data(orphanet_based_silver_standard)
#' head(orphanet_based_silver_standard$table)
#'
"orphanet_based_silver_standard"

#' Silver standard based on rare variant association
#'
#' Here, we focused on rare variant association evidence reported on height and lipid traits:
#'   1. low-density lipid cholesterol (LDL),
#'   2. high-density lipid cholesterol (HDL),
#'   3. triglycerides (TG),
#'   4. total cholesterol levels (high cholesterol)
#' In particular, we collected significant coding/splicing variants reported previously \href{https://www.nature.com/articles/nature21039}{Marouli et al 2017}
#' and kept variants with effect allele frequency < 0.01
#' (Table S6: ExomeChip variants with Pdiscovery <2e-07 in the European-ancestry meta-analysis (N=381,625)).
#' Similarly, we collected significant variants reported by \href{https://www.nature.com/articles/ng.3977}{Liu et al 2017}
#' (Table S12: Association Results for 444 independently associated variants with lipid traits) and filtered out variants with minor allele frequency < 0.01.
#' For the whole-exome sequencing study conducted in Finnish isolates \href{https://www.nature.com/articles/s41586-019-1457-z}{Locke et al 2019},
#' we extracted significant genes identified by a gene-based test using protein truncating variants
#' (Table S9: Gene-based associations from aggregate testing with EMMAX SKAT-O with P<3.88E-6)
#' and significant variants (Table S7: A review of all variants that pass unconditional threshold of P<5E-07 for at least one trait)
#' with gnomAD MAF < 0.01.


#'
#' @docType data
#'
#' @usage data(rare_variant_based_silver_standard)
#'
#' @format a list of:
#'   1. table: a data.frame with
#'      column gene (with Ensembl ID, e.g. ENSG00000084754) ,
#'      column HPO (HPO term in the format HP:0000512),
#'      column EFO (e.g. EFO:0004207),
#'      column trait (trait name, e.g. trait)
#'   2. script_info: meta information on how data set is generated.
#'
#' @keywords datasets
#'
#'
#'
#' @source \href{https://bitbucket.org/yanyul/rotation-at-imlab/src/master/analysis/fdr_power_specificity/data/gen_ewas_rare_variant_gene_list.R}{Source code generating this silver standard dataset}
#'
#' @examples
#' data(rare_variant_based_silver_standard)
#' head(rare_variant_based_silver_standard$table)
#'
"rare_variant_based_silver_standard"

#' Gene annotation (position and gene name etc)
#'
#' Gene annotation obtained from gencode v26 (formatted by Alvaro Barbeira).
#' It uses hg38 genomic coordinate.
#'
#' @docType data
#'
#' @usage data(gene_annotation_gencode_v26_hg38)
#'
#' @format a list of :
#'   1. gene_annotation, data.frame with columns:
#'     1.1. gene_id for Ensambl ID (no dot)
#'     1.2. gene_name for the "conventional" name of the gene
#'     1.3. gene_type for type of gene, see \href{https://www.gencodegenes.org/pages/biotypes.html}{here} for more information
#'     1.4. chromosome
#'     1.5. start (base-1, included)
#'     1.6. end (base-1, included)
#'     1.7. strand: strand of the transcript (optional)
#'   2. meta_information, a list of
#'     2.1 genome_assembly_version
#'     2.2 notes
#'
#'
#' @keywords datasets
#'
#'
#'
#'
#' @examples
#' data(gene_annotation_gencode_v26_hg38)
#'
"gene_annotation_gencode_v26_hg38"


#' LD block annotation (for European population, in hg38)
#'
#' Obtained from \href{https://academic.oup.com/bioinformatics/article/32/2/283/1743626}{Berisa et al}.
#' Formatted by Alvaro Barbeira
#'
#'
#'
#' @docType data
#'
#' @usage data(ld_block_pickrell_eur_b38)
#'
#' @format a list of
#'   1. ld_block, a data.frame with columns (each row is an LD block):
#'     1.1. region_name name of the LD block
#'     1.2. chromosome
#'     1.3. start (base-1, included)
#'     1.4. end (base-1, included)
#'   2. meta_information, a list of
#'     2.1 genome_assembly_version
#'     2.2 notes
#'
#' @keywords datasets
#'
#'
#'
#'
#' @examples
#' data(ld_block_pickrell_eur_b38)
#'
"ld_block_pickrell_eur_b38"

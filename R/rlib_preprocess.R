#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges findOverlaps
#' @importFrom IRanges IRanges
annotate_gwas_loci_with_gene = function(gwas_loci, gene_annotation, gene_column = 'gene_id') {
  gwas_loci_distinct = gwas_loci %>% select(.data$chromosome, .data$start, .data$end, .data$trait) %>% distinct()
  gr = GRanges(
    seqnames = gwas_loci_distinct$chromosome,
    ranges = IRanges(
      start = gwas_loci_distinct$start,
      end = gwas_loci_distinct$end
    )
  )
  if(nrow(gene_annotation) == 0) {
    message('No gene left after filtering by gene type. Please try again')
    quit()
  }
  gr2 = GRanges(
    seqnames = gene_annotation$chromosome,
    ranges = IRanges(
      start = gene_annotation$start,
      end = gene_annotation$end,
    )
  )
  overlap_idx = as.data.frame(findOverlaps(gr, gr2))
  out = gwas_loci_distinct[overlap_idx$queryHits, ]
  out = out %>% mutate(gene = gene_annotation[overlap_idx$subjectHits, gene_column])
  out
}

#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges findOverlaps
#' @importFrom IRanges IRanges
gwas_hit_to_gwas_loci_by_ld_block = function(gwas_hit, ld_block) {
  message('The genome assembly version of ld block is ', ld_block$meta_information$genome_assembly_version, '. Please make sure your gwas_hit matches the version!')
  shared_colnames = intersect(colnames(gwas_hit), colnames(ld_block))
  if(length(shared_colnames) > 0) {
    message('gwas_hit and ld_block share column name: ', paste(shared_colnames), ', please change')
  }
  ld_block = ld_block$ld_block
  gr = GRanges(
    seqnames = gwas_hit$chromosome,
    ranges = IRanges(
      start = gwas_hit$position,
      end = gwas_hit$position
    )
  )
  gr2 = GRanges(
    seqnames = ld_block$chromosome,
    ranges = IRanges(
      start = ld_block$start,
      end = ld_block$end,
    )
  )
  overlap_idx = as.data.frame(findOverlaps(gr, gr2))
  out = gwas_hit[overlap_idx$queryHits, ]
  out = cbind(out, ld_block[overlap_idx$subjectHits, ] %>% select(-.data$chromosome))
  out
}

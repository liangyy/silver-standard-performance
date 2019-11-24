# script to generate data

# data/omim_based_silver_standard.rda
{
  # download from https://github.com/hakyimlab/gold-standard/blob/master/data/silver_phecode_to_mim.rds
  omim_based_silver_standard = readRDS('~/Downloads/silver_phecode_to_mim.rds')
  save(omim_based_silver_standard, file = 'data/omim_based_silver_standard.rda', compress = 'xz')
}

# data/orphanet_based_silver_standard.rda
{
  # download from https://github.com/hakyimlab/gold-standard/blob/master/data/silver_phecode_to_orphanet.rds
  orphanet_based_silver_standard = readRDS('~/Downloads/silver_phecode_to_orphanet.rds')
  save(orphanet_based_silver_standard, file = 'data/orphanet_based_silver_standard.rda', compress = 'xz')
}

# data/rare_variant_based_silver_standard.rda
{
  rare_var = read.table('https://bitbucket.org/yanyul/rotation-at-imlab/raw/85a3fbe8f08df7c67265fed69569b7ea554d4e12/analysis/fdr_power_specificity/data/ewas/trait-to-gene-ewas-maf_0.01.full-list.txt', header = T, stringsAsFactors = F)
  map_dic = list(
    UKB_50_Standing_height = 'EFO:0004339',
    GLGC_Mc_LDL = 'EFO:0004611',
    GLGC_Mc_HDL = 'EFO:0004612',
    GLGC_Mc_TG = 'EFO:0004530',
    UKB_20002_1473_self_reported_high_cholesterol = NA
  )
  map_dic2 = list(
    UKB_50_Standing_height = NA,
    GLGC_Mc_LDL = NA,
    GLGC_Mc_HDL = NA,
    GLGC_Mc_TG = NA,
    UKB_20002_1473_self_reported_high_cholesterol = 'HP:0003119'
  )
  map_dic3 = list(
    UKB_50_Standing_height = 'height',
    GLGC_Mc_LDL = 'LDL',
    GLGC_Mc_HDL = 'HDL',
    GLGC_Mc_TG = 'TG',
    UKB_20002_1473_self_reported_high_cholesterol = 'high cholesterol'
  )
  new = c()
  new2 = c()
  new3 = c()
  for(i in 1 : nrow(rare_var)) {
    new = c(new, map_dic[[rare_var$trait[i]]])
    new2 = c(new2, map_dic2[[rare_var$trait[i]]])
    new3 = c(new3, map_dic3[[rare_var$trait[i]]])
  }
  rare_var$trait_name = new3
  rare_var$HPO = new2
  rare_var$EFO = new
  rare_var = rare_var %>% select(gene, HPO, EFO, trait_name)
  rare_variant_based_silver_standard = list(table = rare_var, script_info = 'gen_ewas_rare_variant_gene_list.R')
  save(rare_variant_based_silver_standard, file = 'data/rare_variant_based_silver_standard.rda', compress = 'xz')
}

# data/gene_annotation_gencode_v26_hg38.rds
{
  gene_annotation_gencode_v26_hg38 = read.table('https://bitbucket.org/yanyul/rotation-at-imlab/raw/85a3fbe8f08df7c67265fed69569b7ea554d4e12/data/annotations_gencode_v26.tsv', header = T, stringsAsFactors = F)
  gene_annotation_gencode_v26_hg38 = list(
    gene_annotation = gene_annotation_gencode_v26_hg38,
    meta_information = list(
      genome_assembly_version = 'hg38',
      notes = 'Gene annotation obtained from gencode v26 (formatted by Alvaro Barbeira). It uses hg38 genomic coordinate.'
    )
  )
  save(gene_annotation_gencode_v26_hg38, file = 'data/gene_annotation_gencode_v26_hg38.rda', compress = 'xz')
}

# data/ld_block_pickrell_eur_b38.rda
{
  ld_block_pickrell_eur_b38 = read.table('https://bitbucket.org/yanyul/rotation-at-imlab/raw/85a3fbe8f08df7c67265fed69569b7ea554d4e12/analysis/allelic_heterogeneity/data/ld_independent_regions.txt', header = T, stringsAsFactors = F)
  ld_block_pickrell_eur_b38 = list(
    ld_block = ld_block_pickrell_eur_b38,
    meta_information = list(
      genome_assembly_version = 'hg38',
      notes = 'LD block annotation (for European population) obtained from Berisa et al (https://academic.oup.com/bioinformatics/article/32/2/283/1743626). Formatted by Alvaro Barbeira'
    )
  )
  save(ld_block_pickrell_eur_b38, file = 'data/ld_block_pickrell_eur_b38.rda', compress = 'xz')
}

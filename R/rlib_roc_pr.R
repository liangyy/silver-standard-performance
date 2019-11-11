#' Calculate power and FDR
#'
#' Given the set of true signals and the set of proposed signals, calculate power and false discovery rate
#'
#' @param discovered_genes a set of proposed signals (their IDs as string)
#' @param true_genes a set of true signals (their IDs as string)
#'
#' @return data.frame with power and fdr as columns along with nTP (number of true positives), nP (number of positives, i.e. proposed signals), and nT (number of true signals)
#'
#' @examples
#' compute_power_and_fdr(
#'   discovered_genes = c('g1', 'g3', 'g4'),
#'   true_genes = c('g1', 'g2')
#' )
#'
#' @export
compute_power_and_fdr = function(discovered_genes, true_genes) {
  true_positive_genes = true_genes[ true_genes %in% discovered_genes ]
  nTP = length(true_positive_genes)
  nT = length(true_genes)
  nP = length(discovered_genes)
  if(nP == 0) {
    nP = 1
  }
  data.frame(power = nTP / nT, fdr = 1 - nTP / nP, nTP = nTP, nP = nP, nT = nT)
}

#' Calculate TPR and FPR
#'
#' Given the set of all candidate genes, true signals, and proposed signals as three separate sets, calculate true positive rate (TPR) and false positive rate (FPR)
#'
#' @param tested_genes a set of candidate genes (their IDs as string; it should include both discovered_genes and true_genes)
#' @param discovered_genes a set of proposed signals (their IDs as string)
#' @param true_genes a set of true signals (their IDs as string)
#'
#' @return data.frame with tpr and fpr as columns along with nTP (number of true positives), nT (number of true signals), nFP (number of false positives, i.e. proposed but not true signals), and nTP (number of proposed and true signals)
#'
#' @examples
#' compute_tpr_and_fpr(
#'   tested_genes = c('g1', 'g2', 'g3', 'g4', 'g5', 'g6'),
#'   discovered_genes = c('g1', 'g3', 'g4'),
#'   true_genes = c('g1', 'g2')
#' )
#'
#' @export
compute_tpr_and_fpr = function(tested_genes, discovered_genes, true_genes) {
  true_positive_genes = true_genes[ true_genes %in% discovered_genes ]
  nTP = length(true_positive_genes)
  nT = length(true_genes)

  all_genes = union(tested_genes, true_genes)
  not_true_genes = all_genes[!all_genes %in% true_genes]
  nNT = length(not_true_genes)
  false_positive_genes = not_true_genes[ not_true_genes %in% discovered_genes ]
  nFP = length(false_positive_genes)
  if(nNT == 0) {
    nNT = 1
  }
  data.frame(tpr = nTP / nT, fpr = nFP / nNT, nTP = nTP, nT = nT, nFP = nFP, nNT = nNT)
}

#' Check if the candidates should be proposed
#'
#' Given the scores of the candidates and a threshold to apply along with method to propose, e.g. 'gt' means propose if score is greater than the threshold and 'lt' means propose if score is less than the threshold. It handle 'tie', when tie = TRUE, the candidate with score equals to threshold will be proposed
#'
#' @param score the scores of candidates
#' @param threshold threshold to propose a signal
#' @param method method = 'gt' to propose score greater than threshold; method = 'lt' to propose score less than threshold
#' @param tie if propose the candidates which have scores equal to threshold (default TRUE)
#'
#' @return a binary vector of the same length as the input score vector indicating if the corresponding score entry is proposed
#'
#' @examples
#' is_sig(score = 1:10, threshold = 3, method = 'gt', tie = TRUE)
#'
#' @export
is_sig = function(score, threshold, method, tie = T) {
  if(method == 'gt') {
    if(tie == T) {
      return(score >= threshold)
    } else {
      return(score > threshold)
    }
  } else if(method == 'lt') {
    if(tie == T) {
      return(score <= threshold)
    } else {
      return(score < threshold)
    }
  } else {
    return(NA)
  }
}

#' Generate precision-recall curve
#'
#' Given the set of true signals and the set of candidate genes along with their scores and proposing method, generate precision-recall curve, which is a sequence of precision and recall measured at a set of cutoffs
#'
#' @param true_genes the set of true signals
#' @param gene the set of all candidate signals
#' @param score the scores of the candidates (the same length as gene)
#' @param method method = 'gt' to propose score greater than threshold; method = 'lt' to propose score less than threshold
#' @param cutoff cutoff = NULL will use all unique scores from true signals as cutoffs (essentially this is the finest grid we can get)
#'
#' @return a data.frame with all columns of compute_power_and_fdr's output along with cutoff column telling which cutoff is used to obtain the results of the row and precision (same as 1 - fdr) and recall (same as power)
#'
#' @examples
#' gen_fdr_power_curve(
#'   true_genes = c('g1', 'g3', 'g4'),
#'   gene = c('g1', 'g2', 'g3', 'g4', 'g5', 'g6'),
#'   score = 1:6,
#'   method = 'gt',
#'   cutoff = NULL
#' )
#'
#' @export
#' @import dplyr
gen_fdr_power_curve = function(true_genes, gene, score, method = 'gt', cutoff = NULL) {
  # set fdr and power to NULL to avoid NOTE in devtools::check()
  fdr <- power <- NULL
  if(is.null(cutoff)) {
    true_cutoffs = sort(unique(score[gene %in% true_genes]))
  } else {
    true_cutoffs = cutoff
  }
  if(method == 'lt') {
    true_cutoffs = sort(true_cutoffs, decreasing = T)
  }
  df_curve = data.frame()
  for(i in true_cutoffs) {
    positive_genes = gene[is_sig(score, i, method)]
    sub = compute_power_and_fdr(positive_genes, true_genes)
    sub$cutoff = i
    sub$include_tie = T
    df_curve = rbind(df_curve, sub)
    positive_genes = gene[is_sig(score, i, method, tie = F)]
    sub = compute_power_and_fdr(positive_genes, true_genes)
    sub$cutoff = i
    sub$include_tie = F
    df_curve = rbind(df_curve, sub)
  }
  df_curve = df_curve %>% mutate(precision = 1 - fdr, recall = power)
  df_curve$precision[df_curve$nP == 0] = 1
  df_curve
}

#' Generate receiver operating characteristic curve
#'
#' Given the set of true signals and the set of candidate genes along with their scores and proposing method, generate ROC curve, which is a sequence of FPR and TPR measured at a set of cutoffs
#'
#' @param true_genes the set of true signals
#' @param gene the set of all candidate signals
#' @param score the scores of the candidates (the same length as gene)
#' @param method method = 'gt' to propose score greater than threshold; method = 'lt' to propose score less than threshold
#' @param cutoff cutoff = NULL will use all unique scores from true signals as cutoffs (essentially this is the finest grid we can get)
#'
#' @return a data.frame with all columns of compute_tpr_and_fpr's output along with cutoff column telling which cutoff is used to obtain the results of the row
#'
#' @examples
#' gen_roc_curve(
#'   true_genes = c('g1', 'g3', 'g4'),
#'   gene = c('g1', 'g2', 'g3', 'g4', 'g5', 'g6'),
#'   score = 1:6,
#'   method = 'gt',
#'   cutoff = NULL
#' )
#'
#' @export
gen_roc_curve = function(true_genes, gene, score, method = 'gt', cutoff = NULL) {
  if(is.null(cutoff)) {
    true_cutoffs = sort(unique(score[gene %in% true_genes]))
  } else {
    true_cutoffs = cutoff
  }
  if(method == 'lt') {
    true_cutoffs = sort(true_cutoffs, decreasing = T)
  }
  df_curve = data.frame()
  for(i in true_cutoffs) {
    positive_genes = gene[is_sig(score, i, method)]
    sub = compute_tpr_and_fpr(gene, positive_genes, true_genes)
    sub$cutoff = i
    sub$include_tie = T
    df_curve = rbind(df_curve, sub)
    positive_genes = gene[is_sig(score, i, method, tie = F)]
    sub = compute_tpr_and_fpr(gene, positive_genes, true_genes)
    sub$cutoff = i
    sub$include_tie = F
    df_curve = rbind(df_curve, sub)
  }
  # sub = compute_tpr_and_fpr(gene, union(gene, true_genes), true_genes)
  # sub$cutoff = NA
  # df_curve = rbind(df_curve, sub)
  df_curve
}

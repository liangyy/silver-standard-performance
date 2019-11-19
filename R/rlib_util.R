#' Calculate AUC for ROC
#'
#' Take the ROC curve output from compute_tpr_and_fpr function and compute the AUC of the curve
#'
#' @param curve output of compute_tpr_and_fpr
#'
#' @return a data.frame with one row and one column named roc_auc
#'
#'
compute_auc = function(curve) {
  # f = curve %>% filter(nTP > 0) %>% group_by(nTP) %>% summarize(y = mean(tpr), x = fpr[1] - fpr[2])
  o1 = rbind(c(1, 1), curve %>% select(.data$tpr, .data$fpr))
  o2 = rbind(curve %>% select(.data$tpr, .data$fpr), c(0, 0))
  # message('e')
  y = data.frame(x1 = o1$fpr, x2 = o2$fpr, y1 = o1$tpr, y2 = o2$tpr)
  e = y %>% mutate(s = (.data$x1 - .data$x2) * (.data$y1 + .data$y2) / 2)
  data.frame(roc_auc = sum(e$s))
}

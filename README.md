# prroc

A light and flexible R package for plotting PR/ROC curves

The initial effort is at [this gist post](https://gist.github.com/liangyy/6d4314dbc238236731e134abef2484f4#file-rlib_roc_pr-r).

# Minimal tutorial

```
library(dplyr)
library(ggplot2)
library(ROCR)
data(ROCR.simple)

source('rlib_roc_pr.R')

# add names for each sample 
ROCR.simple$gene_name = paste('gene', 1 : length(ROCR.simple$predictions))

# extract gene names for real signal
true_genes = ROCR.simple$gene_name[which(ROCR.simple$labels == 1)]

# PR curve
## gen_fdr_power_curve
## true_genes: a vector of gene name of the real signal
## gene: a vector of gene name for all candidate genes
## score: a vector of numeric score for all candidate genes
## method: 'gt' if the higher score means more likely to be signal; otherwise use 'lt'
## cutoff: a vector of numeric score cutoffs. If it is NULL, the function will return the whole curve. Otherwise, it will return FDR/Precision and Power/Recall at the cutoffs 

df = gen_fdr_power_curve(true_genes = true_genes, gene = ROCR.simple$gene_name, score = ROCR.simple$predictions, method = 'gt', cutoff = NULL)
df[-nrow(df), ] %>% ggplot() + geom_path(aes(x = recall, y = precision))

pred <- prediction( ROCR.simple$predictions, ROCR.simple$labels)
perf1 <- performance(pred, "prec", "rec")
plot(perf1)


# ROC curve
## gen_roc_curve
## same input 

df = gen_roc_curve(true_genes = true_genes, gene = ROCR.simple$gene_name, score = ROCR.simple$predictions, method = 'gt', cutoff = NULL)
df %>% ggplot() + geom_path(aes(x = fpr, y = tpr))

pred <- prediction( ROCR.simple$predictions, ROCR.simple$labels)
perf1 <- performance(pred, "tpr", "fpr")
plot(perf1)
```

See a detailed example at [this gist](https://gist.github.com/liangyy/6d4314dbc238236731e134abef2484f4#file-test_roc_pr-r).

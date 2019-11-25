# silver-standard-performance


## Overview 

**SilverStandardPerformance** is an light-weight R package. 
Its main functionality is to evaluate the performance of gene prioritization methods for complex traits, *e.g.* [coloc](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004383), [PrediXcan](https://www.nature.com/articles/ng.3367), [smr](https://www.nature.com/articles/ng.3538), [enloc](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006646), etc. 
The evaluation simply relies on an external "ground truth" and the performance is simply how well the method can prioritize the "ground truth" trait/gene pairs relative to other trait/gene pairs.
In other word, we simply frame the evaluation problem as a typical problem for evaluating the performance of classifier. 
The only difference is that we don't have the ideal "ground truth".   
Instead, we build several "silver standard"s as the imperfect surrogate of "ground truth".
 

# silver-standard-performance


## Overview 

**SilverStandardPerformance** is an light-weight R package. 
Its main functionality is to evaluate the performance of gene prioritization methods for complex traits, *e.g.* [coloc](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004383), [PrediXcan](https://www.nature.com/articles/ng.3367), [smr](https://www.nature.com/articles/ng.3538), [enloc](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006646), etc. 
The evaluation simply relies on an external "ground truth" and the performance is simply how well the method can prioritize the "ground truth" trait/gene pairs relative to other trait/gene pairs.
In other word, we simply frame the evaluation problem as a typical problem for evaluating the performance of a binary classifier. 
The only difference is that we don't have the ideal "ground truth".   
Instead, we build several "silver standard"s as the imperfect surrogate of "ground truth".

## Rationales behind the silver standard

For the purpose of evaluating how well the method prioritizes causal genes for complex traits, the ideal "silver standard" should list the (likely) causal genes for each of the complex traits.
An additional requirement for the "silver standard" is that it should be _independent_ (apparently it is hard to achieve) to GWAS since many of the prioritization methods make use of GWAS results.

Following this rationale, we composed a rare variant association based silver standard which included height and lipid traits (LDL, HDL, TG, and high cholesterol), which can be loaded by typing 
```
load("rare_variant_based_silver_standard")
```
in R console.
However, as you may notice, this silver standard is very limited in terms of the number of both traits and genes being included.
Also, the result is still association-based rather than causality and many of the included trait/gene pairs are lack of experimental validation.
In fact, the gene-to-complex-trait information is currently very sparse (and that's why we care about the prioritization of causal genes for complex traits).


On the other hand, there are databases for causal genes of Mendelian diseases and rare diseases.
Since these phenotypes are mostly monogenic/oligogenic, the causal gene to phenotype mapping is of high confidence. 
To make use of these great resources, we make extra hypothesis that:
**The causal genes of the monogenic/oligogenic rare/Mendelian diseases can have milder effect on the relevant complex traits.** So, if we can find the 'relevant' complex traits for these rare/Mendelian diseases, we can build another trait/gene silver standard by extrapolating to complex trait. 
In this case, the big effect (due to deletion or coding variation) in rare/Mendelian diseases is extrapolate to the mild effect (probably via non-coding variation) in the relevant complex trait (see [GTEx GWAS paper](https://www.biorxiv.org/content/10.1101/814350v1) for more discussion on this topic).  

To do so, we mapped complex traits to rare/Mendelian phenotypes with the use of phenotype ontology.
For complex trait phenotypes, we use

* EFO [paper](https://academic.oup.com/bioinformatics/article/26/8/1112/208992) [website](https://www.ebi.ac.uk/efo/). _Comment: As used in [GWAS Catalog](https://www.ebi.ac.uk/gwas/). It describes quantitative traits, like height._
* HPO [paper](https://academic.oup.com/nar/article/42/D1/D966/1042793) [website](https://hpo.jax.org/app/). _Comment: It mostly describes human abnormal phenotypes._
* phecode [paper](https://academic.oup.com/bioinformatics/article/26/9/1205/201211) [website](https://phewascatalog.org/phecodes). _Comment: it is derived from ICD code._

Besides, the map among EFO, HPO, and phecode terms is achievable.
Phecode has been mapped to HPO (see [this paper](https://science.sciencemag.org/content/359/6381/1233) for more details).
Phecode has also been mapped to GWAS Catalog traits (see [this paper](https://www.nature.com/articles/nbt.2749) for more details) which are further mapped to EFO terms. 

For establishing the map between complex traits and rare/Mendelian phenotypes, HPO provides annotations to rare/Mendelian databases such as OMIM and Orphanet. 
Using all of these resources, we obtained maps among HPO, phecode, and MIM/Orphanet phenotypes and Orphanet to gene map which were curated and provided by Lisa Bastarache (great thanks!). 
We then went through some extra steps like formatting and mapping MIM phenotypes to MIM genes (information on the related scripts are at [here](https://github.com/hakyimlab/gold-standard/tree/master/data#silver-standard-rdss)). 
These gave rise to the OMIM based silver standard and Orphanet based silver standard.
You can load these two datasets by typing 

```
data("omim_based_silver_standard")
data("orphanet_based_silver_standard")
```

As you may notice, these silver standards are far from perfect. 
We should use and interpret the results with cautious while keeping in mind the imperfectness and potential bias of the silver standard. 
For your reference, we discuss some of the limitations of the current silver standard at [here](https://liangyy.github.io/silver-standard-performance/example_with_preprocessing.html#2_about_pre-processing_step_in_analysis_pipeline).



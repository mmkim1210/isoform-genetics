# isoform-genetics

This repository contains code and instructions for heritability and genetic correlation analyses, using the <a href="https://github.com/Hua-Zhou/MultiResponseVarianceComponentModels.jl"><img src="assets/MRVCModels-logo.svg" width="40em"> MultiResponseVarianceComponentModels.jl </a> <a href="https://julialang.org"><img src="https://julialang.org/assets/infra/julia.ico" width="10em"> Julia </a>package and applying it to the [PsychENCODE brain gene and isoform expression data](http://resource.psychencode.org/). There is a great number of primary literature, review papers, and educational materials that explain in broad strokes the heuristics of heritability and genetic correlation analyses. Briefly, heritability captures the degree of genetic effects, while genetic correlation captures the extent of shared genetic influences or genetic overlap and pleiotropy. If you are interested more in technical details with clarity in presentation and mathematical notation, please take a look at the **Methods** section of the associated paper [Kim et al. 2022](#references). All analyses herein were conducted solely using <a href="https://julialang.org"><img src="https://julialang.org/assets/infra/julia.ico" width="10em"> Julia</a>. If you have any questions, let me know via my email minsookim@mednet.ucla.edu.

## Getting started
[To install necessary packages and activate a separate environment](https://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project), open <a href="https://julialang.org"><img src="https://julialang.org/assets/infra/julia.ico" width="10em"> Julia</a> within the directory and type:
```julia
julia> ]
(@v1.10) pkg> activate .
(isoform-genetics) pkg> instantiate
```
## Required data
Some data like [GENCODE](https://www.gencodegenes.org/human/) is automatically downloaded when it is missing, whereas other data like the ones below are not publicly accessible and need to be made available before running any <a href="https://julialang.org"><img src="https://julialang.org/assets/infra/julia.ico" width="10em"> Julia </a>code. See below for some notes on these required data.

    data
    ├── expression
    │   ├── PsychENCODE-EUR-gene.BED.gz                                                 # normalized PsychENCODE gene expression
    │   ├── PsychENCODE-EUR-isoform.BED.gz                                              # normalized PsychENCODE isoform expression
    │   └── PsychENCODE-EUR-covariates.tsv                                              # covariates for mean (or fixed) effects
    ├── genotype
    │   └── Capstone4.HRC.European.unique.frontal.nochr.filter.unrelated.{bed,bim,fam}  # PsychENCODE genotype data
    ├── 1kg
    │   └── kgp.eur.maf0.05.{bed,bim,fam}                                               # 1000 Genomes data subsetted to European individuals
    └── gwas                                                                            # GWAS summary statistics

**Couple notes:**
- `PsychENCODE-EUR-gene.BED.gz` and `PsychENCODE-EUR-isoform.BED.gz`: [RNA-seq reads were previously aligned](https://www.science.org/doi/10.1126/science.aat8127) to the hg19 reference genome with STAR 2.4.2a and gene and isoform-level quantifications calculated using RSEM v1.2.29. Genes and isoforms were filtered to include those with TPM > 0.1 in at least 25% of samples. Gene and isoform expression were separately normalized using TMM normalization in edgeR and log2-transformed. RNA-seq data was also restricted to frontal cortex samples from European individuals as well as genes and isoforms belonging to autosomal chromosomes, resulting in a total of **24,905 genes and 93,293 isoforms** based on GENCODE v19 annotation.
- `PsychENCODE-EUR-covariates.tsv`: The same set of known biological and technical covariates were used for mean or fixed effects, which include age, age<sup>2</sup>, study, sex, diagnosis, RNA integrity number (RIN), RIN<sup>2</sup>, post-mortem interval (PMI), 24 sequencing principal components (PCs), and 5 genetic PCs.
- `Capstone4.HRC.European.unique.frontal.nochr.filter.unrelated.{bed,bim,fam}`: [Genotype data were previously harmonized](https://www.science.org/doi/10.1126/science.aat8464) through phasing and imputation with the Haplotype Reference Consortium (HRC) reference panel. We focused on 860 unique European individuals with matching genotype and frontal cortex RNA-seq data. We started with 5,312,508 HRC imputed SNPs and filtered for SNPs with minor allele frequency (MAF) > 0.01, genotype and individual missingness rate < 0.05, and Hardy-Weinberg equilibrium (HWE) *P* values > 10<sup>-6</sup>. Five pairs of individuals had classic genetic relationship matrix (GRM) values > 0.05 when using all filtered SNPs, while 647 pairs of individuals had GRM values > 0.025. We kept one individual from each of five pairs and only SNPs belonging to autosomal chromosomes, resulting in a total of **855 unrelated European individuals and 4,685,674 SNPs** for downstream analyses.
- `kgp.eur.maf0.05.{bed,bim,fam}`: Linkage disequilibrium (LD) reference panel based on individuals of European ancestry in the [1000 Genomes Project](https://www.internationalgenome.org/) was generated using https://github.com/mmkim1210/1kg. This data was only used for visualization purpose, and not for actual variance components analysis.
- `gwas`: [Multiple GWAS summary statistics](https://github.com/mmkim1210/GeneticsMakie.jl/blob/master/src/gwas.jl) were downloaded and harmonized using `mungesumstats!` function of the <a href="https://github.com/mmkim1210/GeneticsMakie.jl"><img src="assets/GM-logo.svg" width="30em"> GeneticsMakie.jl </a> <a href="https://julialang.org"><img src="https://julialang.org/assets/infra/julia.ico" width="10em"> Julia </a>package. This data was not required for variance components analysis.

## Variance components analysis
<a href="https://hua-zhou.github.io/MultiResponseVarianceComponentModels.jl/dev/"><img src="assets/MRVCModels-logo.svg" width="40em"> MultiResponseVarianceComponentModels.jl </a>is a <a href="https://julialang.org"><img src="https://julialang.org/assets/infra/julia.ico" width="10em"> Julia </a>package for fitting and testing multivariate response variance components linear mixed models of form

$$\text{vec}\ \boldsymbol{Y} \sim \mathcal{N}(\text{vec}(\boldsymbol{X} \boldsymbol{B}), \sum_{i=1}^m \boldsymbol{\Gamma}_i \otimes \boldsymbol{V}_i),$$

where $\boldsymbol{Y}$ and $\boldsymbol{X}$ are $n \times d$ response and  $n \times p$ predictor matrices, respectively, and $\boldsymbol{V}_1, \ldots, \boldsymbol{V}_m$ are $m$ known $n \times n$ positive semidefinite matrices. $\text{vec}\ \boldsymbol{Y}$ creates an $nd \times 1$ vector from $\boldsymbol{Y}$ by stacking its columns and $\otimes$ denotes the Kronecker product. The parameters of the model include $p \times d$ mean effects $\boldsymbol{B}$ and $d \times d$ variance components ($\boldsymbol{\Gamma}_1, \dots, \boldsymbol{\Gamma}_m$), which <a href="https://hua-zhou.github.io/MultiResponseVarianceComponentModels.jl/dev/"><img src="assets/MRVCModels-logo.svg" width="40em"> MultiResponseVarianceComponentModels.jl </a>estimates through either [minorization-maximization (MM)](https://en.wikipedia.org/wiki/MM_algorithm) or [expectation–maximization (EM)](https://en.wikipedia.org/wiki/Expectation–maximization_algorithm) algorithms. Note that univariate and bivariate response models are subsumed under this multivariate response model when $d$ is one and two, respectively.

Herein, we initially modelled human brain gene and isoform-level expression using univariate variance components linear mixed models with `PsychENCODE-EUR-covariates.tsv` as mean effects covariates and specified three variance components, two of which capture *cis*- and *trans*-SNP genetic effects. In other words, $n = 855$ and $m = 3$. We defined *cis*-SNPs as those within 1 Mb window of gene start and gene end sites, and *trans*-SNPs as all the other SNPs. Based on this definition, **24,754 genes and 93,030 isoforms** had non-zero *cis*-SNPs. The same set of *cis*-SNPs was used for a given gene and its constituent isoforms for direct comparison. The mean number of *cis*-SNPs were 3,264 and 3,274 for these genes and isoforms, respectively. Variance components parameters were estimated using both maximum likelihood (ML) and restricted maximum likelihood (REML) estimation. The maximum number of iterations was set to 3,000. In total, **22,965 genes and 89,926 isoforms** had converged estimates.

We fitted multivariate variance components linear mixed models as well for isoform-level expression. We similarly specified three variance components, one of which captures *cis*-SNP genetic effects and the other *trans* effects. We used the same set of *cis*-SNPs that were used for univariate models. To reduce computational burden and the number of variance components parameters that need to be estimated, given limited sample size of the PsychENCODE dataset, we ran the multivariate model for isoforms with significant heritability estimates in a univariate model at *P* value < 0.05. This meant modeling up to 23 isoforms or **$d = 23$**. For isoforms that are perfectly correlated in expression, we included only one isoform of the two. Variance components parameters were estimated using REML. Note that the MM algorithm is numerically stable to fit even the non-heritable isoforms, but we chose not to. Lastly, we also fitted pairwise bivariate variance components models with the same scheme as multivariate models. For each gene with at least two heritable isoforms, the model was fit to all pairwise combinations of isoforms. These analyses were conducted by running `./submit.sh` as follows.

The analyses will require at least $m(nd)^2$ storage space for double-precision floating-point numbers, which in our case corresponds to $3 \cdot (855 \cdot 23)^2 \cdot 8\text{B} \div 10^9 \approx 9.3 \text{GB}$. Accordingly, $48\text{GB}$ memory we requested in `./submit.sh` should have been more than sufficient.

```bash
qsub ./submit.sh
```

## Figures
Once the results are parsed with `./src/parse.jl`, **Figures 1-5** in [Kim et al. 2022](#references) can be reproduced by running `./src/fig{1,2,3,4,5}.jl`, respectively. Then minor edits can be made using `Illustrator` to finalize the figures as follows.

**Figure 1:**
<p align="center"><img width="100%" style="border-radius: 5px;" src="assets/fig1.png"></p>

**Figure 2:**
<p align="center"><img width="100%" style="border-radius: 5px;" src="assets/fig2.png"></p>

**Figure 3:**
<p align="center"><img width="100%" style="border-radius: 5px;" src="assets/fig3.png"></p>

**Figure 4:**
<p align="center"><img width="100%" style="border-radius: 5px;" src="assets/fig4.png"></p>

**Figure 5:**
<p align="center"><img width="100%" style="border-radius: 5px;" src="assets/fig5.png"></p>

<a name="references"/>

## References
- M. Kim, D.D. Vo, C.T. Jops, C. Wen, A. Patowary, A. Bhattacharya, C.X. Yap, H. Zhou, and M.J. Gandal: **Multivariate variance components analysis uncovers genetic architecture of brain isoform expression and novel psychiatric disease mechanisms** (2022) ([link](https://www.medrxiv.org/content/10.1101/2022.10.18.22281204v1))
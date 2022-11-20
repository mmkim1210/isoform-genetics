# isoform-genetics

## Getting started 
```julia
julia> ]
pkg> activate .
pkg> instantiate
```
## Required data
    data
    ├── expression
    │   ├── PsychENCODE-EUR-gene.BED.gz                                                 # normalized gene expression
    │   ├── PsychENCODE-EUR-isoform.BED.gz                                              # normalized isoform expression
    │   └── PsychENCODE-EUR-covariates.tsv                                              # covariates for mean (or fixed) effects
    ├── genotype
    │   └── Capstone4.HRC.European.unique.frontal.nochr.filter.unrelated.{bed,bim,fam}  # genotype data
    ├── 1kg
    │   └── kgp.eur.maf0.05.{bed,bim,fam}                                               # 1000 Genomes data subsetted to European individuals
    └── gwas                                                                            # gwas summary statistics

## Schematic
<p align="center"><img width="100%" style="border-radius: 5px;" src="schematic.png"></p>

## References
- M. Kim, D.D. Vo, C.T. Jops, C. Wen, A. Patowary, A Bhattacharya, C.X. Yap, H. Zhou, and M.J. Gandal: **Multivariate variance components analysis uncovers genetic architecture of brain isoform expression and novel psychiatric disease mechanisms** (2022) ([link](https://www.medrxiv.org/content/10.1101/2022.10.18.22281204v1))

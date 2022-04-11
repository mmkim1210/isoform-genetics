# isoform-genetics

## Getting started 
```julia
julia> ]
pkg> activate .
pkg> instantiate
```
## Required data
    data/
    ├── expression
    │   ├── PsychENCODE-EUR-gene.BED.gz                                                 # normalized gene expression
    │   ├── PsychENCODE-EUR-isoform.BED.gz                                              # normalized isoform expression
    │   └── PsychENCODE-EUR-covariates.tsv                                              # covariates for mean (or fixed) effects
    ├── genotype/
    │   └── Capstone4.HRC.European.unique.frontal.nochr.filter.unrelated.{bed,bim,fam}  # genotype data
    ├── 1kg/
    │   └── kgp.eur.maf0.05.{bed,bim,fam}                                               # 1000 Genomes data subsetted to European individuals
    └── gwas/                                                                           # gwas summary statistics

## References
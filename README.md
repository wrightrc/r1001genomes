# Access and analyze the 1001 genomes *Arabidopsis thaliana* resequencing dataset in R  

By: Morgan O. Hamm and R. Clay Wright

This web app and open source software package for the R programming language was built to facilitate exploration of the naturally occuring sequence diversity for your favorite gene or gene family in Arabidopsis thaliana. We hope this web app will help plant biologists without bioinformatics experience access the 1001 genomes project’s rich catalog of natural variation and formulate hypotheses and identify new alleles in their favorite gene family.

We are developing this package and web application as open source resources (MIT licenced), but plan to publish on the initial use of these resources in early 2019. If you make use of this package and wish to publish prior to our publication please contact us at wright[dot]clay[at]gmail[dot]com and consider establishing a formal collaboration. 

While the package and app are in development, you can install this package in your R instance and run the app on your computer using the following code:
```{r}
# You will need to install several bioconductor packages
install.packages("BiocManager")

library(BiocManager) # shouldn't be necessary
BiocManager::install()

pkgs <- c(
    "biomaRt",
    "Biostrings",
    "BSgenome",
    "BSgenome.Athaliana.TAIR.TAIR9",
    "DECIPHER",
    "GenomicFeatures",
    "VariantAnnotation",
    "vcfR",
    "XVector"
    )

ap.db <- available.packages(contrib.url(BiocManager::repositories()))
ap <- rownames(ap.db)

pkgs_to_install <- pkgs[pkgs %in% ap]

BiocManager::install(pkgs_to_install, update=FALSE, ask=FALSE)

# Then you can install 'r1001genomes' along with the required CRAN packages
require("devtools")
devtools::install_github('wrightrc/r1001genomes')
run1001genomes()
```



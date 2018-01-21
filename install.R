# install.R

url <- "http://bioconductor.org/packages/3.6/bioc"

if ("BiocInstaller" %in% rownames(installed.packages()))
	remove.packages("BiocInstaller")

install.packages("BiocInstaller", repos=url)

builtins <- c("Matrix", "KernSmooth", "mgcv")

for (builtin in builtins)
    if (!suppressWarnings(require(builtin, character.only=TRUE)))
        BiocInstaller::biocLite(builtin)

suppressWarnings(BiocInstaller::biocValid(fix=TRUE, ask=FALSE))


library(BiocInstaller) # shouldn't be necessary


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

ap.db <- available.packages(contrib.url(biocinstallRepos()))
ap <- rownames(ap.db)

pkgs_to_install <- pkgs[pkgs %in% ap]

biocLite(pkgs_to_install)

# just in case there were warnings, we want to see them
# without having to scroll up:
warnings()

if (!is.null(warnings()))
{
    w <- capture.output(warnings())
    if (length(grep("is not available|had non-zero exit status", w)))
        quit("no", 1L)
}

suppressWarnings(BiocInstaller::biocValid(fix=TRUE, ask=FALSE))

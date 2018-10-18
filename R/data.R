#' Collection data for the 1001 genomes project accessions
#'
#' A dataset containing unique identifiers, names, locations, and collectors information for each of the 1135 accessions in the 1001 genomes project.
#'
#' @format A data frame with 1135 rows and 8 variables
#' \describe{
#'   \item{Ecotype.ID}{unique number for each accession which corresponds to
#'   the "Indiv" column in the VCF}
#'   \item{Name}{common name of the accession}
#'   \item{CS.Number}{germplasm identifier for ordering seeds from ABRC or
#'   NASC}
#'   \item{Country}{country in which the accession was initially collected}
#'   \item{Lat}{latitude of initial collection}
#'   \item{Long}{longitude of initial collection}
#'   \item{Collector}{name of collector}
#'   \item{Sequenced.by}{institution which sequenced the accession}
#' }
#'
#' @source <http://1001genomes.org/accessions.html>
"accessions"


#' Araport11 genome annotation
#'
#' @format A granges object with 789890 ranges across 5 nuclear chromosomes as
#' well as the mitochondrial and chloroplastic chromosomes.
"gr"

#' A list of all accessions Ecotype ID's
#'
#' @format vector
"strains"

#' Rank order of effects called by SNPeff
#'
#' @format A data frame with 19 rows and 2 variables
#' \describe{
#'   \item{strength}{the strength of the effect on function, where higher
#'   strengths are more likely to be deleterious}
#'   \item{effect}{the name of the effect}
#' }
"SNPeff_order"
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
#' @source \url<https://www.googleapis.com/fusiontables/v2/query?sql=SELECT%20id%20AS%20tg_ecotypeid,%20name,%20CS_number,%20country,%20latitude,%20longitude,%20collector,%20seq_by%20FROM%201oabIAuVSTkoG3qTAqM2sj-qiurysSWc5JEMq-TWm%20ORDER%20by%20id&key=AIzaSyB1GbyVSIOK12RJbFMkaIJjwhVNG-b8fjc&alt=csv>
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
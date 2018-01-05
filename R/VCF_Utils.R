# clean up object names to be consistent with bioconductor style (camelCaps)
# diversity_calc -> calcDiversity
# GT_freq -> GTfreq
# get_coding_diversity -> getCodingDiv
# plot_coding_diversity -> plotCodingDiv
# unique_coding_variants -> uniqueCodingVars
# add_ecotype_details -> addAccDetails
# Ecotype_column -> ecotypeColumn
# label_bySNPs -> labelBySNPs
# label_by_SNPs_kernel -> labelBySNPsKernel




### LIBRARIES ==================================================================

#library(VariantAnnotation)
#library(plyr)
#library(stringr)
#library(biomaRt)
#library(ggplot2)
#library(reshape2)
#library(vcfR)
#library(dplyr)
#library(ggmap)
#library(ggthemes)

### APP FUNCTIONS =============================================================


#' Run 1001 genomes browser shiny app
#'
#' @return Runs the app
#' @export
#'
#' @examples
run1001genomes <- function() {
  appDir <- system.file("shiny-app", package = "r1001genomes")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mypackage`.",
         call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}

### UTILITY FUNCTIONS =========================================================


#' Make a string containing the desired genomic ranges
#' of the format "chrom:start-stop"
#' @param data a data.frame with columns `chromosome_name`, `transcript_start`, and `transcript_end`
#'
#' @return a vector of strings with the format "chrom:start-stop"
#' @export
#'
#' @examples
#'
makeRegionString <- function (data) {
  # format "chrom:start-stop"
  return(paste(c(as.character(data$chromosome_name),
                 ":",
                 as.character(data$transcript_start),
                 "-",
                 as.character(data$transcript_end)), collapse=""))
}


#' Download and save a .VCF file from 1001genomes.org.
#'
#' @param fName the file name the downloaded vcf will be saved as eg.
#'    "data1.vcf"
#' @param strainStr A character string of comma separated strain or
#'    "Ecotype ID"s used in the construction of the URL
#' @param regionStr A character string of the format "[chrom]:[start]-[stop]"
#'    where [chrom], [start], and [stop] are numbers. this string is used
#'    directly in the construction of the URL
#' @param download A logical, if TRUE (default) the file will be downloaded, if
#'     FALSE, the file will not be downloaded, the URL will still be returned
#'
#' @return The URL used for the download
#' @export
#'
#' @examples
downloadData <- function (fName, strainStr, regionStr,
                          download=TRUE) {
  url <- c("http://tools.1001genomes.org/api/v1/vcfsubset/strains/", strainStr, "/regions/",
           regionStr, "/type/snpeff/format/vcf")
  url <- paste(url, collapse="")
  if (download == TRUE) {
    download.file(url, fName)
  }
  return (url)
}



#' Download two vcf files from 1001genomes.org and merge their contents
#'
#' @param fName a character string of the file name to save the final .vcf file to
#' @param strainVect numeric vector list of the strains
#' @param regionStr A character string of the format "[chrom]:[start]-[stop]"
#'    where [chrom], [start], and [stop] are numbers. this string is used
#'    directly in the construction of the URL
#'
#' @return a vcfR object of teh combined vcf file
#' @export
#'
#' @examples my.VCF <- downloadMerge("myFullVCF.vcf.gz", strains, "1:4368760-4371298")
downloadMerge <- function (fName, strainVect, regionStr) {
  strains <- as.character(strainVect)
  # the URL to download the VCF for all 1135 strains is too long so we have to split it into two sets

  splitPoint <- round(length(strains) / 2) #where to split the strain list

  # first file
  strainString <- paste(as.character(strains[1:splitPoint]), collapse=",")
    tempFile1 <- tempfile(fileext=".vcf")
  downloadData(tempFile1, strainString, regionStr)
  # second file
  strainString <- paste(as.character(strains[(splitPoint + 1):length(strains)]), collapse=",")
  tempFile2 <- tempfile(fileext=".vcf")
  downloadData(tempFile2, strainString, regionStr)
  #
  #load the two temporary vcf files then delete the temp files
  data1 <- vcfR::read.vcfR(tempFile1, verbose=FALSE, convertNA=TRUE)
  if (file.exists(tempFile1)) file.remove(tempFile1)

  data2 <- vcfR::read.vcfR(tempFile2, verbose=FALSE, convertNA=TRUE)
  if (file.exists(tempFile2)) file.remove(tempFile2)

  #combine the two gt fields, and write a new combined vcf file.
  data1@gt <- cbind(data1@gt, data2@gt[,-1])
  vcfR::write.vcf(data1, fName)

  return (data1)  #return the vcfR object of the combined .vcf
}



#' Parse the EFF field of the VCF files from 1001genomes.org
#'
#' @param tidyVCF tidyVCF$dat tibble to be parsed. should have an 'EFF' field
#' @param Transcript_ID Character string of the transcript ID. if NULL, the
#' function will look for it as the attribute "transcript_ID" of the dataframe
#'
#' @return new tidyVCF object with added columns for the parsed effect fields
#' @export
#' @import plyr
#'
#' @examples
parseEFF <- function (tidyVCF){
  EFFColNames = c("Effect", "Effect_Impact", "Functional_Class", "Codon_Change",
                  "Amino_Acid_Change", "Amino_Acid_Length", "Gene_Name", "Transcript_BioType",
                  "Gene_Coding", "Transcript_ID", "Exon_Rank", "Genotype_Number")
  data <- tidyVCF
  # stop if there is no "transcript_ID" attribute
  if (is.null(attr(tidyVCF, "transcript_ID"))) stop("Can not parse EFF field without transcript ID")
  transcript_ID <- attr(tidyVCF, "transcript_ID")
  output <- ddply(data, "POS", .fun=parseEFFKernel, transcript_ID, EFFColNames)
  attr(output, "transcript_ID") <- attr(tidyVCF, "transcript_ID")
  return (output)
}

#' Kernel of the parseEFF funciton, operates on a single row of data and parses the 'EFF' field
#'
#' @param data Single row of tidyVCF with "EFF" field
#' @param transcript_ID the transcript ID to be used to parse the 'EFF' field
#' @param EFFColNames the column names to be used for the effect fields of the 'EFF' column
#'
#' @return a data frame row consisting of the original data row with the parsed effect columns appended to it.
#' @export
#'
#' @examples
parseEFFKernel <- function (data, transcript_ID, EFFColNames){
  if (length(unique(data$EFF)) > 1) {
    warning("warning multiple effects found")
  }

  effect <- unique(data$EFF)
  # split by comma to generate a vector of different effects
  effect <- stringr::str_split(effect, pattern=",", simplify=TRUE)
  # split by "(", "|", and ")" to separate fields
  effect <- data.frame(stringr::str_split(effect, pattern="\\(|\\||\\)", simplify=TRUE), stringsAsFactors = FALSE)
  # remove last column of empty strings ""found after ")"
  effect <- effect[, 1:12, drop=FALSE]
  # add column names to effects
  colnames(effect) <- c(EFFColNames)
  # only keep effects that match the transcript ID
  effect <- effect[effect$Transcript_ID == transcript_ID, ]

  if (nrow(effect) > 0){   # if there are some effects remaining:
    #create a "gt_GT" column in the effect dataframe that matches the format of the VCF$dat
    effect$gt_GT <- paste(effect$Genotype_Number, "|", effect$Genotype_Number, sep = "")
    # merge the effect df with the original data df by the gt_GT field
    output <- merge(data, effect, by="gt_GT", all.x=TRUE)
  }

  else{  #if there are no effects matching the transcript ID, return the data unaltered
    output <- data
  }

  return(output)
}

#' download VCF, optionally in tidyVCF format
#'
#' @param geneInfo single row of geneInfo dataframe
#' @param strains numeric vector of strain/ecotypes to include
#' @param tidy logical, if true VCF will be provided in tidyVCF format, see
#' vcfR documentation
#' @param dataOnly logical, if true only return the $dat field of the VCF, (drop the metadata)
#'
#' @return the VCF, in the format as requested by the parameters of the function
#' @export
#'
#' @examples
VCFByTranscript <- function (geneInfo, strains, tidy=TRUE, dataOnly=TRUE){
  transcript_ID <- as.character(geneInfo$transcript_ID)
  regionString <- as.character(geneInfo$regionString)

  fName <- tempfile(fileext=".vcf.gz")
  VCF.out <- downloadMerge(fName, strains, regionString)

  if (tidy == TRUE){
    VCF.out <- vcfR::vcfR2tidy(VCF.out, single_frame = TRUE, info_fields = c("AC", "EFF"), format_fields = ("GT"))
    VCF.out$dat <- VCF.out$dat[!(is.na(VCF.out$dat$gt_GT)), ]
  }

  if (dataOnly){
    VCF.out <- VCF.out$dat
  }

  attr(VCF.out, "transcript_ID") <- transcript_ID # add transcript_ID attribute to the output object

  return (VCF.out)
}



#' download and store several VCFs in a list structure, named by their transcirpt ID
#'
#' @param geneInfo a geneInfo dataframe as created by the getGeneInfo() function
#' @param by
#'
#' @return list of VCF objects
#' @export
#' @import plyr
#'
#' @examples
VCFList <- function (geneInfo, by="transcript", tidy=TRUE) {
  output <- alply(geneInfo,.margins=1, .fun=VCFByTranscript, strains=strains)
  for (i in 1:length(output)){
    names(output)[i] <- attr(output[[i]], "transcript_ID")
  }
  return(output)
}




#' Get gene information
#'
#' @description Get the start and end position of each transcript of a set of
#'  genes using `bioMart` package to query the TAIR10 database via plants.ensembl.org
#'
#' @param genes a character vector of tair IDs of the genes to retrieve
#' @param firstOnly logical, if true only return transcript IDs containing ".1"
#' @param inputType
#' @param useCache logical, read from and write to a file of cached genes?
#'
#' @returna table containing fields from the Tair database on the provided genes
#' including "transcript_ID" and "regionString" columns required for other fuctions in this code
#' @export
#' @import plyr
#'
#' @examples
getGeneInfo <- function (genes, firstOnly=TRUE, inputType="tair_locus", useCache=TRUE) {

  retrievedInfo <- NULL
  genes2 <- genes
  cacheFile <- system.file("shiny-app", "geneInfoCache.txt", package="r1001genomes")

  if (useCache == TRUE){
    geneInfoCache <- read.table(file=cacheFile, header=TRUE, stringsAsFactors=FALSE)
    retrievedInfo <- geneInfoCache[geneInfoCache$tair_locus %in% genes, ]
    genes2 <- genes[!(genes %in% geneInfoCache$tair_locus)] #remove genes present in the cache from genes list
    print("new genes:")
    print(genes2)   # list new genes, not found in cache
  }

  output <- NULL
  if (length(genes2) > 0){
    tair10 <- useMart("plants_mart", host="plants.ensembl.org", dataset="athaliana_eg_gene")
    output <- getBM(attributes=c("tair_locus", "tair_symbol","ensembl_transcript_id", "chromosome_name", "start_position",
                                 "end_position", "strand", "transcript_start", "transcript_end"
                                  ), filters=inputType, values=genes2, mart=tair10)
    # create a list of strings encoding the chromosome and start and end position of all transcript IDs to be analyzed
    output$regionString <- as.character(alply(output, .fun=makeRegionString, .margins=1, .expand=FALSE))
    names(output)[names(output) == "ensembl_transcript_id"] <- "transcript_ID"
    output$transcript_length <- abs(output$transcript_end - output$transcript_start)

  }
  if (useCache == TRUE) {
    # append cache
    geneInfoCache <- unique(rbind(geneInfoCache, output))
    write.table(geneInfoCache, file=cacheFile, row.names=FALSE)
  }

  output <- rbind(retrievedInfo, output)

  if (firstOnly == TRUE) {
    # if firstOnly is TRUE, only return transcript IDs containing ".1"
    output <- output[(grepl(".1", output$transcript_ID, fixed=TRUE)), ]
  }
  return (output)
}

#' Calculate Nei's nucleotide diversity statistic for a single position,
#' given a vector of counts of uninque genotypes at that location
#'
#' @param GTfreq numeric vector listing counts of uninque genotypes
#'
#' @return numeric value nucleotide diversity
#' @export
#'
#' @examples
calcDiversity <- function (GTfreq){
  result <- (sum(GTfreq)**2 - sum(GTfreq**2))/(sum(GTfreq)**2)
}

#' Calculate nucleotide diversity for each position in the coding sequence
#'
#' @param tidyVCF the $dat field of a tidyVCF object
#'
#' @return input with Diversity field appended
#' @export
#'
#' @examples
Nucleotide_diversity <- function (tidyVCF){
  data <- unique(tidyVCF[, c("POS", "gt_GT", "Indiv")])
  GT_Frequencies <- plyr::count(data, c("POS", "gt_GT"))
  GT_Frequencies <- dplyr::group_by(GT_Frequencies, POS)
  diversityByPOS <- dplyr::summarise(GT_Frequencies, Diversity = calcDiversity(freq))
  output <- dplyr::full_join(tidyVCF, diversityByPOS, by="POS")

  if (!is.null(attr(tidyVCF, "transcript_ID"))) {
    # if the input had the "transcript_ID" attribute, pass it to the output
    attr(output, "transcript_ID") <- attr(tidyVCF, "transcript_ID")
  }

  return(output)
}

#' Add codon numbering to a row of data from a VCF object with parsed 'EFF' field
#' used in the nucleotide diversity kernel.
#' use with ddply()
#'
#' @param SNPData
#'
#' @return original data row with 'codon_Number' field appended
#' @export
#'
#' @examples
codonNumberKernel <- function (SNPData) {
  # add the codon number to a dataframe containing a "Amino_Acid_Change"
  changeStr <- SNPData$Amino_Acid_Change[grepl( "p.", SNPData$Amino_Acid_Change)][1]
  codonNumber <- stringr::str_extract_all(stringr::str_extract_all(changeStr, "p.[A-z]{3}[0-9]*")[[1]], "[0-9]+")[[1]]
  codonNumber <- as.numeric(codonNumber)
  rows <- nrow(SNPData)
  return (cbind(SNPData, "Codon_Number"=codonNumber))
}

#' Title
#'
#' @param geneInfo
#' @param strains
#'
#' @return
#' @export
#'
#' @examples
polymorphTable <- function (geneInfo, strains) {
  effects <- c("5_prime_UTR_variant",
               "intron_variant",
               "3_prime_UTR_variant",
               "synonymous_variant",
               "missense_variant",
               "upstream_gene_variant")
  tableData <- matrix(nrow=length(geneInfo$transcript_ID), ncol=length(effects) + 6)
  rownames(tableData) <- geneInfo$transcript_ID
  colnames(tableData) <- c(effects, "coding_total", "Pi_coding", "Pi_non_syn", "Pi_syn", "Pi_NS_Ratio", "Pi_transcript")
  tableData <- data.frame(tableData)

  #for each transcript
  for (i in 1:length(geneInfo$transcript_ID)) {
    tidyVCF <- VCFByTranscript(geneInfo[i, ], strains)
    data <- parseEFF(tidyVCF)
    data <- Nucleotide_diversity(data)

    #fill in the first part of the table
    variant_counts <- plyr::count(data, "Effect")
    for (j in 1:length(effects)){
      if (effects[j] %in% variant_counts$Effect){
        tableData[i,j] <- variant_counts[variant_counts$Effect %in% effects[j], "freq"]
      } else {
        tableData[i,j] <- 0
      }
    }
    tableData[i, "coding_total"] <- tableData[i, "missense_variant"] + tableData[i, "synonymous_variant"]

    #nucleotide diversity sums:
    reducedData <- unique(data[, c("POS", "Effect", "Diversity")])
    AA_Length <- unique(data$Amino_Acid_Length)
    AA_Length <- as.numeric(AA_Length[!is.na(AA_Length)])

    tableData[i, "Pi_non_syn"] <- sum(reducedData[reducedData$Effect %in% "missense_variant", "Diversity"]) / (3*AA_Length)
    tableData[i, "Pi_syn"] <- sum(reducedData[reducedData$Effect %in% "synonymous_variant", "Diversity"]) / (3*AA_Length)
    tableData[i, "Pi_NS_Ratio"] <- tableData[i, "Pi_non_syn"] / tableData[i, "Pi_syn"]

    tableData[i, "Pi_coding"] <- sum(unique(reducedData[reducedData$Effect %in% c("synonymous_variant","missense_variant") , c("POS", "Diversity")])$Diversity) / (3*AA_Length)
    tableData[i, "Pi_transcript"] <- sum(unique(reducedData[, c("POS", "Diversity")])$Diversity) / geneInfo[i, "transcript_length"]

  }
  return(tableData)

}




#' create a row of statistics for a given tidyVCF dataframe including non unique
#' variant counts and nucleotide diversity values.
#'
#' @param data tidyVCF
#' @param geneInfo geneInfo dataframe containing a row for the provided tidyVCF gene
#'
#' @return
#' @export
#'
#' @examples
polymorphRow <- function (data, geneInfo=NULL) {
  effects <- c("5_prime_UTR_variant",
               "intron_variant",
               "3_prime_UTR_variant",
               "synonymous_variant",
               "missense_variant",
               "upstream_gene_variant")

  # create table
  tableData <- data.frame(row.names=attr(data,"transcript_ID"))

  # polymorphism counts (non-unique)
  variant_counts <- plyr::count(data, "Effect")
  for (effect in effects){
    if (effect %in% variant_counts$Effect){
      tableData[effect] <- variant_counts[variant_counts$Effect %in% effect, "freq"]
    } else {
      tableData[effect] <- 0
    }
  }
  tableData$coding_total <- tableData$missense_variant + tableData$synonymous_variant

  #nucleotide diversity sums:
  reducedData <- unique(data[, c("POS", "Effect", "Diversity")])
  AA_Length <- unique(data$Amino_Acid_Length)
  AA_Length <- as.numeric(AA_Length[!is.na(AA_Length)])
  tableData$Pi_non_syn <- sum(reducedData[reducedData$Effect %in% "missense_variant", "Diversity"]) / (3*AA_Length)
  tableData$Pi_syn <- sum(reducedData[reducedData$Effect %in% "synonymous_variant", "Diversity"]) / (3*AA_Length)
  tableData$Pi_NS_Ratio <- tableData$Pi_non_syn / tableData$Pi_syn
  tableData$Pi_coding <- sum(unique(reducedData[reducedData$Effect %in% c("synonymous_variant","missense_variant"), c("POS", "Diversity")])$Diversity) / (3*AA_Length)

  if (!is.null(geneInfo)){
     transcript_length <- geneInfo$transcript_length[geneInfo$transcript_ID == attr(data,"transcript_ID")]
     tableData$Pi_transcript <- sum(unique(reducedData[, c("POS", "Diversity")])$Diversity) / transcript_length
  }

  return(tableData)
}



#' Title
#'
#' @param data tidyVCF with 'EFF' field parsed
#'
#' @return
#' @export
#'
#' @examples
getCodingDiv <- function(data){
  #input is tidy tibble/df with EFF field parsed and diversity calculated:
  # eg.
  # myvcf <- VCFByTranscript(geneInfo[1, ], strains)
  # mydata <- parseEFF(myvcf)
  # mydata <- Nucleotide_diversity(mydata)
  # coding_Diversity_Plot(mydata)

  coding_variants <- data[data$Effect %in% c("missense_variant", "synonymous_variant"), ]   #"stop_gained", "frameshift_variant"
  #extract uniuqe position and effect
  uniqueCodingVars <- unique(coding_variants[ , c("POS", "Effect",
                                                        "Amino_Acid_Change",
                                                        "Diversity") ])
  #add codon number to uniqueCodingVars
  uniqueCodingVars <-plyr::ddply(uniqueCodingVars, .fun=codonNumberKernel,
                                 .variables=c("POS", "Amino_Acid_Change"))
  return(uniqueCodingVars)
}

#' Title
#'
#' @param uniqueCodingVars
#'
#' @return
#' @export
#' @import ggplot2
#' @import ggthemes
#'
#' @examples
plotCodingDiv <- function(uniqueCodingVars){
  #plot the diversity
  plot <- ggplot(uniqueCodingVars, aes(x=Codon_Number,y=Diversity, colour=Effect, shape = Effect)) +
    geom_point(size = 4) +
    scale_y_log10(breaks=c(0.0001, 0.001, 0.01, 0.1),limits=c(0.0001, 1)) +
    #scale_colour_manual(values=c(synonymous_diversity="blue", missense_diversity="red")) +
    ylab("nucleotide diversity, log scale") + theme_few(base_size = 18) +
    scale_color_colorblind()
  return(plot)
}


#' Add accession metadata to a dataset containing ecotype IDs
#'
#' @param data
#' @param ecotypeColumn
#'
#' @return
#' @export
#'
#' @examples
addAccDetails <- function(data, ecotypeColumn="Indiv") {
  # add ecotype details (location, collector, sequencer) to any df containing an "Indiv" column
  ecoIDs <- r1001genomes::accessions
  return(merge(data, ecoIDs, by.x=ecotypeColumn, by.y="Ecotype.ID", all.y=TRUE))
}


# buildGT <- function(indivData) {
#   # use with ddply, chunk by "Indiv"
#   indivGT <- unique(indivData[,c("POS", "gt_GT")])
#   Indiv <- indivData[1,"Indiv"]
#   rownames(indivGT) <- indivGT[,1]
#   indivGT <- t(indivGT[,2, drop=FALSE])
#   rownames(indivGT) <- Indiv
#   return(indivGT)
# }

#' Label accessions with the variants they contain
#'
#' @param data a data frame of variants with a column named "Indiv"
#' @param collapse should each accession will be a single line
#'
#' @return a df with a single row per accession, with a new column "SNPs" that
#'  has a single text string detailing all the variants from data within that
#'  accession
#' @export
#' @import plyr
#'
#' @examples
labelBySNPs <- function(data, collapse=TRUE) {
  # creates a df with a single row per individual, with a new column "SNPs" that
  # has a single text string detailing the

  output <- ddply(data, .variables="Indiv", .fun=labelBySNPsKernel, collapse=collapse)

  output <- addAccDetails(output)
  return(output)
}



#' Title
#'
#' @param indivData
#' @param collapse
#'
#' @return
#' @export
#'
#' @examples
labelBySNPsKernel <- function(indivData, collapse=TRUE) {
  #if collapse == TRUE, each accession will be a single line.

  # store ecotypeID as a single value
  Indiv <- indivData[1,"Indiv"]

  # filter only rows with an effect (ie not reference or NA)
  data <- indivData[!is.na(indivData$Effect), ]
  # order the rows by transcript ID frist, then Amino_Acid_change field
  # note: Amino_Acid_change should be always present even on non coding UTRs and introns
  data <- data[order(data[, "Transcript_ID"], data[, "Amino_Acid_Change"]), ]

  if (collapse == TRUE) {
    collapseString <- ","
  } else {
    collapseString <- NULL
  }

  SNPstring <- paste("[", data[, "Transcript_ID"],"|", data[, "Amino_Acid_Change"], "]", collapse=collapseString)

  output <- data.frame(Indiv, SNPs=SNPstring, stringsAsFactors=FALSE)
  return(output)

}




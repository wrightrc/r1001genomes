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
#' @param tidyVCF tidyVCF$dat tibble to be parsed. should have an "EFF" field in
#' its $dat section
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

#' Title
#'
#' @param data
#' @param transcript_ID
#' @param EFFColNames
#'
#' @return
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
#'
#' @return
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

  attr(VCF.out, "transcript_ID") <- transcript_ID

  return (VCF.out)
}



#' download and store several VCFs in a list structure, named by their transcirpt ID
#'
#' @param geneInfo
#' @param by
#'
#' @return
#' @export
#' @import plyr
#'
#' @examples
VCFList <- function (geneInfo, by="transcript", tidy=TRUE) {
  output <- alply(geneInfo,.margins=1, .fun=VCFByTranscript, strains=strains)
  #names(output) <- attr(output, "split_labels")$transcript_ID
  for (i in 1:length(output)){
    #attr(output[[i]], "transcript_ID") <- names(output)[i]
    names(output)[i] <- attr(output[[i]], "transcript_ID")
  }
  return(output)
}




#' Get gene information
#'
#' @description Get the start and end position of each transcript of a set of
#'  genes using `bioMart` package to query the TAIR10 database via plants.ensembl.org
#'
#' @param genes
#' @param firstOnly
#' @param inputType
#' @param useCache
#'
#' @return
#' @export
#' @import plyr
#' @import biomaRt
#'
#' @examples
getGeneInfo <- function (genes, firstOnly=TRUE, inputType="tair_locus", useCache=TRUE) {
  # input "genes" should be a character vector of tair IDs
  # use bioMart to find the start and end position and transcript IDs for a given set of genes
  # outputs a table containing "transcript_ID" and "regionString" columns required for other fuctions in this code

  retrievedInfo <- NULL
  genes2 <- genes
  if (useCache == TRUE){
    geneInfoCache <- read.table(file=system.file("shiny-app",
                                                 "geneInfoCache.txt",
                                                 package = "r1001genomes"),
                                header=TRUE, stringsAsFactors=FALSE)
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
    write.table(geneInfoCache, file="geneInfoCache.txt", row.names=FALSE)
  }


  output <- rbind(retrievedInfo, output)

  if (firstOnly == TRUE) {
    # if firstOnly is TRUE, only return transcript IDs containing ".1"
    output <- output[(grepl(".1", output$transcript_ID, fixed=TRUE)), ]
  }
  return (output)
}

#' Title
#'
#' @param GTfreq
#'
#' @return
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
#' @return
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

#' Add codon numbering to the nucleotide diversity kernel
#'
#' @param SNPData
#'
#' @return
#' @export
#'
#' @examples
codonNumberKernel <- function (SNPData) {
  # add the codon number to a dataframe containing a "Amino_Acid_Change"
  # use with ddply()
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




#' Title
#'
#' @param data
#' @param geneInfo
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

  variant_counts <- plyr::count(data, "Effect")

  tableData <- data.frame(row.names=attr(data,"transcript_ID"))

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

  tableData$Pi_coding <- sum(unique(reducedData[reducedData$Effect %in% c("synonymous_variant","missense_variant") , c("POS", "Diversity")])$Diversity) / (3*AA_Length)

  if (!is.null(geneInfo)){
     transcript_length <- geneInfo$transcript_length[geneInfo$transcript_ID == attr(data,"transcript_ID")]
     tableData$Pi_transcript <- sum(unique(reducedData[, c("POS", "Diversity")])$Diversity) / transcript_length
  }

  return(tableData)

}



#' Title
#'
#' @param data
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
                                                  "Transcript_ID",
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


#' Make an alignment of Coding sequences
#'
#' @param IDs a vector of TAIR transcript IDs, e.g. "AT3G26890.1"
#'
#' @return aligned CDS and amino acid sequences as a list of XStringSet objects
#' @export
#' @import GenomicFeatures
#' @import S4Vectors
#' @importFrom XVector "subseq"
#'
#' @examples
#' IDs <- c("AT3G62980.1", "AT3G26810.1")
#' alignment <- alignCDS(IDs)
#' browseSeqs(alignment[[2]])
alignCDS <- function(IDs) {
  #use_package("BSgenome.Athaliana.TAIR.TAIR9", "imports")
  Athaliana <- BSgenome.Athaliana.TAIR.TAIR9::BSgenome.Athaliana.TAIR.TAIR9
  #use_package("GenomicFeatures", "imports")
  #use_package("rtracklayer", "imports")
  # gr is now built in
  # gr <- rtracklayer::import(system.file("extdata",
  #           "Araport11_GFF3_genes_transposons.201606.gff.gz",
  #           package = "r1001genomes"))
  gr_sub <-   gr[which(grepl(pattern = paste(
    gsub(pattern = "\\..*", replacement = "", x = IDs), collapse = "|"),
                             x = mcols(gr)$ID)),]
  txdb <- GenomicFeatures::makeTxDbFromGRanges(gr_sub) # maybe use exclude.stop
  CDSseqs <- GenomicFeatures::extractTranscriptSeqs(Athaliana,
                GenomicFeatures::cdsBy(txdb, by = 'tx', use.names = TRUE))
  #devtools::use_package("XVector", "imports")
  CDSseqs.xstop <- XVector::subseq(CDSseqs, start = rep(1,length(CDSseqs)),
                          end = nchar(CDSseqs)-3)
  CDSseqs.xstop
  #devtools::use_package("DECIPHER", "depends")
  # Update to imports once DECIPHER has fixed environment issue
  CDSAlignment <- DECIPHER::AlignTranslation(CDSseqs.xstop, type = "both")
  return(CDSAlignment)
}

#' Make an alignment data frame for plotting
#'
#' @param alignment an AAStringSet or DNAStringSet resulting from an alignment made with DECIPHER.
#'
#' @return a tidy long data frame of the alignment with columns `seq_name`, `letter`, `aln_pos`, and `seq_pos`, where `aln_pos` is the position of the letter in the alignment of all sequences and `seq_pos` is the position of the letter in the native sequence.
#' @importFrom magrittr "%>%"
#' @import DECIPHER
#' @export
#'
#' @examples
#' IDs <- c("AT3G62980.1", "AT3G26810.1")
#' alignment <- alignCDS(IDs)
#' makeAlnDF(alignment[[2]])
#'
makeAlnDF <- function(alignment){
  seqs <- strsplit(as.character(alignment), split = NULL)
  #devtools::use_package("plyr")
  aln_df <- plyr::ldply(.data = seqs, .id = "seq_name", .fun = function(sequence){
    data.frame(
      "letter" = sequence,
      "aln_pos" = 1:length(sequence)
    )
  })
  #devtools::use_package("magrittr")
  #devtools::use_package("dplyr")
  aln_df <- aln_df %>% dplyr::group_by(seq_name) %>% dplyr::mutate(
    seq_pos = {
      matches <- which(letter %in% c(LETTERS, letters, "*"))
      gaps <- which(letter %in% c("-"))
      seqPos <- vector(length = length(letter))
      seqPos[gaps] <- "-"
      seqPos[matches] <- 1:length(matches)
      return(seqPos)
    })
  return(aln_df)
}


#' Add SNPs and indels to an alignment data frame
#'
#' @param aln_df an alignment data frame resulting from
#'  \link[r1001genomes]{makeAlnDF}
#' @param SNPs a data frame of SNPs
#' @param by_aln_SNPs a named list of character objects with each equivalency
#' representing matching columns in `aln_df` (on the LHS) and `SNPs`
#' (on the RHS), e.g. `"seq_name" = "transcript_id"`
#'
#' @return aln_df with the addition of a `variants` column containing a string
#'  listing the variant types at each position
#'
#' @importFrom magrittr "%>%"
#' @export
#'
#' @examples
#' # make an alignment
#' IDs <- c("AT3G62980.1", "AT3G26810.1")
#' alignment <- alignCDS(IDs)
#' # make an alignment data frame
#' aln_df <- makeAlnDF(alignment[[2]])
#' # load a VCF list
#' vcf <- readRDS(file = system.file("shiny-app", "demo_VCFs.rds",
#'    package = "r1001genomes"))
#' vcf <- ldply(.data = vcf, .fun = subset,
#'   !is.na(Transcript_ID) & gt_GT != "0|0")
#' coding_vcf <- getCodingDiv(vcf)
#' addSNPsToAlnDF(aln_df, coding_vcf, seq_name = Transcript_ID,
#' seq_pos = Codon_Number)
#'
addSNPsToAlnDF <- function(aln_df, SNPs, seq_name = Transcript_ID,
                           seq_pos = Codon_Number, effect = Effect){
  seq_name <- dplyr::enquo(seq_name)
  seq_pos <- dplyr::enquo(seq_pos)
  effect <- dplyr::enquo(effect)
  temp <- SNPs %>%
    dplyr::group_by(!!seq_name, !!seq_pos) %>%
    dplyr::summarise(variants = {switch(as.character(length(unique(!!effect))),
                                 "0" = NA,
                                 "1" = unique(!!effect),
                                 paste(sort(unique(!!effect)),
                                       collapse = " & "))})
  temp[,quo_text(seq_pos)] <- as.character(temp[,quo_text(seq_pos)])
  aln_df <- dplyr::left_join(x = aln_df, y = temp,
                        by = c("seq_name" = quo_text(seq_name),
                               "seq_pos" = quo_text(seq_pos)))
  return(aln_plot)
}

#' Add alignment positions to an annotation data frame
#'
#' @param anno_df
#' @param aln_df
#'
#' @return
#' @export
#'
#' @examples
addAlnPosToAnno <- function(anno_df, aln_df){

}

#' String alignment geom for plotting XStringSet Alignments
#'
#' @param data an alignment data.frame made from an XStringSet. In the future,
#'  we plan to add pretty defaults allowing XStringSets to be plotted.
#' @param ... additional parameters passed on to geoms like aes()
#'
#' @return
#' @export
#'
#' @examples
geom_str_align <- function(mapping = NULL, data = NULL){
  ggplot(aln_plot, aes(x, as.integer(seqname), group = seqPos)) +
    geom_rect(data = na.omit(aln_plot), aes(xmin = x - 0.5, xmax = x + 0.5,
                                            ymin = as.integer(seqname) - 0.5,
                                            ymax = as.integer(seqname) + 0.5,
                                            fill = variants), alpha = 0.8) +
    geom_text(aes(label=letter), alpha= 1,
              check_overlap = TRUE) +
    scale_x_continuous(breaks=seq(1,max(aln_plot$x), by = 10)) +
    scale_y_continuous(breaks = c(1:length(levels(aln_plot$seqname))),
                       labels = levels(aln_plot$seqname)) +
    # expand increases distance from axis
    xlab("") +
    ylab("") +
    #scale_size_manual(values=c(5, 6)) + # does nothing unless 'size' is mapped
    theme_logo(base_family = "Courier") +
    theme(panel.grid = element_blank(), panel.grid.minor = element_blank())
}

#' make DNAStrings of sequences for each gene of each accession
#'
#' @param SNPs a VRanges object containing the SNPs you would like to make
#'  sequences for
#' @param genome a BSgenomes object corresponding to the genome used to call
#'  the SNPs and also corresponding to the ranges set in 'ranges'
#' @param ranges a named GRanges object containing the ranges of the sequence
#'  elements you would like to extract
#'
#' @return
#' @export
#'
#' @examples
promoterVariantToString<-function(SNPs, genome, ranges, files_out = TRUE){
  #get reference sequences
  #mcols(ranges)$seqs <- getSeq(x = genome, ranges)
  #ranges(SNPs) <- ranges(mapToTranscripts(SNPs, ranges))
  #mcols(SNPs)$GENEID <- seqnames(mapToTranscripts(SNPs, ranges))
  #loop to generate accession sequences for each seq of each gene
  llply(mcols(ranges)$names,function(gene) {
    #pull all coding snps mapped to gene
    genes<-subset(SNPs, mcols(SNPs)$GENEID == gsub("\\..*$","",gene))
    #for each unique accession in the new SNP_list
    geneList<-llply(unique(sampleNames(SNPs)),function(acc) {
      #get variants in the strain
      variants <- unique(subset(genes, sampleNames(genes) == acc))
      #create sequence
      strainVar <- replaceAt(x=mcols(ranges)$seqs[
        which(mcols(ranges)$names==gene)],
        at = ranges(variants), value = alt(variants))
      #name it
      names(strainVar)<-paste(acc,names(strainVar),sep='-')
      strainVar
    })
    #stitch this list together for fasta creation
    geneSet<-do.call(c,geneList)
    geneSet<-c(mcols(ranges)$seqs[
      which(mcols(ranges)$names==gene)],geneSet)
    if(files_out) writeXStringSet(geneSet,
                                  paste0(gsub("\\..*$","",gene),"seqs.fasta"))
    geneSet
  })
}

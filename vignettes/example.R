## ---- echo = FALSE, message=FALSE, warning=FALSE-------------------------
knitr::opts_chunk$set(warning = FALSE, collapse = TRUE, comment = "#>", tidy = TRUE, out.width = "90%", fig.height=4, fig.width=7, out.extra='style="margin: auto; display: block; padding-top: 15px;"')


## ---- echo = FALSE-------------------------------------------------------
knitr::include_graphics("Screenshots/App_Initial.PNG")
knitr::asis_output("
                   ") # This newline required for headings to render after image

## ---- eval=FALSE---------------------------------------------------------
#  run1001genomes()

## ---- echo=FALSE---------------------------------------------------------
knitr::include_graphics("Screenshots/Tab_1_Select_Genes.PNG")
knitr::asis_output("
                   ") # This newline required for headings to render after image

## ------------------------------------------------------------------------
library(r1001genomes)
library(DT)
library(dplyr)
library(biomaRt)

## ------------------------------------------------------------------------

tair_ids <- c("AT3G62980", "AT3G26810")
geneInfo <- getGeneInfo(genes = tair_ids, firstOnly = TRUE, inputType = "tair_locus", useCache = TRUE)

## ------------------------------------------------------------------------
datatable(geneInfo[, -c(5,6,9)], colnames = c("tair locus", "symbol", "transcript", 
                                              "Chr", "transcript \nstart", "transcript \nend", 
                                              "transcript \nlength"), rownames = FALSE, 
          options=list(paging=FALSE, searching=FALSE))

## ------------------------------------------------------------------------
AFB_loci <- read.csv(system.file("extdata", "AFB_gene_ids.csv", package = "r1001genomes"))
paste(AFB_loci$tair_locus, collapse = ",")

## ---- eval=FALSE---------------------------------------------------------
#  YFG_VCF <- VCFList(geneInfo = geneInfo, by = "transcript", tidy = FALSE)

## ----eval = FALSE--------------------------------------------------------
#  YFG_VCF <- llply(YFG_VCF, parseEFF)
#  YFG_VCF <- llply(YFG_VCF, Nucleotide_diversity)

## ------------------------------------------------------------------------
formatRound(datatable(TIR1_AFB2_nonunique_variants,
                  colnames = c("transcript", "5' UTR", "intron", "3' UTR",
                               "coding \n synonymous", "coding \n missense",
                               "coding \n total"), 
                  escape = FALSE, rownames = FALSE,
    options=list(paging=FALSE, searching=FALSE)), columns = 2:7, digits = 0)

## ------------------------------------------------------------------------
TIR1_AFB2_unique_variants <- ldply(YFG_VCF, variantCounts, unique = TRUE, .id = "Transcript_ID")

formatRound(datatable(TIR1_AFB2_unique_variants,
                  colnames = c("transcript", "5' UTR", "intron", "3' UTR",
                               "coding \n synonymous", "coding \n missense",
                               "coding \n total"), 
                  escape = FALSE, rownames = FALSE,
    options=list(paging=FALSE, searching=FALSE)), columns = 2:7, digits = 0)

## ------------------------------------------------------------------------
TIR1_AFB2_diversity_table <- ldply(YFG_VCF, diversityStats, geneInfo = geneInfo, .id = "Transcript_ID")

formatRound(datatable(TIR1_AFB2_diversity_table,
                  colnames = c("transcript",
                               "&pi; transcript",
                               "&pi; coding",
                               "&pi;<sub>N</sub>",
                               "&pi;<sub>S</sub>",
                               "&pi;<sub>N</sub>/&pi;<sub>S</sub>"), 
                  escape = FALSE,
                  options = list(paging=FALSE, searching=FALSE)),
                columns = 1:6, digits = 6)

## ---- echo=FALSE---------------------------------------------------------
imagePath <- "Screenshots/Tab_2_Select_Gene.PNG"
knitr::include_graphics(imagePath)
knitr::asis_output("
                   ") # This newline required for headings to render after image

## ------------------------------------------------------------------------
YFG_Coding_Diversity <- llply(YFG_VCF, getCodingDiv)
lapply(YFG_Coding_Diversity, plotCodingDiv)


## ------------------------------------------------------------------------
AFB2_missense <- YFG_VCF$AT3G26810.1 %>% 
  dplyr::filter(Diversity > 3*10^(-3) & gt_GT != "0|0" & 
                  Effect == "missense_variant")
head(AFB2_missense)

## ------------------------------------------------------------------------
AFB2_missense <- labelBySNPs(AFB2_missense)
head(AFB2_missense)

## ------------------------------------------------------------------------
#load the leaflet package 
library(leaflet)
AFB2_missense <- rbind(AFB2_missense[is.na(AFB2_missense$SNPs), ],
                       AFB2_missense[!is.na(AFB2_missense$SNPs), ])

    # make a field with text to be displayed when clicking on a marker
    AFB2_missense$popup <- paste("EcoID:",  AFB2_missense$Indiv,"Name:", AFB2_missense$Name, " SNPs:", AFB2_missense$SNPs)

    # create the color pallet for the map points
    pal <- RColorBrewer::brewer.pal(8, "Set1")
    pallet <- colorFactor(palette=pal, domain=AFB2_missense$SNPs)
map <- leaflet()
    map <- addProviderTiles(map, providers$Stamen.TonerLite,
                     options = providerTileOptions(noWrap = TRUE))

    # groupnames to be used by draw groups of points as separate layers below
    groupnames <- unique(AFB2_missense$SNPs)
    groupnames <- groupnames[!is.na(groupnames)]

    # add markers for NA points first so they are furthest back layer
    map <- addCircleMarkers(map, data=AFB2_missense[is.na(AFB2_missense$SNPs), ], color= "#9b9b9b", group="NA",
                            radius=6, stroke=FALSE, fillOpacity=0.6)

    # for each of the group names, add a set of markers
    for (SNP in groupnames){
          map <- addCircleMarkers(map, data=AFB2_missense[AFB2_missense$SNPs == SNP, ], color= ~pallet(SNPs), group= SNP,
                            radius=6, popup= ~popup, stroke=FALSE, fillOpacity=0.85)
    }

    # add the legend to the map
    map <- addLegend(map, position="bottomright", pal=pallet,
                     values=AFB2_missense$SNPs, title="Marker Colors", opacity=1)

    # add layer control to map to turn on or off groups of points
    map <- addLayersControl(map, overlayGroups=c(groupnames, "NA"),
                            options = layersControlOptions(collapsed = TRUE),
                            position="bottomleft")

    map

## ------------------------------------------------------------------------
library(ggplot2)
library(ggmap)
library(ggrepel)

europe <- get_map(location = 'europe', zoom = 4)
michigan <- get_map(location = 'michigan', zoom = 5)


ggmap(europe) + xlab('longitude') + ylab('latitude') + geom_point(data = AFB2_missense, aes(x = Long, y = Lat), alpha = 0.5, shape = 18, color = "darkgrey") + geom_label_repel(data = subset(AFB2_missense, grepl(pattern = "204", AFB2_missense$SNPs)), mapping = aes(label = Name, x = Long, y = Lat), fill = "white", point.padding = unit(0.5, "lines"), segment.size = 1) 

ggmap(michigan) + xlab('longitude') + ylab('latitude') + geom_point(data = AFB2_missense, aes(x = Long, y = Lat), alpha = 0.5, shape = 18, color = "darkgrey") + geom_label_repel(data = subset(AFB2_missense, grepl(pattern = "516", AFB2_missense$SNPs)), mapping = aes(label = Name, x = Long, y = Lat), point.padding = unit(0.5, "lines"), segment.size = 1) 



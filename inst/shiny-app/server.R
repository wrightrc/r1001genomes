
# Server =================================================================

library(shiny)
library(biomaRt)
library(leaflet)
library(RColorBrewer)
library(r1001genomes)
library(knitr)
library(stringr)
library(DECIPHER)
library(ggseqlogo)
library(shinyBS)
library(ggplot2)
library(ggpmisc)
library(dplyr)
library(cowplot)


parseInput <- function (textIn) {
  names <- str_extract_all(textIn, "AT[1-5]G[0-9]{5}")
  return (names[[1]])
}

parseFilterText <- function (textIn) {
  inputList <- strsplit(textIn, ", ")
  inputList <- gsub(" ", "", inputList[[1]]) # remove extra spaces
  return(inputList)
}


enableBookmarking(store = "url")

server <- function(input, output, session){

  ##
  ##  ------------------------------------------------------------------------
  ##
  ##  INPUT AREA #############


#### tab1.buttons ####
  tab1.buttons <- reactiveValues(last_button="none pressed", total_presses=0)
  observeEvent(input$STATS_submit,{
    if (input$STATS_submit > 0){
      tab1.buttons$last_button <- "STATS_submit"
      tab1.buttons$total_presses <- tab1.buttons$total_presses + 1
    }
  })
  observeEvent(input$file_submit,{
    if (input$file_submit > 0){
      tab1.buttons$last_button <- "file_submit"
      tab1.buttons$total_presses <- tab1.buttons$total_presses + 1
    }
  })
#### all.Genes ####
  all.Genes <- eventReactive({tab1.buttons$total_presses},{
    req(tab1.buttons$last_button!="none pressed")
    if (tab1.buttons$last_button == "file_submit"){
      genes <- geneInfoFromFile(input$genesFile$datapath)
      req(genes != FALSE)
      return(genes)
    }
    if (input$STATS_quick_demo){
      names <- c("AT3G62980", "AT3G26810")
      genes <- getGeneInfo(names)
      req(genes != FALSE)
      return(genes)
    }
    # list of genes for tab 1, updated on pressing submit button
    names <- parseInput(input$gene_ids)
    genes <- getGeneInfo(names)
    req(genes != FALSE)
    return(genes)
  })
#### all.GeneChoices ####
  all.GeneChoices <- reactive({
    # displayNames <- paste(all.Genes()$transcript_ID, " (", all.Genes()$tair_symbol, ")", sep="" )
    # displayNames <- gsub(" \\(\\)", displayNames, replacement="")  # if no tair symbol, remove empty parens.
    displayNames <- paste(all.Genes()$tair_symbol, " (", all.Genes()$transcript_ID, ")", sep="")
    output <- all.Genes()$transcript_ID
    names(output) <- displayNames
    return(output)
  })
#### anno_df ####
  anno_df <- eventReactive(input$annoSubmit,{
    anno_df <- readAnnotationFile(input$annoFile$datapath, gene_info = all.Genes())
    return(anno_df)
  })
#### annoTemplateDownload ####
  output$annoTemplateDownload <- downloadHandler(
    filename="annotations_template.csv",
    content = function(file) {
      file.copy("annotations_template.csv", file)
    }
  )
#### all.VCFList ####
  all.VCFList <- reactive({
    if(isolate(input$STATS_quick_demo) & (tab1.buttons$last_button == "STATS_submit")) {
      all.Genes() # DO NOT DELETE this is here to make all.VCFList update after unchecking quickdemo
      return(readRDS(file = system.file("shiny-app", "demo_VCFs.rds",
                                        package = "r1001genomes")))
    }
    withProgress(message="downloading data from 1001genomes.org",
                   detail="this will take a while, progress bar will not move",
                   value=0.3, {
                     output <- VCFList(all.Genes())
                     setProgress(value=0.7, message="downloading complete, processing data...",
                                 detail="Parsing EFF field")
                     output <- llply(output, parseEFF)
                     setProgress(value=0.9, message=NULL,
                                 detail="Calculating nucleotide diversity")
                     output <- llply(output, Nucleotide_diversity)
                     output <- llply(output, addAccDetails)
                     setProgress(value=1)
    })
    return(output)
  })

  ##   _________
  ##  /  tab1   \
  ##             --------------------------------------------------
  ## Tab 1 ####################

#### tab1.genes_table ####
  output$tab1.genes_table <- DT::renderDataTable(DT::datatable(all.Genes()[, -c(5,6,7,10)], colnames = c("tair locus", "symbol", "transcript", "Chr", "transcript \nstart", "transcript \nend", "transcript \nlength"), rownames = FALSE, options=list(paging=FALSE, searching=FALSE)))
#### tab1.nonUniqueVariants ####
  tab1.nonUniqueVariants <- eventReactive({all.VCFList()},{
    req(isolate(tab1.buttons$last_button)!="none pressed")
    ldply(all.VCFList(), variantCounts, unique=FALSE, .id="transcript_ID")
  })
#### tab1.uniqueVariants ####
  tab1.uniqueVariants <- eventReactive({all.VCFList()},{
    req(isolate(tab1.buttons$last_button)!="none pressed")
    ldply(all.VCFList(), variantCounts, unique=TRUE, .id="transcript_ID")
  })
#### tab1.divStats ####
  tab1.divStats <- eventReactive({all.VCFList()},{
    req(isolate(tab1.buttons$last_button)!="none pressed")
    ldply(all.VCFList(), diversityStats, geneInfo=isolate(all.Genes()), .id="transcript_ID")
  })
#### SNPStats ####
  SNPStats <- reactive({
    req(isolate(tab1.buttons$last_button)!="none pressed")
    # rename column names on unique variant counts.
    uniqueVariantsRenamed <- tab1.uniqueVariants()
    colnames(uniqueVariantsRenamed) <- paste(colnames(uniqueVariantsRenamed),
                                             "unique", sep="_")
    cbind(tab1.nonUniqueVariants(), uniqueVariantsRenamed[, -1], tab1.divStats()[, -1])
  })
#### tab1.SNPcounts ####
  output$tab1.SNPcounts <- DT::renderDataTable({
    table <- tab1.nonUniqueVariants()
    colnames(table) <- c("transcript", "symbol", "5' UTR", "intron", "3' UTR",
                         "coding \n synonymous", "coding \n missense",
                         "stop\ngained", "frameshift\nvariant",
                         "upstream", "coding \n total")
    table <-  table[,c(TRUE,TRUE, colSums(table[,3:11])!=0)] # remove columns with all zeros
    table <- DT::datatable(table,rownames = FALSE, options=list(paging=FALSE, searching=FALSE))
    return(table)
  })
#### tab1.SNPcountsUnique ####
  output$tab1.SNPcountsUnique <- DT::renderDataTable({
    table <- tab1.uniqueVariants()
    colnames(table) <- c("transcript", "symbol", "5' UTR", "intron", "3' UTR",
                         "coding \n synonymous", "coding \n missense",
                         "stop\ngained", "frameshift\nvariant",
                         "upstream", "coding \n total")
    table <-  table[,c(TRUE,TRUE, colSums(table[,3:11])!=0)] # remove columns with all zeros
    table <- DT::datatable(table,rownames = FALSE, options=list(paging=FALSE, searching=FALSE))
    return(table)
  })
#### tab1.Diversity_table ####
  output$tab1.Diversity_table <- DT::renderDataTable(
    DT::formatRound(DT::datatable(tab1.divStats(),
                  #
                  colnames = c("transcript",
                               "symbol",
                               "&pi;<sub>N</sub>",
                               "&pi;<sub>S</sub>",
                               "&pi;<sub>N</sub>/&pi;<sub>S</sub>",
                               "&pi; coding",
                               "&pi; transcript"),
                  rownames = FALSE, escape = FALSE,
                  options = list(paging=FALSE, searching=FALSE)),
                columns = 2:7, digits = 6))

  output$tab1.downloadStats <- downloadHandler(
    filename=function(){
      paste("SNPStats-", Sys.time(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(SNPStats(), file, row.names=FALSE)
    }
  )
#### tab1.downloadGeneInfo ####
  output$tab1.downloadGeneInfo <- downloadHandler(
    filename=function(){
      paste("GeneInfo-", Sys.time(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(all.Genes(), file, row.names=FALSE)
    }
  )
  ##                 _________
  ##                /  tab2   \
  ## ---------------           -------------------------------------
  ## Tab 2 ###################
#### tab2.selectGene ####
  output$tab2.selectGene <- renderUI({
    tagList(
      selectInput("tab2.transcript_ID", label=NULL, choices=all.GeneChoices()),
      actionButton(inputId="tab2.Submit", label = "Submit")
    )
  })
#### tab2.Genes ####
  tab2.Genes <- eventReactive(input$tab2.Submit, {
      #gene Info for gene on tab 2, updates on 'submit' button press
    # names <- parseInput(input$plotGene)
    # genes <- getGeneInfo(names[1])
    # return(genes)
    return(all.Genes()[ all.Genes()$transcript_ID == input$tab2.transcript_ID,])
  })
#### tab2.gene_table ####
  output$tab2.gene_table <- DT::renderDataTable(DT::datatable(tab2.Genes()[, -c(5,6,7,10)], colnames = c("tair locus", "symbol", "transcript", "Chr", "transcript \nstart", "transcript \nend", "transcript \nlength"), rownames = FALSE, options=list(paging=FALSE, searching=FALSE)))
    #rendered table of Gene info
#### tab2.tableData ####
  #tab2.tableData <- reactive({load_tab_2_Data(tab2.Genes())})
    #SNP reactive data
  tab2.tableData <- eventReactive(input$tab2.Submit, {
    tab2data <- all.VCFList()[[input$tab2.transcript_ID]]
    coding_variants <- getCodingDiv(tab2data)
    return(coding_variants)
  })
#### Diversity_Table ####
  output$Diversity_Table <- DT::renderDataTable(tab2.tableData())
    #render table of diversity data

  output$tab2.downloadSNPData <- downloadHandler(
    filename=function(){
      paste("SNPData-", Sys.time(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(tab2.tableData(), file, row.names=FALSE)
    }
  )

#### diversityPlot ####
  output$diversityPlot <- renderPlot({
    p <- plotCodingDiv(uniqueCodingVars = tab2.tableData())
    if(!is.null(input$annoFile)){
    p <- append_layers(p,list(
      geom_rect(data = subset(anno_df()$domains,
                              transcript_ID == input$tab2.transcript_ID),
                mapping = aes(xmin = as.integer(start),
                              xmax = as.integer(end),
                              fill = annotation),
                ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha = 0.2),
        geom_rect(data = subset(anno_df()$positions,
                                 transcript_ID == input$tab2.transcript_ID),
                   mapping = aes(xmin = as.integer(position)-0.5,
                                 xmax = as.integer(position)+0.5,
                                 fill = annotation),
                  ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha = 0.8)),
      position = "bottom")
    }
    return(p)
  })

#### tab2.hover ####

  output$tab2.hover <- renderUI({
    hover <- input$div_plot_hover
    point <- nearPoints(tab2.tableData(), hover, "Codon_Number", "Diversity", maxpoints=1)
    if (nrow(point) == 0) return(NULL)

    # calculate point position INSIDE the image as percent of total dimensions
    # from left (horizontal) and from top (vertical)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - log10(hover$y)) / (hover$domain$top - hover$domain$bottom)

    # calculate distance from left and bottom side of the picture in pixels
    left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
    top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)

    # create style property fot tooltip
    # background color is set so tooltip is a bit transparent
    # z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100;
                    background-color: rgba(245, 245, 245, 0.85); ",
                    "left:", left_px + 2, "px; top:", top_px + 2, "px;")

    wellPanel(style=style,
      p(HTML(paste0("<b>Codon: </b>", point$Codon_Number, "<br/>",
                    "<b>A_A_Change: </b>", point$Amino_Acid_Change, "<br/>",
                    "<b>Effect: </b>", point$Effect, "<br/>",
                    "<b>Diversity: </b>", point$Diversity)))

    )
  })



#### info ####
  output$info <- renderPrint({
    brushedPoints(tab2.tableData(), input$plot_brush, "Codon_Number", "Diversity")
  })

#### annotations

  ##                            _________
  ##                           /  tab3   \
  ## --------------------------           ----------------------------
  ## Tab 3 ##################################
#### tab3.selectGene ####
  output$tab3.selectGene <- renderUI({
    tagList(
      checkboxGroupInput("tab3.transcript_ID", label=NULL, choices=all.GeneChoices()),
      actionButton(inputId="tab3.Submit", label = "Submit")
    )
  })
#### tab3.Genes ####
  tab3.Genes <- eventReactive(input$tab3.Submit, {
    #gene Info for gene on tab 3, updates on 'submit' button press
    return(all.Genes()[ all.Genes()$transcript_ID %in% input$tab3.transcript_ID,])
  })

#### tab3.tidyData ####
  tab3.tidyData <- eventReactive(input$tab3.Submit, {
    data <- ldply(all.VCFList()[tab3.Genes()$transcript_ID])
    # remove 0|0 genotypes
    data <- data[data$gt_GT != "0|0",]
    return(data)
  })
#### tab3.EffectValues ####
  tab3.EffectValues <- reactive({
    # effects <- c("5_prime_UTR_variant",
    #              "intron_variant",
    #              "3_prime_UTR_variant",
    #              "synonymous_variant",
    #              "missense_variant",
    #              "upstream_gene_variant",
    #              "downstream_gene_variant")

    effects <- unique(tab3.tidyData()$Effect)

    return( switch(input$tab3.SNPtype,
            "All"=effects,
            "Missense"="missense_variant",
            "Coding"= c("missense_variant", "synonymous_variant"))
    )
  })
#### tab3.debug ####
  output$tab3.debug <- renderPrint({
    # temporary debug output
      print(input$tab3.filter_value)
  })
#### tab3.filteredByDiv ####
  tab3.filteredByDiv <- reactive({
    # filter by diversity slider and SNP type radio button then add SNPs column
    data <- tab3.tidyData()
    # filter by effect type (all, coding, or missense)
    data2 <- data[data$Effect %in% tab3.EffectValues(), ]
    # filter on positions with diversity greater than or equal to the 10^slider value
    keyPOS <- unique(data2[which(data2$Diversity >= 10^input$tab3.filter_value[1] &
                                data2$Diversity <= 10^input$tab3.filter_value[2]), "POS"])
    keydata <- data[data$POS %in% keyPOS, ]
    return(keydata)
  })
#### tab3.mutationList ####
  tab3.mutationList <- reactive({
    mutList <- labelBySNPs(tab3.filteredByDiv(), collapse=FALSE)$SNPs
    mutList <- unique(mutList[!is.na(mutList)])
    return(mutList)
  })
#### tab3.mutation_checkbox ####
  output$tab3.mutation_checkbox <- renderUI({
    tagList(
      tags$div(class="input-format",
          tags$h3("Allele selection"),
          tags$h5("Select the alleles you want to see on the map by clicking the checkboxes"),
          tags$div(class="checkbox-format",
                   checkboxGroupInput("tab3.allele_select", "select_alleles to display", choices=tab3.mutationList())
          ),
          actionButton(inputId="tab3.update_map", label = "Update Map")
      )
    )
  })
#### tab3.labeled ####
  tab3.labeled <- eventReactive(input$tab3.update_map, {
    # a dataframe with a single row per accession, containing accession info,
    # start with the data filtered by the diversity slider and type buttons
    data <- tab3.filteredByDiv()
    # label by SNPs creates column SNPs with text strings formatted [transcriptID|AA_Change]
    data <- labelBySNPs(data, collapse=FALSE)
    # filter on selected SNPs
    data <- data[data$SNPs %in% input$tab3.allele_select, ]
    # combine mutations to single row (this is slow)
    data <- ddply(data, "Indiv", summarise, SNPs=paste(SNPs, collapse=","))
    # add back ecotype details
    data <- addAccDetails(data, allAccs=TRUE)
    return(data)
  })
#### tab3.map ####
  output$tab3.map <- renderLeaflet({
    mapdata <- tab3.labeled()
    # Reorganize to plot NA's underneath non NA's
    mapdata <- rbind(mapdata[is.na(mapdata$SNPs), ], mapdata[!is.na(mapdata$SNPs), ])
    # make a field with text to be displayed when clicking on a marker
    mapdata$popup <- paste("EcoID:",  mapdata$Indiv,"Name:", mapdata$Name, " SNPs:", mapdata$SNPs)
    # create the color pallet for the map points
    pal <- brewer.pal(8, "Set1")
    pallet <- colorFactor(palette=pal, domain=mapdata$SNPs)
    # create a new leaflet map
    map <- leaflet()
    map <- addProviderTiles(map, providers$Stamen.TonerLite,
                     options = providerTileOptions(noWrap = TRUE))
    # groupnames to be used by draw groups of points as separate layers below
    groupnames <- unique(mapdata$SNPs)
    groupnames <- groupnames[!is.na(groupnames)]
    # add markers for NA points first so they are furthest back layer
    map <- addCircleMarkers(map, data=mapdata[is.na(mapdata$SNPs), ], color= "#9b9b9b", group="NA",
                            radius=6, popup= ~popup, stroke=FALSE, fillOpacity=0.6)
    # for each of the group names, add a set of markers
    for (SNP in groupnames){
          map <- addCircleMarkers(map, data=mapdata[mapdata$SNPs == SNP, ], color= ~pallet(SNPs), group= SNP,
                            radius=6, popup= ~popup, stroke=FALSE, fillOpacity=0.85)
    }
    # add the legend to the map
    map <- addLegend(map, position="bottomright", pal=pallet,
                     values=mapdata$SNPs, title="Marker Colors", opacity=1)
    # add layer control to map to turn on or off groups of points
    map <- addLayersControl(map, overlayGroups=c(groupnames, "NA"),
                            options = layersControlOptions(collapsed = TRUE),
                            position="bottomleft")
    return(map)
  })
#### tab3.dataTable ####
  output$tab3.dataTable <- DT::renderDataTable(tab3.labeled())
#### tab3.downloadMapData ####
  output$tab3.downloadMapData <- downloadHandler(
    filename=function(){
      paste("MapData-", Sys.time(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(tab3.labeled(), file, row.names=FALSE)
    }
  )



  ##                                        _________
  ##                                       /  tab4   \
  ## --------------------------------------           ----------------

  ## Tab 4 #####################

#### tab4.selectGene ####
  output$tab4.selectGene <- renderUI({
    tagList(
      checkboxGroupInput("tab4.transcript_ID", label=NULL, choices=all.GeneChoices()),
      actionButton(inputId="tab4.Submit", label = "Submit")
    )
  })
#### tab4.Genes ####
  tab4.Genes <- eventReactive(input$tab4.Submit, {
    #gene Info for gene on tab 3, updates on 'submit' button press
    return(all.Genes()[ all.Genes()$transcript_ID %in% input$tab4.transcript_ID,])
  })
#### tab4.tidyData ####
  tab4.tidyData <- eventReactive(input$tab4.Submit, {
    data <- ldply(all.VCFList()[tab4.Genes()$transcript_ID])
    data <- subset(data, select=-c(EFF, Transcript_ID, ID, FILTER ))
    data <- data[,filterTab.allCols]
    return(data)
  })
#### tab4.textFilters ####
  tab4.textFilters <- reactive({
    textFilters <- data.frame("filterID" = c("filter1", "filter2"),
                              "column" = c(input$tab4.filter1.column, input$tab4.filter2.column),
                              "values" = I(list(parseFilterText(input$tab4.filter1.textIn),
                                              parseFilterText(input$tab4.filter2.textIn))),
                              stringsAsFactors=FALSE)
  })
#### tab4.numFilters ####
  tab4.numFilters <- reactive({
    numFilters <- data.frame("filterID" = c("filter3", "filter4"),
                             "column" = c(input$tab4.filter3.column, input$tab4.filter4.column),
                             "max" = c(input$tab4.filter3.max, input$tab4.filter4.max),
                             "min" = c(input$tab4.filter3.min, input$tab4.filter4.min),
                             "missing" = c(input$tab4.filter3.missing, input$tab4.filter4.missing),
                             stringsAsFactors=FALSE)
  })
#### tab4.filteredVariants ####
  tab4.filteredVariants <- eventReactive(input$tab4.updateFilter,{
    # add all filtering here.
    data <- tab4.tidyData()
    if (input$tab4.filterRef) {
      # remove 0|0 genotypes
      data <- data[data$gt_GT != "0|0",]
    }
    for (i in 1:nrow(tab4.textFilters())){
      if (length(tab4.textFilters()[i,"values"][[1]]) > 0) {
        data <- data[as.character(data[, tab4.textFilters()[i, "column"]]) %in% tab4.textFilters()[i, "values"][[1]] , ]
      }
    }
    for (i in 1:nrow(tab4.numFilters())){
      naRows <- data[is.na(data[, tab4.numFilters()[i, "column"]]) , ]
      # remove NA rows to avoid issues with logical operators
      data <- data[!is.na(data[, tab4.numFilters()[i, "column"]]) , ]
      if (!is.na(tab4.numFilters()[i, "max"])){
        data <- data[  data[, tab4.numFilters()[i, "column"]] <=  tab4.numFilters()[i, "max"], ]
      }
      if (!is.na(tab4.numFilters()[i, "min"])){
        data <- data[  data[, tab4.numFilters()[i, "column"]] >=  tab4.numFilters()[i, "min"], ]
      }
      if (tab4.numFilters()[i,"missing"]){
        # add back NA rows if checkbox checked
        data <- rbind(data, naRows)
      }
    }
    return(data)
  })
#### tab4.debug ####
  output$tab4.debug <- renderPrint({
    print(tab4.numFilters())
  })
#### tab4.variantTable ####
  output$tab4.variantTable <- DT::renderDataTable(tab4.filteredVariants())
#### tab4.downloadVariantTable ####
  output$tab4.downloadVariantTable <- downloadHandler(
    filename=function(){
      paste("VariantTable-", Sys.time(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(tab4.filteredVariants(), file, row.names=FALSE)
    }
  )


  ##                                        _________
  ##                                      /   tab5   \
  ## --------------------------------------           ----------------
  ## Tab 5 #########################


#### tab5.selectGene ####
  output$tab5.selectGene <- renderUI({
    tagList(
      checkboxGroupInput("tab5.transcript_ID",
                         label=NULL, choices=all.GeneChoices()),
      actionButton(inputId="tab5.Submit", label = "Submit"),
      checkboxInput(inputId = "tab5.primary_transcript",
                         label = "Primary transcripts only?",
                         value = TRUE),
      radioButtons(inputId = "tab5.type",
                   label = "Alignment type:",
                   choices = c("DNA", "AA"),
                   selected = "AA", inline = TRUE)
    )
  })
#### tab5.Genes ####
  tab5.Genes <- eventReactive(input$tab5.Submit, {
    #gene Info for gene on tab 5, updates on 'submit' button press
    return(input$tab5.transcript_ID)
  })
#### debug ####
  output$tab5.debug <- renderPrint({
    aln_df()})
#### type ####
  type <- reactive({
    return(switch(input$tab5.type, "AA" = 2, "DNA" = 1))
  })
#### alignment ####
  alignment <- eventReactive(input$tab5.Submit, {
    alignment <- alignCDS(IDs = tab5.Genes(), primary_only = input$tab5.primary_transcript, all = {if(input$tab5.primary_transcript) FALSE else TRUE})
    return(alignment)
  })
#### tab5.alignment ####
  # output$tab5.alignment <- renderMsaR({
  #   msaR(alignment()[[type]], alignmentHeight = 100,
  #        colorscheme = {if(type) "taylor" else "nucleotide"})
  #   })
#### tab5.BrowseSeqs ####
  # output$tab5.BrowseSeqs <- reactive({
  #   file <- BrowseSeqs(alignment()[[type + 1]],
  #                      openURL = FALSE)
  #   html <- paste(readLines(file), collapse="\n")
  #   return(html)
  #   })
#### aln_df ####
  aln_df <- reactive({
    aln_df <- makeAlnDF(alignment()[[type()]])
    vcf <- ldply(.data = all.VCFList()[input$tab5.transcript_ID],
                 .fun = subset, !is.na(Transcript_ID) & gt_GT != "0|0")
    vcf <- getCodingDiv(vcf)
    aln_df <- addSNPsToAlnDF(aln_df, vcf)
    aln_df <- left_join(aln_df, dplyr::select(all.Genes(), "tair_locus",
                                       "tair_symbol", "transcript_ID"),
                        by = c("transcript_ID" = "transcript_ID"))
    ## chunk up aln_df
    aln_df <- chunkAlnDF(aln_df, chunk_width = 80)
    aln_df$seq_name <- as.character(aln_df$seq_name)
    aln_df$seq_name[!is.na(aln_df$tair_symbol)] <- aln_df$tair_symbol[!is.na(aln_df$tair_symbol)]
    aln_df$seq_name <- as.factor(aln_df$seq_name)
    print(aln_df)
    return(aln_df)
  })
#### tab5.aln_anno ####
  tab5.aln_anno <- reactive({
    ## read in annotation
    anno_df <- anno_df()

    anno_df <- addAlnPosToAnno(anno_df, aln_df())
    #print(anno_df)
    ## make chunks from aln_df
    chunks <- makeChunksDF(aln_df())
    ## chunk up annotations
    print(chunks)
    anno_df <- chunkAnnotation(anno_df, chunks)
    if(is.null(input$tab5.primary_transcript)) {
      anno_df$domains$seq_name <- factor(anno_df$domains$transcript_ID,
                       levels = levels(as.factor(aln_df()$transcript_ID)))
      anno_df$positions$seq_name <-
        factor(anno_df$positions$transcript_ID,
                        levels = levels(as.factor(aln_df()$transcript_ID)))}
    else {
      anno_df$domains$seq_name <- factor(anno_df$domains$tair_symbol,
                        levels = levels(as.factor(aln_df()$tair_symbol)))
      anno_df$positions$seq_name <- factor(anno_df$positions$tair_symbol,
                       levels = levels(as.factor(aln_df()$tair_symbol)))
    }
    print(anno_df)
    return(anno_df)
  })

#### aln_plot_height ####
  aln_plot_height <- reactive({
      N <- length(unique(aln_df()$seq_name))
      chunks <- length(unique(aln_df()$chunk))
      height <- 262 + 1.14*N + 19*chunks + 10*N*chunks
      return(ceiling(height))
    }
  )

#### tab5.aln_plot ####
  tab5.aln_plot <- reactive({
    p <-ggplot(aln_df(), aes(x = aln_pos, y = seq_name,
                             group = seq_pos, text = variants))
    if(!is.null(input$annoFile)) p <- p +
      geom_rect(data = tab5.aln_anno()$domains,
                mapping = aes(xmin = start_aln_pos - 0.5,
                              xmax = end_aln_pos + 0.5,
                              color = annotation,
                              ymin = as.numeric(seq_name)-0.5,
                              ymax = as.numeric(seq_name)+0.5),
                inherit.aes = FALSE, fill = NA, size = 1.2, alpha = 0.5) +
      geom_tile(data = tab5.aln_anno()$positions,
                mapping = aes(x = aln_pos, y = seq_name, color = annotation),
                width = 1, height = 1,
                fill = NA, size = 1.2, alpha = 0.5, inherit.aes = FALSE)
    p <- p +
      geom_tile(data = na.omit(aln_df()), mapping = aes(fill = effects),
                width = 1, height = 1, alpha = 0.8) +
      geom_text(aes(label=letter), alpha= 1, family = "Courier") +
      scale_fill_brewer(type = "qual", palette = 1, direction = -1) +
      scale_x_continuous(breaks=seq(1,max(aln_df()$aln_pos), by = 10)) +
      scale_y_discrete() +
      # expand increases distance from axis
      xlab("") +
      ylab("") +
      theme_logo(base_family = "Helvetica") +
      theme(panel.grid = element_blank(), panel.grid.minor = element_blank()) +
      facet_wrap(facets = ~chunk, ncol = 1, scales = "free") +
      theme(strip.background = element_blank(),
            strip.text.x = element_blank(),
            legend.box = "vertical")
    p})
  output$tab5.aln_plot <- renderPlot(expr = tab5.aln_plot() +
                                       theme(legend.position = "none"),
                                     res = 100)
  tab5.aln_plot_legend <- reactive({
    get_legend(tab5.aln_plot())
  })
  output$tab5.aln_plot_legend <- renderPlot(plot_grid(tab5.aln_plot_legend()),
                                            res = 100)

#### plot.ui ####
  output$plot.ui <- renderUI({
    plotOutput('tab5.aln_plot', height = aln_plot_height(),
               hover = hoverOpts("plot_hover", delay = 100,
                                 delayType = "debounce"))
  })

  #### aln_plot_hover ####
  output$aln_plot_hover <- renderUI({
    hover <- input$plot_hover
    point <- nearPoints(aln_df(), coordinfo = hover, xvar = "aln_pos",
                        yvar = "seq_name", panelvar1 = "chunk", threshold = 8,
                        maxpoints = 1, addDist = TRUE)
    if (nrow(point) == 0) return(NULL)

    # calculate point position INSIDE the image as percent of total dimensions
    # from left (horizontal) and from top (vertical)
    left_pct <- (hover$x - hover$domain$left) /
      (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) /
      (hover$domain$top - hover$domain$bottom)

    # calculate distance from left and bottom side of the picture in pixels
    left_px <- hover$range$left + left_pct *
      (hover$range$right - hover$range$left)
    right_px <- (1-left_pct) *
      (hover$range$right - hover$range$left)
    top_px <- hover$range$top + top_pct *
      (hover$range$bottom - hover$range$top)

    # create style property fot tooltip
    # background color is set so tooltip is a bit transparent
    # z-index is set so we are sure are tooltip will be on top
    if(left_pct < .70)
    style <- paste0("position:absolute; z-index:100;
                    background-color: rgba(245, 245, 245, 0.85); ",
                    "left:", left_px + 2, "px; top:", top_px + 2, "px;") else
          style <- paste0("position:absolute; z-index:100;
                    background-color: rgba(245, 245, 245, 0.85); ",
                    "right:", right_px + 10, "px; top:", top_px + 2, "px;")

    # actual tooltip created as wellPanel
    wellPanel(
      style = style,
      p(HTML(paste0("<b>symbol: </b>", point$seq_name, "<br/>",
                    "<b>transcript: </b>", point$transcript_ID, "<br/>",
                    "<b>seq_pos: </b>", point$seq_pos, "<br/>",
                    "<b>variants: </b>", point$variants)))
    )
  })
}


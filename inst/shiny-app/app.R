library(shiny)
library(biomaRt)
library(leaflet)
library(RColorBrewer)
library(shinythemes)
library(r1001genomes)
library(knitr)
library(stringr)

CSSCode <- tags$head(tags$style(
   HTML("
      .checkbox-format {
         -webkit-column-width: 350px;
         -moz-column-width: 350px;
         column-width: 350px;
      }

      .input-format {
         background-color: #dddddd;
         border: 1px solid #dddddd;
         border-radius: 12px;
         padding:1px 15px 10px 10px;
      }

      .output-format {
         border: 5px solid #dddddd;
         border-radius: 12px;
         padding:1px 15px 15px 20px;
      }

      .btn-default{
         color: #333;
         background-color: #eeeeee;
         border-color: #ccc;
      }

      .form-control{
         color: #333;
         background-color: #eeeeee;
         border-color: #ccc;
      }

      h1 {
         font-family: Helvetica;
         font-weight: 500;
         line-height: 1.1;
      }



   ")

))


filterTab.allCols <- c("Gene_Name", ".id", "Indiv", "POS", "Codon_Number", "gt_GT", "REF",
                       "gt_GT_alleles", "AC", "Effect", "Effect_Impact",
                       "Codon_Change", "Amino_Acid_Change", "Diversity")

filterTab.numericCols <- c("Indiv", "POS", "Codon_Number", "AC", "Diversity")


ui <- function(request){ fluidPage(

  tags$head(tags$style(
    HTML("
      .checkbox-format {
         -webkit-column-width: 350px;
         -moz-column-width: 350px;
         column-width: 350px;
      }
    ")
  )),


  #CSSCode,
  headerPanel("Arabidopsis Natural Variation Webtool"),
  "This app provides an interface to examine the natural variation of specified genes of interest in the 1001 Genomes project dataset. To save or share a state of this app, use the bookmark button.", HTML("</br>"),
  bookmarkButton(),
  tags$h5('style'="color:red", "This app is currently a work in progress."),
  themeSelector(),

  tags$br(),
  tags$div(class="input-format",
      fluidRow(
        column(6,
           tags$h3("Select Genes"),
           tags$h5("Type a list of gene loci in the box below, separated by commas. "),
           textAreaInput(inputId = "gene_ids", label = NULL,
                         width = "375px", height = 75, value = "AT3G62980, AT3G26810"),
           checkboxInput("STATS_quick_demo", label="Quick Demo"),
           actionButton(inputId="STATS_submit", label = "Submit")

        ),
        column(6,
           tags$h3("OR Upload File"),
           tags$h5("brows to a .csv file containing 'tair_locus' and 'name' fields.
                   The name field should be the TAIR symbol or moniker you would like to identify your genes by."),
           fileInput("genesFile", label=NULL),
           actionButton(inputId="file_submit", label = "Submit")

        )
      ),
      tags$br()
  ),


  tags$br(),



  tabsetPanel(
    tabPanel("SNP Stats",
        ## Tab 1 ###############################################################
      tags$hr(),
      tags$div(class="output-format",
               tags$h3("Gene Information"),
               tags$h5("This table provides details on the gene(s) input above, including transcript IDs and chromosomal locations."),
               downloadButton("tab1.downloadGeneInfo","Download Content of Table Below"),
               DT::dataTableOutput("tab1.genes_table")


      ),
      tags$br(),
      tags$div(class="output-format",
      tags$h3("Summary of Sequence Diversity"),
      downloadButton("tab1.downloadStats","Download Content of Tables Below"),
          tags$h4("Total Polymorphism Counts"),
          HTML("<h5>
                   This table provides counts of the total, non-unique polymorphisms present in the given genes by gene structure. These numbers can be quite high if the reference (Col-0) has a minor allele. \"coding_total\" is the sum of missense, nonsense and synonymous variants.
                </h5>"),
          DT::dataTableOutput("tab1.SNPcounts"),
          tags$h4("Unique Allele Counts"),
          HTML("<h5>
              This table provides counts of unique alleles by gene structure.
               </h5>"),
          DT::dataTableOutput("tab1.SNPcountsUnique"),
          tags$h4("Nucleotide Diversity Statistics"),
          HTML("<h5>
               This table provides for each given gene the nucleotide diversity as Nei and Li's <i>&pi;</i> (the average number of nucleotide differences per site between all possible pairs of sequence) at synonymous (pi_s) and missense (pi_n) sites.
                </br>
               In the future we plan to add
               Fixation index (<i>F<sub>ST</sub></i>),
               Tajima's <i>D</i> and
               Watterson's <i>&theta;</i> to this table."),
          DT::dataTableOutput("tab1.Diversity_table")
      )
    ),

    tabPanel("Diversity Plot",
        ## Tab 2 ###############################################################
      tags$br(),
      tags$div(class="input-format",
               tags$h3("Select a Gene"),
               tags$h5("Select a transcript ID in the box below"),
               uiOutput("tab2.selectGene")
               # textInput(inputId = "plotGene", label =NULL,
               #           value = "AT1G80490"),
               #actionButton(inputId="tab2.Submit", label = "Submit")
      ),

      tags$hr(),
      tags$div(class="output-format",
          tags$h3("Selected Gene Information"),
          tableOutput("tab2.GeneInfo")
      ),
      tags$br(),

      tags$div(class="output-format",
          tags$h3("Plot of Nucleotide Diversity Statistic by Codon"),
          tags$h5("To see details on specific points, click and drag to create a box selecting points."),
          plotOutput("diversityPlot", brush="plot_brush", click="plot_click", height = 400),
          verbatimTextOutput("info")
      ),
      tags$br(),
      tags$div(class="output-format",

          tags$h3("Diversity Plot Data"),
          tags$h5("This table provides the raw data from the plot. \"POS\" is the chromosomal position of the SNP,
                  the Amino_Acid_Change field provides both the amino acid as well as the base change"),
          downloadButton("tab2.downloadSNPData","Download Content of Table Below"),
          DT::dataTableOutput("Diversity_Table")
      )
    ),

    tabPanel("SNP Mapping",
             ## Tab 3 ##########################################################
             tags$br(),
             tags$div(class="input-format",
                 tags$h3("Select Genes and Filter Diversity Parameter"),
                 tags$h5("Select one or more transcript IDs below and use the slider to select a minimum sitewise nucleotide diversity"),
                 # textInput(inputId="tab3.Gene", label=NULL,
                 #           value="AT1G80490"),
                 uiOutput("tab3.selectGene"),
                 # actionButton(inputId="tab3.Submit", label="Submit"),
                 sliderInput(inputId="tab3.filter_value", label="Log Nucleotide diversity filter limit",
                             min=-4, max=0, value=-2, step=0.05),
                 radioButtons("tab3.SNPtype", "Type of SNP to mark",
                              choices=c("All", "Coding", "Missense"))
                 #verbatimTextOutput("tab3.debug")
             ),

             tags$br(),

             uiOutput("tab3.mutation_checkbox"),

             tags$hr(),
             tags$div(class="output-format",
             tags$h3("Accession Map"),
             tags$h5("Zoom with scroll wheel, click and drag to pan, click on individual point to see details.
                     Use the layers pop-out to the lower left of the map to hide or show sets of accessions with the same variant."),
             leafletOutput("tab3.map", width = "95%")
             ),
             tags$br(),
             tags$div(class="output-format",
                 tags$h3("Map Data"),
                 downloadButton("tab3.downloadMapData","Download Content of Table Below"),
                 DT::dataTableOutput("tab3.dataTable")
             )
    ),


    tabPanel("Accessions and Mutations",
              ## Tab 4 #########################################################
             tags$br(),
             tags$div(class="input-format",
                      tags$h3("Gene Select"),
                      tags$h5("select one or more transcipt IDs below"),
                      uiOutput("tab4.selectGene")
             ),
             tags$br(),
             tags$div(class="input-format",
                        tags$h3("Filters"),
                        checkboxInput("tab4.filterRef", "hide 0|0 genotypes?", FALSE),

                      fluidRow(

                        column(3, wellPanel(
                          tags$h4("Filter 1"),
                          tags$br(),
                          tags$h5("select a column to filter on"),
                          selectInput("tab4.filter1.column", label="column select",
                                      choices=filterTab.allCols),
                          tags$br(),
                          tags$h5("values to match. separate values with a comma followed by a space \n(ie. \"a, b\") "),
                          textAreaInput("tab4.filter1.textIn", NULL)
                        )),

                        column(3, wellPanel(
                          tags$h4("Filter 2"),
                          tags$br(),
                          tags$h5("select a column to filter on"),
                          selectInput("tab4.filter2.column", label="column select",
                                      choices=filterTab.allCols),
                          tags$br(),
                          tags$h5("values to match. separate values with a comma followed by a space \n(ie. \"a, b\") "),
                          textAreaInput("tab4.filter2.textIn", NULL)
                        )),

                        column(3, wellPanel(
                          tags$h4("Filter 3 (Numeric)"),
                          tags$br(),
                          tags$h5("select a column to filter on"),
                          selectInput("tab4.filter3.column", label="column select",
                                      choices=filterTab.numericCols),
                          tags$br(),
                          tags$h5("Max value"),
                          numericInput("tab4.filter3.max", NULL, NA),
                          tags$h5("Min Value"),
                          numericInput("tab4.filter3.min", NULL, NA),
                          checkboxInput("tab4.filter3.missing", "keep rows with missing values?")
                        )),

                        column(3, wellPanel(
                          tags$h4("Filter 4 (Numeric)"),
                          tags$br(),
                          tags$h5("select a column to filter on"),
                          selectInput("tab4.filter4.column", label="column select",
                                      choices=filterTab.numericCols),
                          tags$br(),
                          tags$h5("Max value"),
                          numericInput("tab4.filter4.max", NULL, NA),
                          tags$h5("Min Value"),
                          numericInput("tab4.filter4.min", NULL, NA),
                          checkboxInput("tab4.filter4.missing", "keep rows with missing values?")
                        ))
                      ),
                      actionButton(inputId="tab4.updateFilter", label = "Apply Filters")


             ),
             tags$hr(),
             verbatimTextOutput("tab4.debug"),
             tags$div(class="output-format",
                      tags$h3("Filtered Variants"),
                      tags$h5("This table provides ..."),
                      downloadButton("tab4.downloadVariantTable","Download Content of Table Below"),
                      DT::dataTableOutput("tab4.variantTable")


             )

    ),

    tabPanel("About",
             ## About Tab ######################################################
             tags$br(),
             column(6,
                    tags$div(class="output-format",
                             includeHTML("Glossary.html")
                    )
             ),
             column(6,
                    tags$div(class="output-format",
                             includeHTML("Bibliography.html")
                    )
             )


    )

  ) #end of tabset panel

  # "THIS IS THE FOOTER"
)}


# Server =================================================================

#source("VCF_Utils.R")
#source("Strains_and_Gene_Families.R")

parseInput <- function (textIn) {
  names <- str_extract_all(textIn, "AT[1-5]G[0-9]{5}")
  return (names[[1]])
}

parseFilterText <- function (textIn) {
  inputList <- strsplit(textIn, ", ")
  inputList <- gsub(" ", "", inputList[[1]]) # remove extra spaces
  return(inputList)
}




#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



server <- function(input, output){

  ##   _________
  ##  /  tab1   \
  ##             --------------------------------------------------
  ## Tab 1 stuff:

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

  all.GeneChoices <- reactive({
    displayNames <- paste(all.Genes()$transcript_ID, " (", all.Genes()$tair_symbol, ")", sep="" )
    displayNames <- gsub(" \\(\\)", displayNames, replacement="")  # if no tair symbol, remove empty parens.
    output <- all.Genes()$transcript_ID
    names(output) <- displayNames
    return(output)
  })

  output$tab1.genes_table <- DT::renderDataTable(DT::datatable(all.Genes()[, -c(5,6,7,10)], colnames = c("tair locus", "symbol", "transcript", "Chr", "transcript \nstart", "transcript \nend", "transcript \nlength"), rownames = FALSE, options=list(paging=FALSE, searching=FALSE)))

  all.VCFList <- reactive({

    if(isolate(input$STATS_quick_demo)) {
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
                     setProgress(value=1)
    })

    return(output)
  })

  tab1.nonUniqueVariants <- eventReactive({all.VCFList()},{
    req(isolate(tab1.buttons$last_button)!="none pressed")
    ldply(all.VCFList(), variantCounts, unique=FALSE, .id="transcript_ID")
  })

  tab1.uniqueVariants <- eventReactive({all.VCFList()},{
    req(isolate(tab1.buttons$last_button)!="none pressed")
    ldply(all.VCFList(), variantCounts, unique=TRUE, .id="transcript_ID")
  })

  tab1.divStats <- eventReactive({all.VCFList()},{
    req(isolate(tab1.buttons$last_button)!="none pressed")
    ldply(all.VCFList(), diversityStats, geneInfo=isolate(all.Genes()), .id="transcript_ID")
  })

  SNPStats <- reactive({
    req(isolate(tab1.buttons$last_button)!="none pressed")
    # rename column names on unique variant counts.
    uniqueVariantsRenamed <- tab1.uniqueVariants()
    colnames(uniqueVariantsRenamed) <- paste(colnames(uniqueVariantsRenamed),
                                             "unique", sep="_")
    cbind(tab1.nonUniqueVariants(), uniqueVariantsRenamed[, -1], tab1.divStats()[, -1])
  })

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
                columns = 2:6, digits = 6))

  output$tab1.downloadStats <- downloadHandler(
    filename=function(){
      paste("SNPStats-", Sys.time(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(SNPStats(), file, row.names=FALSE)
    }
  )

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
  ## Tab 2 stuff:

  output$tab2.selectGene <- renderUI({
    tagList(
      selectInput("tab2.transcript_ID", label=NULL, choices=all.GeneChoices()),
      actionButton(inputId="tab2.Submit", label = "Submit")
    )
  })

  tab2.Genes <- eventReactive(input$tab2.Submit, {
      #gene Info for gene on tab 2, updates on 'submit' button press
    # names <- parseInput(input$plotGene)
    # genes <- getGeneInfo(names[1])
    # return(genes)
    return(all.Genes()[ all.Genes()$transcript_ID == input$tab2.transcript_ID,])
  })

  output$tab2.GeneInfo <- renderTable(tab2.Genes())
    #rendered table of Gene info

  #tab2.tableData <- reactive({load_tab_2_Data(tab2.Genes())})
    #SNP reactive data
  tab2.tableData <- eventReactive(input$tab2.Submit, {
    tab2data <- all.VCFList()[[input$tab2.transcript_ID]]
    coding_variants <- getCodingDiv(tab2data)
    return(coding_variants)
  })


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

  output$diversityPlot <- renderPlot(plotCodingDiv(tab2.tableData()))
    #plot output

  output$info <- renderPrint({
    brushedPoints(tab2.tableData(), input$plot_brush, "Codon_Number", "Diversity")
  })


  ##                            _________
  ##                           /  tab3   \
  ## --------------------------           ----------------------------
  ## Tab 3 stuff:

  output$tab3.selectGene <- renderUI({
    tagList(
      checkboxGroupInput("tab3.transcript_ID", label=NULL, choices=all.GeneChoices()),
      actionButton(inputId="tab3.Submit", label = "Submit")
    )
  })





  tab3.Genes <- eventReactive(input$tab3.Submit, {
    #gene Info for gene on tab 3, updates on 'submit' button press
    return(all.Genes()[ all.Genes()$transcript_ID %in% input$tab3.transcript_ID,])
  })


  tab3.tidyData <- eventReactive(input$tab3.Submit, {
    data <- ldply(all.VCFList()[tab3.Genes()$transcript_ID])

    # remove 0|0 genotypes
    data <- data[data$gt_GT != "0|0",]

    return(data)
  })

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

  output$tab3.debug <- renderPrint({
    # temporary debug output
      print(paste("last button =", tab1.buttons$last_button))
      print(paste("total presses =", tab1.buttons$total_presses))
  })

  tab3.filteredByDiv <- reactive({
    # filter by diversity slider and SNP type radio button then add SNPs column

    data <- tab3.tidyData()

    # filter by effect type (all, coding, or missense)
    data2 <- data[data$Effect %in% tab3.EffectValues(), ]

    # filter on positions with diversity greater than or equal to the 10^slider value
    keyPOS <- unique(data2[which(data2$Diversity >= 10^input$tab3.filter_value), "POS"])
    keydata <- data[data$POS %in% keyPOS, ]

    return(keydata)
  })

  tab3.mutationList <- reactive({
    mutList <- labelBySNPs(tab3.filteredByDiv(), collapse=FALSE)$SNPs
    mutList <- unique(mutList[!is.na(mutList)])
    return(mutList)
  })

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
    data <- addAccDetails(data)

    return(data)
  })

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

  output$tab3.dataTable <- DT::renderDataTable(tab3.labeled())

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
  ## Tab 4 stuff:

  output$tab4.selectGene <- renderUI({
    tagList(
      checkboxGroupInput("tab4.transcript_ID", label=NULL, choices=all.GeneChoices()),
      actionButton(inputId="tab4.Submit", label = "Submit")
    )
  })

  tab4.Genes <- eventReactive(input$tab4.Submit, {
    #gene Info for gene on tab 3, updates on 'submit' button press
    return(all.Genes()[ all.Genes()$transcript_ID %in% input$tab4.transcript_ID,])
  })

  tab4.tidyData <- eventReactive(input$tab4.Submit, {
    data <- ldply(all.VCFList()[tab4.Genes()$transcript_ID])
    data <- subset(data, select=-c(EFF, Transcript_ID, ID, FILTER ))
    data <- data[,filterTab.allCols]
    return(data)
  })

  tab4.textFilters <- reactive({
    textFilters <- data.frame("filterID" = c("filter1", "filter2"),
                              "column" = c(input$tab4.filter1.column, input$tab4.filter2.column),
                              "values" = I(list(parseFilterText(input$tab4.filter1.textIn),
                                              parseFilterText(input$tab4.filter2.textIn))),
                              stringsAsFactors=FALSE)
  })

  tab4.numFilters <- reactive({
    numFilters <- data.frame("filterID" = c("filter3", "filter4"),
                             "column" = c(input$tab4.filter3.column, input$tab4.filter4.column),
                             "max" = c(input$tab4.filter3.max, input$tab4.filter4.max),
                             "min" = c(input$tab4.filter3.min, input$tab4.filter4.min),
                             "missing" = c(input$tab4.filter3.missing, input$tab4.filter4.missing),
                             stringsAsFactors=FALSE)


  })


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

  output$tab4.debug <- renderPrint({
    print(tab4.numFilters())
  })



  output$tab4.variantTable <- DT::renderDataTable(tab4.filteredVariants())

  output$tab4.downloadVariantTable <- downloadHandler(
    filename=function(){
      paste("VariantTable-", Sys.time(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(tab4.filteredVariants(), file, row.names=FALSE)
    }
  )


}

enableBookmarking(store = "url")
shinyApp(ui = ui, server = server)


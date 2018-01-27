library(shiny)
library(biomaRt)
library(leaflet)
library(RColorBrewer)
library(shinythemes)
library(r1001genomes)
library(knitr)
library(stringr)
library(msaR)
library(DECIPHER)
library(plotly)
library(ggseqlogo)
library(shinyBS)

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

      .wrapper {
        background:#EFEFEF;
        box-shadow: 1px 1px 10px #999;
        margin: auto;
        text-align: center;
        position: relative;
        -webkit-border-radius: 5px;
        -moz-border-radius: 5px;
        border-radius: 5px;
        margin-bottom: 20px !important;
        width: 800px;
        padding-top: 5px;
        }

      .scrolls {
        overflow-x: scroll;
        overflow-y: hidden;
        height: 80px;
        white-space:nowrap
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
  theme = shinytheme(theme = "flatly"),
  tags$head(tags$style(
    HTML("
      .checkbox-format {
         -webkit-column-width: 350px;
         -moz-column-width: 350px;
         column-width: 350px;
      }
    ")
  )),
  #### Header ####
  #CSSCode,
  headerPanel("Arabidopsis Natural Variation Webtool"),
  "This app provides an interface to examine the natural variation of specified genes of interest in the 1001 Genomes project dataset. To save or share a state of this app, use the bookmark button.", HTML("</br>"),
  bookmarkButton(),
  tags$h5('style'="color:red", "This app is currently a work in progress."),
  # themeSelector(),
  tags$br(),
  bsCollapse(id = "collapse 1", multiple=TRUE, open=c("Gene Select", "Annotation Files"),
    bsCollapsePanel("Gene Select",
                    fluidRow(
                      column(5,
                             tags$h3("Select Genes"),
                             tags$h5("Type a list of gene loci in the box below, separated by commas. "),
                             textAreaInput(inputId = "gene_ids", label = NULL,
                                           width = "375px", height = 75, value = "AT3G62980, AT3G26810"),
                             checkboxInput("STATS_quick_demo", label="Quick Demo"),
                             actionButton(inputId="STATS_submit", label = "Submit")

                      ),
                      column(2, align="center", tags$h3("OR")),
                      column(5,
                             tags$h3("Upload File"),
                             tags$h5("brows to a .csv file containing 'tair_locus' and 'name' fields.
                                     The name field should be the TAIR symbol or moniker you would like to identify your genes by."),
                             fileInput("genesFile", label=NULL),
                             actionButton(inputId="file_submit", label = "Submit")

                             )
                    )
    ),
    bsCollapsePanel("Annotation Files",
      fluidRow(
        column(5,
            tags$h3("Upload an annotation file"),
            tags$h5("Browse to a '.csv' file containing a 'gene' field matching the tair loci or gene symbols of your genes of interest.
                    Or create a new annotation file by downloading the empty template file and adding annotations to it.")
        ),
        column(4,
               fileInput("annoFile", label="Upload Annotation File:"),
               tags$h5(tags$strong("Download Empty Annotation File Template:")),
               downloadButton("annoTemplateDownload",label="Download Template")
        ),
        column(3,
               tags$h5(tags$strong("Submit Uploaded Annotation File:")),
               actionButton(inputId="annoSubmit", label = "Submit")
        )
      )
    )

  ),


  tags$br(),
  tabsetPanel(
    tabPanel("SNP Stats",
             ## Tab 1 ###############################################################
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
                      tags$hr(),
                      tags$h4("Unique Allele Counts"),
                      HTML("<h5>
              This table provides counts of unique alleles by gene structure.
               </h5>"),
                      DT::dataTableOutput("tab1.SNPcountsUnique"),
                      tags$hr(),
                      tags$h4("Nucleotide Diversity Statistics"),
                      HTML("<h5>
               This table provides for each given gene the nucleotide diversity as Nei and Li's <i>&pi;</i> (the average number of nucleotide differences per site between all possible pairs of sequence) at synonymous (<i>&pi;<sub>S</sub></i>) and nonsynonymous (<i>&pi;<sub>N</sub></i>) sites.
                </br>
               In the future we plan to add
               Fixation index (<i>F<sub>ST</sub></i>),
               Tajima's <i>D</i> and
               Watterson's <i>&theta;</i> to this table."),
                      DT::dataTableOutput("tab1.Diversity_table")
             )
    ),
    ## Tab 2 - Diversity Plot #####################################################
    tabPanel("Diversity Plot",
             tags$br(),
             tags$div(class="input-format",
                      tags$h3("Select a Gene"),
                      tags$h5("Select a transcript ID in the box below"),
                      uiOutput("tab2.selectGene")
             ),
             tags$hr(),
             tags$div(class="output-format",
                      tags$h3("Selected Gene Information"),
                      DT::dataTableOutput("tab2.gene_table")
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
    ## Tab 3  - SNP Mapping #######################################################
    tabPanel("SNP Mapping",
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
    ## Tab 4 #########################################################
    tabPanel("SNP Browser",
             tags$br(),
             tags$div(class="input-format",
                      tags$h3("Gene Select"),
                      tags$h5("select one or more transcipt IDs below"),
                      uiOutput("tab4.selectGene")
             ),
             tags$br(),
             tags$div(class="input-format",
                      tags$h3("Filters"),
                      tags$h5("NOTE: all filters are combined by a logical AND.
                                So for a row to be displayed, it must satisfy the requirements of ALL the filters."),
                      checkboxInput("tab4.filterRef", "hide 0|0 genotypes?", FALSE),
                      #### Filter 1 ####
                      wellPanel(fluidRow(
                        column(2,tags$h4("Filter 1")),
                        column(3,
                               selectInput("tab4.filter1.column", label="column select",
                                           choices=filterTab.allCols)
                        ),
                        column(5,
                               textInput("tab4.filter1.textIn", "values to match")
                        ),
                        column(2,
                               tags$h5("Separate values with a comma followed by a space \n(ie. \"a, b\"). ")
                        )
                      )),
                      #### Filter 2 ####
                      wellPanel(fluidRow(
                        column(2, tags$h4("Filter 2")),
                        column(3,
                               selectInput("tab4.filter2.column", label="column select",
                                           choices=filterTab.allCols)
                        ),
                        column(5,
                               textInput("tab4.filter2.textIn", "values to match")
                        ),
                        column(2,
                               tags$h5("Separate values with a comma followed by a space \n(ie. \"a, b\"). ")
                        )
                      )),
                      #### filter 3 ####
                      wellPanel(fluidRow(
                        column(2, tags$h4("Filter 3"),tags$h4("(Numeric)")),
                        column(3,
                               selectInput("tab4.filter3.column", label="column select",
                                           choices=filterTab.numericCols)
                        ),
                        column(2, numericInput("tab4.filter3.min", "MIN", NA)),
                        column(2, numericInput("tab4.filter3.max", "MAX", NA)),
                        column(3,
                               checkboxInput("tab4.filter3.missing", "keep rows with missing values?")
                        )
                      )),
                      #### Filter 4 ####
                      wellPanel(fluidRow(
                        column(2, tags$h4("Filter 4"),tags$h4("(Numeric)")),
                        column(3,
                               selectInput("tab4.filter4.column", label="column select",
                                           choices=filterTab.numericCols)
                        ),
                        column(2, numericInput("tab4.filter4.min", "MIN", NA)),
                        column(2, numericInput("tab4.filter4.max", "MAX", NA)),
                        column(3,
                               checkboxInput("tab4.filter4.missing", "keep rows with missing values?")
                        )
                      )),

                      actionButton(inputId="tab4.updateFilter", label = "Apply Filters")

             ),
             tags$hr(),
             # verbatimTextOutput("tab4.debug"),  un-comment to debug
             tags$div(class="output-format",
                      tags$h3("Filtered Variants"),
                      tags$h5("This table provides ..."),
                      downloadButton("tab4.downloadVariantTable","Download Content of Table Below"),
                      DT::dataTableOutput("tab4.variantTable")

             )

    ),

    ## Tab 5 - Alignments #########################################################
    tabPanel("Alignments",
             tags$br(),
             tags$div(class="input-format",
                      tags$h3("Select Genes and Type"),
                      tags$h5("Select one or more transcript IDs below and the type of alignment to show"),
                      uiOutput("tab5.selectGene")),
             tags$br(),
             tags$hr(),
             tags$div(class = "output-format",
                      tags$h3("Sequence Alignment"),
                      tags$h5("Click and drag to pan. The x-axis is the position within the alignment. Hover over the alignment to see details. 'seq_pos' is the position in the sequence with name 'seq_name' of the type chosen above. Use the pop-up menu in the upper right for zoom and other plotly functionalities. Made with",
                              tags$a(href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0749-z", target = "_blank", "DECIPHER")),
                      plotlyOutput('tab5.aln_plot', height = "auto"),
                      # verbatimTextOutput("event")
                      tags$br())#,
             #tags$h5("Click and drag to pan. Made with", tags$a(href="https://zachcp.github.io/msaR/", "msaR")),
             #msaROutput(outputId = "tab5.alignment"), height = "auto")
             # Remove DECIPHER BrowseSeqs
             # but perhaps the side scrolling div will come in handy again
             # tags$div(class = "wrapper",
             #          tags$div(class = "scrolls",
             #            htmlOutput("tab5.BrowseSeqs", inline = TRUE)
             #   )
             # ),
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
#### Annotation Template Download ####
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
    if(!is.null(input$anno_df)){
    p <- append_layers(p,
      geom_rect(data = subset(input$anno_df, gene == tab2.selectGene),
                mapping = aes(xmin = as.integer(start),
                              xmax = as.integer(end),
                              fill = annotation),
                ymin = -Inf, ymax = Inf, inherit.aes = FALSE),
      position = "bottom")
    }
    return(p)
  })
    #plot output
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
      print(paste("last button =", tab1.buttons$last_button))
      print(paste("total presses =", tab1.buttons$total_presses))
  })
#### tab3.filteredByDiv ####
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
    data <- addAccDetails(data)
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
      checkboxGroupInput(inputId = "tab5.transcript_ID",
                         label=NULL, choices=all.GeneChoices()),
      radioButtons(inputId = "tab5.type",
                   label = "Alignment type:",
                   choices = c("DNA", "AA"),
                   selected = "AA", inline = TRUE),
      actionButton(inputId="tab5.Submit", label = "Submit")
    )
  })
#### tab5.Genes ####
  tab5.Genes <- eventReactive(input$tab5.Submit, {
    #gene Info for gene on tab 5, updates on 'submit' button press
    return(input$tab5.transcript_ID)
  })
#### debug ####
  output$debug <- renderPrint({tab5.Genes()})
#### type ####
  type <- reactive({
    return(switch(input$tab5.type, "AA" = 2, "DNA" = 1))
  })
#### alignment ####
  alignment <- eventReactive(input$tab5.Submit, {
    alignment <- alignCDS(IDs = tab5.Genes())
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
    ## chunk up aln_df
    chunk_width <- 80
    chunk_num <- round(max(aln_df$aln_pos)/chunk_width, 1)
    aln_df$chunk <- cut_number(aln_df$aln_pos, n = chunk_num)
    return(aln_df)
  })
#### tab5.aln_anno ####
  # tab5.aln_anno <- reactive({
  #   ## read in annotation
  #   anno
  #   ## subset to domains
  #   anno.domain
  #   ## consolidate to aln_pos
  #   anno.domain
  #   ## convert annotations to type
  #   anno.domain
  #   ## chunk up annotations
  #   ### make chunks df
  #   chunks <- data.frame("chunk" = levels(aln_df()$chunk)) %>%
  #     separate(col = "chunk", into = c("start", "end"), sep = ",", remove=FALSE)
  #   chunks$start <- if_else(condition = str_detect(chunks$start, "\\("), true = as.numeric(str_extract(chunks$start, "\\d+")) + 1, false = as.numeric(str_extract(chunks$start, "\\d+")))
  #   chunks$end <- if_else(condition = str_detect(chunks$end, "\\)"), true = as.numeric(str_extract(chunks$end, "\\d+")) - 1, false = as.numeric(str_extract(chunks$end, "\\d+")))
  #   ### chunk annotations
  #   chunks.anno.domain <- adply(anno.domain, 1, function(domain) {
  #     #check which chunks it spans
  #     rows <- (domain$start < chunks$start & domain$end > chunks$end) |
  #       (domain$start < chunks$end & domain$start > chunks$start) |
  #       (domain$end > chunks$start & domain$end < chunks$end)
  #     print(rows)
  #     if(sum(rows) != 0){
  #       #create copies of the row for each chunk
  #       new_rows <- domain[rep(1, sum(rows)),]
  #       new_rows[,c("chunk", "start", "end")] <-
  #         chunks[rows, c("chunk", "start","end")]
  #       print(new_rows)
  #       new_rows[(domain$start < chunks$end & domain$start > chunks$start),
  #                "start"] <- domain[, "start"]
  #       new_rows[(domain$end > chunks$start & domain$end < chunks$end),
  #                "end"] <- domain[, "end"]
  #       new_rows
  #     }
  #     # or if mono-chunk-ular add chunk info
  #     else{
  #       domain[,"chunk"] <-
  #         chunks[which(chunks$start <= domain$start & chunks$end >= domain$end),
  #                c("chunk")]
  #       domain
  #     }
  #   })
  #   return(chunks.anno.domain)
  # })
#### tab5.aln_plot ####
  output$tab5.aln_plot <- renderPlot({
    p <-
      ggplot(aln_df(), aes(x = aln_pos, y = seq_name,
                           group = seq_pos, text = variants)) +
      # geom_rect(data = tab5.aln_anno(),
      #           mapping = aes(xmin = start - 0.5,
      #                         xmax = end + 0.5,
      #                         fill = annotation),
      #           ymin = -Inf, ymax = Inf, inherit.aes = FALSE) +
      geom_tile(data = na.omit(aln_df()), mapping = aes(fill = effects),
                width = 1, height = 1, alpha = 0.8) +
      geom_text(aes(label=letter), alpha= 1) +
      scale_fill_brewer(type = "qual", palette = 1, direction = -1) +
      scale_x_continuous(breaks=seq(1,max(aln_df()$aln_pos), by = 10)) +
      # expand increases distance from axis
      xlab("") +
      ylab("") +
      #scale_size_manual(values=c(5, 6)) + # does nothing unless 'size' is mapped
      theme_logo(base_family = "Courier") +
      theme(panel.grid = element_blank(), panel.grid.minor = element_blank()) +
      facet_wrap(facets = ~chunk, ncol = 1, scales = "free") +
      theme(strip.background = element_blank(),
            strip.text.x = element_blank())
    return(p)
    # ggplotly(p, tooltip = c("seq_name", "seq_pos", "variants"),
    #          height = length(unique(aln_df()$seq_name)) * 20 + 110) %>%
    #   # Total height = 125 = N*50 + 25
    #   # 125 = N*30 + 65 better still a bit long with 6 sequences
    #   # 125 = N*20 + 85 looks really good but maybe a bit
    #   config(collaborate = FALSE) %>%
    #   layout(margin = list(l = 100, r = 0, t = 20, b = 0),
    #          legend = list(yanchor = "bottom", y = -1, orientation = "h"),
    #          dragmode = "pan", yaxis = list(fixedrange = TRUE),
    #          xaxis = list(range = c(0,70)))
  })
}

enableBookmarking(store = "url")
shinyApp(ui = ui, server = server)


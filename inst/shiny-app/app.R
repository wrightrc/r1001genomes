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
         color: #48ca3b;
         background-color: #dce4f2;
         border: 10px solid #dce4f2;
         border-radius: 12px;
      }



   ")

))


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


  #CSSCode,
  headerPanel("Arabidopsis Natural Variation Webtool"),
  "This app provides an interface to examine the natural variation of specified genes of interest in the 1001 Genomes project dataset. To save or share a state of this app, use the bookmark button.", HTML("</br>"),
  bookmarkButton(),
  tags$h5('style'="color:red", "This app is currently a work in progress."),
  #themeSelector(),
  tabsetPanel(
    tabPanel("SNP Stats",
## Tab 1 - SNP Stats ##########################################################
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
               tags$h5("Browse to a '.csv' file containing 'tair_locus' and 'name' fields.
                       The 'name' field should be the TAIR symbol or moniker you would like to identify your genes by."),
               fileInput("genesFile", label=NULL),
               actionButton(inputId="file_submit", label = "Submit")

            )
          ),
          tags$br()
      ),

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
               This table provides for each given gene the nucleotide diversity as Nei and Li's <i>&pi;</i> (the average number of nucleotide differences per site between all possible pairs of sequence) at synonymous (<i>&pi;<sub>S</sub></i>) and nonsynonymous (<i>&pi;<sub>N</sub></i>) sites.
                </br>
               In the future we plan to add
               Fixation index (<i>F<sub>ST</sub></i>),
               Tajima's <i>D</i> and
               Watterson's <i>&theta;</i> to this table."),
          DT::dataTableOutput("tab1.Diversity_table")
      )
    ),

    tabPanel("Diversity Plot",
## Tab 2 - Diversity Plot #####################################################
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
## Tab 3  - SNP Mapping #######################################################
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


    # tabPanel("Accessions and Mutations",
    ### Tab 4  - Accessions and Mutations #####################################
    #          tags$br(),
    #          tags$div(class="input-format",
    #                   tags$h3("Gene Select"),
    #                   tags$h5("select one or more transcipt IDs below"),
    #                   uiOutput("tab4.selectGene"),
    #                   tags$hr(),
    #                   tags$br()
    #
    #          ),
    #          tags$br(),
    #          tags$div(class="input-format",
    #                     tags$h3("Filters"),
    #                     checkboxInput("tab4.filterRef", "hide 0|0 genotypes?", FALSE),
    #
    #                   fluidRow(
    #
    #
    #                     column(3, wellPanel(
    #                            tags$h4("Filter 1"),
    #                            tags$br(),
    #                            tags$h5("select a column to filter on"),
    #                            selectInput("tab4.filter1.column", label="column select", choices=c("POS", "gt_GT", "...")),
    #                            tags$br(),
    #                            tags$h5("values to match, separated by commas"),
    #                            textAreaInput("tab4.filter1.textIn", NULL)
    #                     )),
    #
    #                     column(3, wellPanel(
    #                       tags$h4("Filter 2"),
    #                       tags$br(),
    #                       tags$h5("select a column to filter on"),
    #                       selectInput("tab4.filter2.column", label="column select", choices=c("POS", "gt_GT", "...")),
    #                       tags$br(),
    #                       tags$h5("values to match, separated by commas"),
    #                       textAreaInput("tab4.filter2.textIn", NULL)
    #                     )),
    #
    #                     column(3, wellPanel(
    #                       tags$h4("Filter 3 (Numeric)"),
    #                       tags$br(),
    #                       tags$h5("select a column to filter on"),
    #                       selectInput("tab4.filter3.column", label="column select", choices=c("POS", "gt_GT", "...")),
    #                       tags$br(),
    #                       tags$h5("Max value"),
    #                       textInput("tab4.filter3.max", NULL),
    #                       tags$br(),
    #                       tags$h5("Min Value"),
    #                       textInput("tab4.filter3.min", NULL)
    #                     )),
    #
    #                     column(3, wellPanel(
    #                       tags$h4("Filter 4 (Numeric)"),
    #                       tags$br(),
    #                       tags$h5("select a column to filter on"),
    #                       selectInput("tab4.filter4.column", label="column select", choices=c("POS", "gt_GT", "...")),
    #                       tags$br(),
    #                       tags$h5("Max value"),
    #                       textInput("tab4.filter4.max", NULL),
    #                       tags$br(),
    #                       tags$h5("Min Value"),
    #                       textInput("tab4.filter4.min", NULL)
    #                     ))
    #                   )
    #          )
    # ),

## Tab 5 - Alignments #########################################################
    tabPanel("Alignments",
             tags$br(),
             tags$div(class="input-format",
                      tags$h3("Select Genes and Type"),
                      tags$h5("Select one or more transcript IDs below and the type of alignment to show"),
                      # textInput(inputId="tab3.Gene", label=NULL,
                      #           value="AT1G80490"),
                      uiOutput("tab5.selectGene")),
                      # actionButton(inputId="tab3.Submit", label="Submit"),
             tags$br(),
             tags$hr(),
             tags$div(class = "output-format",
                      tags$h3("Sequence Alignment"),
                      tags$h5("Click and drag to pan. The x-axis is the position within the alignment. Hover over the alignment to see details. 'seq_pos' is the position in the sequence with name 'seq_name' of the type chosen above. Use the pop-up menu in the upper right for zoom and other plotly functionalities. Made with",
                              tags$a(href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0749-z", target = "_blank", "DECIPHER")),
               plotlyOutput('tab5.aln_plot', height = "auto"),
              # verbatimTextOutput("event")
              tags$br(),
              tags$h5("Click and drag to pan. Made with", tags$a(href="https://zachcp.github.io/msaR/", "msaR")),
               msaROutput(outputId = "tab5.alignment"), height = "auto")
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
  )
  # "THIS IS THE FOOTER"
)}


# Server =================================================================

#source("VCF_Utils.R")
#source("Strains_and_Gene_Families.R")

parseInput <- function (textIn) {
  names <- str_extract_all(textIn, "AT[1-5]G[0-9]{5}")
  return (names[[1]])
}

# load_tab_2_Data <- function (geneInfo){
#   tab2VCF <- VCFByTranscript(geneInfo[1, ], strains)
#   tab2data <- tab2VCF$dat
#   tab2data <- parseEFF(tab2data)
#   tab2data <- Nucleotide_diversity(tab2data)
#
#   coding_variants <- coding_Diversity_Plot(tab2data)
#
#   return(coding_variants)
# }



# plotPi <- function(uniqueCodingVars) {
#   plot <- ggplot(uniqueCodingVars, aes(x=Codon_Number,y=Diversity, colour=Effect)) +
#     geom_point() +
#     scale_y_log10(breaks=c(0.0001, 0.001, 0.01, 0.1),limits=c(0.0001, 1)) +
#     #scale_colour_manual(values=c(synonymous_diversity="blue", missense_diversity="red")) +
#     ylab("nucleotide diversity, log scale")
#   return(plot)
#
# }




#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



server <- function(input, output){

  ##   _________
  ##  /  tab1   \
  ##             --------------------------------------------------
  ## Tab 1 ####################

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

  #### all.GeneChoices ####
  all.GeneChoices <- reactive({
    displayNames <- paste(all.Genes()$transcript_ID, " (", all.Genes()$tair_symbol, ")", sep="" )
    displayNames <- gsub(" \\(\\)", displayNames, replacement="")  # if no tair symbol, remove empty parens.
    output <- all.Genes()$transcript_ID
    names(output) <- displayNames
    return(output)
  })
#### tab1.genes_table ####
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
#### tab2.GeneInfo ####
  output$tab2.GeneInfo <- renderTable(tab2.Genes())
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
  output$diversityPlot <- renderPlot(plotCodingDiv(tab2.tableData()))
    #plot output
#### info ####
  output$info <- renderPrint({
    brushedPoints(tab2.tableData(), input$plot_brush, "Codon_Number", "Diversity")
  })


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


  ##                                        _________
  ##                                      /   tab5   \
  ## --------------------------------------           ----------------
  ## Tab 5 #########################
#### tab3.selectGene ####
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
  output$type <- eventReactive(input$tab5.Submit, {
    return(input$tab5.type)
  })
#### alignment ####
  alignment <- eventReactive(input$tab5.Submit, {
    alignment <- alignCDS(IDs = tab5.Genes())
    return(alignment)
  })
#### tab5.alignment ####
  output$tab5.alignment <- renderMsaR({
    type <- switch(EXPR = input$tab5.type, "DNA" = 0, "AA" = 1)
    msaR(alignment()[[type+1]], alignmentHeight = 100,
         colorscheme = {if(type) "taylor" else "nucleotide"})
    })
#### tab5.BrowseSeqs ####
  output$tab5.BrowseSeqs <- reactive({
    type <- switch(EXPR = input$tab5.type, "DNA" = 0, "AA" = 1)
    file <- BrowseSeqs(alignment()[[type + 1]],
                       openURL = FALSE)
    html <- paste(readLines(file), collapse="\n")
    return(html)
    })
#### aln_df ####
  aln_df <- reactive({
    type <- switch(EXPR = input$tab5.type, "DNA" = 0, "AA" = 1)
    aln_df <- makeAlnDF(alignment()[[type + 1]])
    vcf <- ldply(.data = all.VCFList()[input$tab5.transcript_ID],
                 .fun = subset, !is.na(Transcript_ID) & gt_GT != "0|0")
    vcf <- getCodingDiv(vcf)
    aln_df <- addSNPsToAlnDF(aln_df, vcf)
  })
#### tab5.aln_plot ####
  output$tab5.aln_plot <- renderPlotly({
    p <-
      ggplot(aln_df(), aes(x = aln_pos, y = seq_name,
                           group = seq_pos, text = variants)) +
      geom_tile(data = na.omit(aln_df()), mapping = aes(fill = effects),
                width = 1, height = 1, alpha = 0.8) +
      geom_text(aes(label=letter), alpha= 1,
                check_overlap = TRUE) +
      scale_x_continuous(breaks=seq(1,max(aln_df()$aln_pos), by = 10)) +
      # expand increases distance from axis
      xlab("") +
      ylab("") +
      #scale_size_manual(values=c(5, 6)) + # does nothing unless 'size' is mapped
      theme_logo(base_family = "Courier") +
      theme(panel.grid = element_blank(), panel.grid.minor = element_blank())
    ggplotly(p, tooltip = c("seq_name", "seq_pos", "variants"),
             height = length(unique(aln_df()$seq_name)) * 20 + 110) %>%
      # Total height = 125 = N*50 + 25
      # 125 = N*30 + 65 better still a bit long with 6 sequences
      # 125 = N*20 + 85 looks really good but maybe a bit
      config(collaborate = FALSE) %>%
      layout(margin = list(l = 100, r = 0, t = 20, b = 0),
             legend = list(yanchor = "bottom", y = -1, orientation = "h"),
             dragmode = "pan", yaxis = list(fixedrange = TRUE),
             xaxis = list(range = c(0,70)))
  })
}

enableBookmarking(store = "url")
shinyApp(ui = ui, server = server)


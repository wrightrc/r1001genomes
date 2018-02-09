library(leaflet)
library(shinyBS)
library(shinythemes)
library(knitr)
library(ggplot2)

ui <- function(request){ fluidPage(
  theme = shinytheme(theme = "flatly"),
  #### Header ####
  tags$head(
    # Custom CSS code
    tags$link(rel = "stylesheet", type = "text/css", href = "simple_blue_theme.css"),
    tags$link(rel = "stylesheet", type = "text/css", href = "general.css")
  ),
  titlePanel(
    tags$div(class= "title-panel",
             "Arabidopsis Natural Variation Webtool"
    )
  ),
  "This app provides an interface to examine the natural variation of specified genes of interest in the 1001 Genomes project dataset. To save or share a state of this app, use the bookmark button.", HTML("</br>"),
  bookmarkButton(),
  tags$h5('style'="color:red", "This app is currently a work in progress."),
  # themeSelector(),
  tags$br(),
  bsCollapse(id = "collapse 1", multiple=TRUE, open=c("Gene Select"),
             bsCollapsePanel(
               title=tagList("Gene Select", tags$h6(style= "display:inline", "(click to expand/collapse)")),
               value="Gene Select",
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
                                      tags$h5("browse to a .csv file containing 'tair_locus' and 'name' fields.
                                              The name field should be the TAIR symbol or moniker you would like to identify your genes by."),
                                      fileInput("genesFile", label=NULL),
                                      actionButton(inputId="file_submit", label = "Submit")

                                      )
                             )
             ),
             bsCollapsePanel(
               title=tagList("Annotation Files", tags$h6(style= "display:inline", "(click to expand/collapse)")),
               value="Annotation Files",
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
  tags$div(class="navbar-margin",navbarPage(title="TABS:",
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
               #tags$br(),
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
                        tags$div(style = "position:relative",
                          plotOutput("diversityPlot", brush="plot_brush", hover = hoverOpts("div_plot_hover", delay = 100,
                                                                                         delayType = "debounce"), height = 400),
                          uiOutput("tab2.hover")
                        ),
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
               #tags$br(),
               tags$div(class="input-format",
                        tags$h3("Select Genes and Filter Diversity Parameter"),
                        tags$h5("Select one or more transcript IDs below and use the slider to select a minimum sitewise nucleotide diversity"),
                        # textInput(inputId="tab3.Gene", label=NULL,
                        #           value="AT1G80490"),
                        uiOutput("tab3.selectGene"),
                        # actionButton(inputId="tab3.Submit", label="Submit"),
                        sliderInput(inputId="tab3.filter_value", label="Log Nucleotide diversity filter limit",
                                    min=-4, max=0, value=c(-3, -1), step=0.05),
                        radioButtons("tab3.SNPtype", "Type of SNP to mark",
                                     choices=c("All", "Coding", "Missense"))
                        # verbatimTextOutput("tab3.debug")
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
               #tags$br(),
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
               #tags$br(),
               tags$div(class="input-format",
                        tags$h3("Select Genes and Type"),
                        tags$h5("Select one or more transcript IDs below and the type of alignment to show"),
                        uiOutput("tab5.selectGene")
               ),
               tags$br(),
               tags$div(class = "output-format",
                        tags$h3("Sequence Alignment"),
                        tags$h5("Alignment made with",tags$a(href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0749-z", target = "_blank", "DECIPHER"), ". The x-axis is the position within the alignment. Hover over the alignment to see details (ggplot2 tooltip by", tags$a(href = "https://gitlab.com/snippets/16220", target = "_blank", "Pawel"), ". 'seq_pos' is the position in the sequence with name 'seq_name' of the type chosen above."),
                        tags$div(
                          style = "position:relative",
                          uiOutput("plot.ui"),
                          uiOutput("aln_plot_hover")),
                        plotOutput('tab5.aln_plot_legend', height = "200px")
                        ),
               # verbatimTextOutput("event")
               tags$br()#,
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
               #tags$br(),
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
  )) #end of tabset panel

  # "THIS IS THE FOOTER"
)}


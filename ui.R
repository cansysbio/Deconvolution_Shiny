library(shiny)
library(shinyjs)
library(shinyBS)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(plotly)
library(htmlwidgets)
library(BBmisc)
library(GSVA)
library(singscore)
library(plyr)
library(shinythemes)
library(rdrop2)
library(readr)
library(shinycssloaders)
library(parallel)
library(promises)
library(future)
library(ipc)
library(shinyWidgets)
library(future.apply)
library(heatmaply)
library(ggbeeswarm)
library(data.table)
library(crosstalk)


plan(multiprocess)

cancerTypes = c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC",
                "KICH", "KIRC", "KIRP","LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD",
                "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", 
                "UCEC", "UCS", "UVM")

cancerList = as.list(cancerTypes)
radioTooltip <- function(id, choice, title, placement = "bottom", trigger = "hover", options = NULL){
  
  options = shinyBS:::buildTooltipOrPopoverOptionsList(title, placement, trigger, options)
  options = paste0("{'", paste(names(options), options, sep = "': '", collapse = "', '"), "'}")
  bsTag <- shiny::tags$script(shiny::HTML(paste0("
                                                 $(document).ready(function() {
                                                 setTimeout(function() {
                                                 $('input', $('#", id, "')).each(function(){
                                                 if(this.getAttribute('value') == '", choice, "') {
                                                 opts = $.extend(", options, ", {html: true});
                                                 $(this.parentElement).tooltip('destroy');
                                                 $(this.parentElement).tooltip(opts);
                                                 }
                                                 })
                                                 }, 500)
                                                 });
                                                 ")))
  htmltools::attachDependencies(bsTag, shinyBS:::shinyBSDep)
}


# UI for Application for Benchmarking Gene Sets
navbarPage("TME Cell Estimation Benchmarking", collapsible = TRUE, theme = "new.css", inverse = T,
           
  #### Home Page ####           
           tabPanel("Home", icon = icon("fas fa-home"),
                    HTML("<h2> Benchmarking of Tumour Microenvrionment Cell Type Estimation From Bulk RNA </h2>" ),
                    sidebarPanel(
                      HTML("<h3> Immune Estimation / Deconvolution Tools </h3>" ),
                      tags$hr(),
                      HTML(paste0("The tumour microenvironment (TME) is comprised of complex cellular ",
                                  "compositions with the interactions between the cancer, immune ",
                                  "and stromal components playing a crucial role in determining tumour ",
                                  "progression, metastasis and response to treatment. Thus, being able ",
                                  "to quantify the presence of different cell types is of paramount importance. ",
                                  "The unprecedented availability of RNA sequencing data has motivated the ",
                                  "development of computational tools to predict the presence of immune cells ",
                                  "from bulk transcriptomic profiles. ")),
                      tags$br(),
                      tags$br(),
                      HTML("<center><img src='Deconvolution_Methods.png' width='75%' height='auto' ></center>"),
                      tags$br(),
                      HTML(paste0("The majority of the tools tend to fall into two categories to solove this ",
                                  "problem: <br> <br> <b>1) Gene Set Scoring Based Methods</b> - For these methods ",
                                  "presence of immune cells is represented by specific gene sets then different ",
                                  "statistical frameworks can be leveraged to generate enrichment scores. <br> <br> <b>2) ",
                                  "Regression Based Deconvolution Methods</b> - These methods use an additional ",
                                  "signature matrix containing gene expression profiles of immune cells. Then, using ",
                                  "various regression methods, attempt to resolve to what degree the gene expression ",
                                  "profiles of each immune cell are present in the bulk tumour profile. <br><br> A ",
                                  "thorough review of the different approaches currently used was carried out by ",
                                  "<a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6006237/'>Finotello & Trajanoski (2018) </a>.")),
                      width = 12
                      ),
                    sidebarPanel(
                      HTML("<h3> Benchmarking Approach </h3>" ),
                      tags$hr(),
                      HTML(paste0("An issue faced when comparing different tools is that performance of a tool ",
                                  "in one biological context may not be seen in another. Subsequently here we aim ",
                                  "to provide a series of un-baised and objective benchmarks. These fall into two ",
                                  "categories: <br><br> <b> TCGA Benchmarks </b> - This set of benchmarks leverages",
                                  "the large scale pan-cancer RNA data. With no ground truth immune estimation ",
                                  "available, a series of orthogonal inferences were used. <br><br> <b> Cell Type ",
                                  "Specific Benchmarks </b> - These benchmarks allow cell type specific comparisons using ",
                                  "various \"ground truth estimates\". However, these represent only a few anatomical contexts ",
                                  "and often have small sample sizes.")),
                      tags$br(),
                      HTML("<center><img src='Benchmarking_Overview.png' width='75%' height='auto' ></center>"),
                      width = 12)
                    ),
  #### Benchmark New Genesets ####
           navbarMenu("Benchmark New Geneset",
                      tabPanel("Input Gene Signatures",
                               fluidPage(
                                 useShinyjs(),
                                 sidebarLayout(
                                   sidebarPanel(
                                     selectInput(
                                       "consCanIn",
                                       HTML("<h4> Use Consensus<sup>TME</sup> Gene List \n for Cancer Type: </h4>"),
                                       choices = cancerList),
                                     actionButton("loadCons",
                                                  label = HTML("Load Consensus<sup>TME</sup> Geneset")),
                                     tags$hr(),
                                     h4("Use New Gene Set \n "),
                                     fileInput("geneFile", "Choose Gene Signatures File:",
                                               accept = c(
                                                 "text/csv",
                                                 "text/comma-separated-values,text/plain",
                                                 ".csv",
                                                 ".rds"
                                                 )
                                               ),
                                     radioButtons(inputId = "fileType", "Gene Set Type",
                                                  choices = c("Single Set of Gene Signatures - Dataframe" = "fileSingle",
                                                              "Cancer Specific Signatures List - List (.rds)" = "fileList"
                                                              ),
                                                  selected = "fileSingle"
                                                  ),
                                     radioTooltip(id = "fileType",
                                                  title =  "Upload .txt or .csv ",
                                                  placement = "bottom", trigger = "hover",
                                                  choice = "fileSingle",
                                                  options = NULL),
                                     radioTooltip(id = "fileType",
                                                  title =  "Upload .rds File",
                                                  placement = "bottom", trigger = "hover",
                                                  choice = "fileList",
                                                  options = NULL),
                                     tags$hr(),
                                     uiOutput("inSettings"),
                                     tags$hr()
                                     ),
                                   mainPanel(
                                     uiOutput("previews")
                                     )
                                   )
                                 )
                               ),
                      tabPanel("Produce Cell Type Estimates",
                               fluidPage(
                                 useShinyjs(),
                                 useSweetAlert(),
                                 div(id ="estSidebar",
                                     sidebarPanel(
                                       HTML("<center> <h3> Produce Estimates For Benchmarking Data Sets </h3> </center>"),
                                       tags$hr(),
                                       div(radioButtons("statMethod",
                                                        "Choose Statistical Framework:",
                                                        choices = c("ssGSEA - single sample Gene Set Enrichment Analysis" = "statSsgsea",
                                                                     "singScore - alternative single sample Scoring approach" = "statSingscore"),
                                                        selected = "statSsgsea"
                                                        ),
                                           style="font-size:100%"
                                           ),
                                       tags$hr(),
                                       uiOutput("isSelector"),
                                       tags$hr(),
                                       textInput(
                                         "geneSetName",
                                         "Enter Name For New Gene Set",
                                         placeholder = "e.g: MyGeneSet"
                                         ),
                                       fluidRow(
                                         column(4, offset = 4,
                                                actionBttn("runEstimation",
                                                           "Generate Estimates For Benchmarking Datasets",
                                                           style = "gradient",
                                                           color = "primary",
                                                           block = TRUE,
                                                           size = "md"
                                                           )
                                                )
                                         ),
                                       width = 12
                                       )
                                     ),
                                 div(id ="progressPanel",
                                     sidebarPanel(
                                       HTML("<center> <h3>  Progress  </h3> </center>"),
                                       tags$hr(),
                                       fluidRow(
                                         column(
                                           htmlOutput("progressToDo"),
                                           width = 6
                                           ),
                                         column(
                                           htmlOutput("progressFinished"),
                                           width = 6
                                         )
                                        ),
                                       tags$hr(),
                                       progressBar(id = "estProgressBar",
                                                   value = 0,
                                                   display_pct = TRUE,
                                                   striped = T,
                                                   total = length(cancerTypes),
                                                   title = "Estimation in Progress"),
                                       width = 12
                                     )
                                 ),
                                 div(id ="downloadPanel",
                                     sidebarPanel(
                                       HTML("<center> <h3>  Download Generated Estimates  </h3> </center>"),
                                       tags$br(),
                                       div(id = "estRdsDownButton",
                                           downloadButton(
                                             outputId = "estRdsDownload",
                                             label = HTML("<h4 style='color:#F0F1F2'> Download Single .rds File </h4>"),
                                             ),
                                           style = "text-align: center"
                                           ),
                                       bsTooltip("estRdsDownButton",
                                                 "For Use Within The App",
                                                 placement = "bottom",
                                                 trigger = "hover",
                                                 options = NULL),
                                       tags$br(),
                                       div(id = "estTxtDownButton",
                                           downloadButton(
                                             outputId = "estTxtDownload",
                                             label = HTML("<h4 style='color:#F0F1F2'> Download Multiple .txt Files <br> (1 Per Dataset - Zipped) </h4>"),
                                             ),
                                           style = "text-align: center"
                                           ),
                                       tags$br(),
                                       bsTooltip("estTxtDownButton",
                                                 "For Other Uses",
                                                 placement = "bottom",
                                                 trigger = "hover",
                                                 options = NULL),
                                     width = 12
                                     )
                                 ),
                                 div(id ="loadEstSidebar",
                                   sidebarPanel(
                                   HTML("<center> <h3> Load Previously Produced Estimates For Benchmarking Data Sets </h3> </center>"),
                                   tags$hr(),
                                   fileInput("estFile", "Upload Previously Generated .rds File:",
                                             accept = ".rds",
                                             width = "50%"
                                             ),
                                   width = 12
                                   )
                                 )
                               )
                      ),
    tabPanel("Explore Estimates",
             mainPanel(
               uiOutput("estPreviews"),
               width = 12
               )
             )
           ),

  #### Tumour Purity Benchmark #### 
  navbarMenu("TCGA Benchmarking",
             tabPanel("Tumour Purity",
                      HTML("<h2> Tumour Purity Benchmark </h2>" ),
                      tags$hr(),
                      HTML(paste0("<font color=\"#101011\"> This benchmark assesses the ability ",
                                  "of the immune estimation tools to capture the global level of ",
                                  "immune infiltration into the tumour. This is assessed through the " ,
                                  "use of copy number and mutation data available from TCGA to derive tumour purity. ",
                                  "Immune estimation tools with a good ability to estimate immune infiltation ",
                                  "show a strong negative correlation with tumour purity. </font>")),
                      tags$br(),
                      tags$br(),
                      withSpinner(
                        plotlyOutput("purCorrPlotly"),
                        type = 6
                        ),
                      tags$br(),
                      tags$br(),
                      withSpinner(
                        plotlyOutput("purCorrBoxPlotly"),
                        type = 6
                        )
                      ),
  #### Leukocyte Methylation Benchmark ####
             tabPanel("Leukocyte Methylation Benchmark",
                      HTML("<h2> Leukocyte Methylation Benchmark </h2>" ),
                      tags$hr(),
                      HTML(paste0("<font color=\"#101011\"> This benchmark assesses the ability ",
                                  "of the immune estimation tools to accurately predict leukocyte ",
                                  "infiltration in the tumour. This is assessed through the use of " ,
                                  "methylation data from TCGA samples. Previous work (Hoadley et al. Cell 2008) ",
                                  "identified loci across the epigenome with differential methylation in immune ",
                                  "vs tumour populations and used this information to estimation a leukocyte ",
                                  "fraction for each sample within TCGA. We assess method accuracy using the ",
                                  "assumption that variation in leukocyte fraction should be able to be explained ",
                                  "by variation in the cell type estimates in the category of being leukocytes. ",
                                  "Multiple linear regression is performed for each cancer and various goodness of fit ",
                                  "metrics are used to assess performance. </font>")),
                      tags$br(),
                      tags$br(),
                      fluidRow(
                        column(
                          HTML("<center><img src='Leukocyte_Methylation.png' width='100%' height='auto' ></center>"),
                          width = 6
                          ),
                        column(
                          HTML("<center><img src='MLR_Equation.png' width='100%' height='auto' ></center>"),
                          width = 6
                        )
                      ),
                      tags$br(),
                      tags$br(),
                      HTML("<h4> Fit Metric:  </h4>"),
                      tabsetPanel(
                          tabPanel("Adjusted R-Squared",
                                   tags$br(),
                                   HTML("<h4> Adjusted R-Squared: Higher Value = Better Model </h4>"),
                                   tags$br(),
                                   withSpinner(
                                     plotlyOutput("leukR2FitHeatmap"),
                                     type = 6
                                   ),
                                   tags$br(),
                                   tags$br(),
                                   withSpinner(
                                     plotlyOutput("leukR2FitBoxplot"),
                                     type = 6
                                   )
                          ),
                          tabPanel("AIC",
                                   tags$br(),
                                   HTML("<h4> Akaike Information Criterion: Lower Value = Better Model </h4>"),
                                   tags$br(),
                                   withSpinner(
                                     plotlyOutput("leukAicFitHeatmap"),
                                     type = 6
                                   ),
                                   tags$br(),
                                   tags$br(),
                                   withSpinner(
                                     plotlyOutput("leukAicFitBoxplot"),
                                     type = 6
                                   )
                          ),
                          tabPanel("BIC",
                                   tags$br(),
                                   HTML("<h4> Bayesian Information Criterion: Lower Value = Better Model </h4>"),
                                   tags$br(),
                                   withSpinner(
                                     plotlyOutput("leukBicFitHeatmap"),
                                     type = 6
                                   ),
                                   tags$br(),
                                   tags$br(),
                                   withSpinner(
                                     plotlyOutput("leukBicFitBoxplot"),
                                     type = 6
                                   )
                          )
                      )
             ),
             tabPanel("Image Analysis",
                      HTML("<h2> Image Analysis Benchmark </h2>" ),
                      tags$hr(),
                      HTML(paste0("<font color=\"#101011\"> This benchmark assesses the ability ",
                                  "of the immune estimation tools to accurately predict lymphocyte ",
                                  "infiltration in the tumour. This is assessed through the use of " ,
                                  "H&E pathology slides availble from TCHA. Previous work by the Thorsson group ",
                                  "(Saltz et al. Cell Reports 2018) used a convolutional neural network to generate ",
                                  "tumour infiltration lymphocyte estimates by image analysis. In this benchmark we ",
                                  "use the lymphocyte score as the response variable for multiple linear regression ",
                                  "with each of cell type estimates fitting into the category of being lymphocytes as ",
                                  "explanatory variables. Models were assessed using three goodness of fit metrics to ",
                                  "account for varying numbers of terms (i.e cell types) in the models. A smaller number ",
                                  "of cancer types have H&E slides available. </font>")),
                      tags$br(),
                      tags$br(),
                      HTML("<center><img src='Image_Analysis_MLR.png' width='80%' height='auto' ></center>"),
                      tags$br(),
                      tags$br(),
                      HTML("<h4> Fit Metric:  </h4>"),
                      tabsetPanel(
                        tabPanel("Adjusted R-Squared",
                                 tags$br(),
                                 HTML("<h4> Adjusted R-Squared: Higher Value = Better Model </h4>"),
                                 tags$br(),
                                 withSpinner(
                                   plotlyOutput("imageR2FitHeatmap"),
                                   type = 6
                                 ),
                                 tags$br(),
                                 tags$br(),
                                 withSpinner(
                                   plotlyOutput("imageR2FitBoxplot"),
                                   type = 6
                                 )
                        ),
                        tabPanel("AIC",
                                 tags$br(),
                                 HTML("<h4> Akaike Information Criterion: Lower Value = Better Model </h4>"),
                                 tags$br(),
                                 withSpinner(
                                   plotlyOutput("imageAicFitHeatmap"),
                                   type = 6
                                 ),
                                 tags$br(),
                                 tags$br(),
                                 withSpinner(
                                   plotlyOutput("imageAicFitBoxplot"),
                                   type = 6
                                 )
                        ),
                        tabPanel("BIC",
                                 tags$br(),
                                 HTML("<h4> Bayesian Information Criterion: Lower Value = Better Model </h4>"),
                                 tags$br(),
                                 withSpinner(
                                   plotlyOutput("imageBicFitHeatmap"),
                                   type = 6
                                 ),
                                 tags$br(),
                                 tags$br(),
                                 withSpinner(
                                   plotlyOutput("imageBicFitBoxplot"),
                                   type = 6
                                 )
                        )
                      )
             )
  ),
  #### Cell-Type Specific Benchmarking ####
  navbarMenu("Cell-Type Specific Benchmarking",
             tabPanel("MCP-Counter Dataset",
                      HTML("<h2> MCP-Counter Colon Cancer IHC Benchmark </h2>" ),
                      tags$hr(),
                      HTML(paste0("<font color=\"#101011\"> This benchmark replicates a validation experiment ",
                                  "carried out in the orginal MCP-Counter manuscript (Becht et al. Genome Biology 2016). ",
                                  "Here the authors used a collection of colon cancer samples with matching RNA and ",
                                  "immunohistochemistry (IHC) staining. Performance of methods was assessed by ",
                                  "generating estimates for each of the methods from RNA and correlating these against ",
                                  "estimates from IHC. Due to methods estimating different cell types, where appropriate ",
                                  "cell sub-types are combined together in order to match IHC marker categories. ",
                                  "While only assessing the performance of three cell-types the advantages of this benchmark ",
                                  "is to samples are <i>in-vivo</i> tumour samples recapitulating the transcriptome in ",
                                  "which the immune estimation tools are intended to be used. </font>")),
                      tags$br(),
                      tags$br(),
                      withSpinner(
                        plotlyOutput("mcpScatter",
                                     height = 600,
                                     width = "100%"),
                        type = 6
                      ),
                      tags$hr(),
                      withSpinner(
                        plotlyOutput("mcpBar"),
                        type = 6
                      )
                      ),
             tabPanel("TIMER Dataset",
                      HTML("<h2> TIMER Bladder Carcinoma Pathologist Estimation Benchmark </h2>" ),
                      tags$hr(),
                      HTML(paste0("<font color=\"#101011\"> This benchmark replicates a validation experiment ",
                                  "carried out in the orginal TIMER manuscript (Li et al. Genome Biology 2016). ",
                                  "Here the authors used the heamatoxylin and eosin (H&E) slides from the Bladder ",
                                  "Carcinoma study (BLCA) in TCGA. A pathologist reviewed 404 slides and classified ",
                                  "each slide into having \"Low\" \"Medium\" or \"High\" levels of neutrophil abundance. ",
                                  "Performance was then assessed by generating neutrophil estimates for all corresponding ",
                                  "samples from RNA and comparing values from each category. Analysis of variance (ANOVA) ",
                                  "with Tukey Honest Significant Difference (HSD) post-hoc tests were used to assess ",
                                  "differences between RNA derived scores across pathologist estimation categories. </font>")),
                      tags$br(),
                      tags$br(),
                      withSpinner(
                        plotlyOutput("timerBoxplots",
                                     height = 600,
                                     width = "100%"),
                        type = 6
                      )
                      ),
             tabPanel("xCell Dataset",
                      HTML("<h2> xCell PBMC Benchmark </h2>" ),
                      tags$hr(),
                      HTML(paste0("<font color=\"#101011\"> This benchmark replicates validation experiments ",
                                  "carried out in the orginal xCell manuscript (Aran et al. Genome Biology 2017). ",
                                  "Here the authors used data collected from ImmPort, accession SDY311 & SDY420; ",
                                  "these contained PBMC samples from 61 and 104 healthly individuals respectively.  ",
                                  "Performance of immune estimation methods was assessed using matching RNA-Seq and ",
                                  "CyTOF data for each of the samples. Due to methods estimating different cell types",
                                  ", where appropriate cell sub-types are combined together in order to match CyTOF ",
                                  "categories. While having the benefits of giving an estimate of cell type specific performance ",
                                  "this form of benchmark suffers from using immune cells from peripheral blood which doesn't ",
                                  "reflect the complexity of the tumour microenvironment. </font>")),
                      tags$br(),
                      tags$br(),
                      withSpinner(
                        plotlyOutput("xCellBoxplots420",
                                     height = 450,
                                     width = "100%"),
                        type = 6
                      ),
                      tags$hr(),
                      withSpinner(
                        plotlyOutput("xCellBoxplots311",
                                     height = 450,
                                     width = "100%"),
                        type = 6
                      )
                      ),
             tabPanel("CIBERSORT Dataset",
                      HTML("<h2> CIBERSORT PBMC Benchmark </h2>" ),
                      tags$hr(),
                      HTML(paste0("<font color=\"#101011\"> This benchmark replicates a validation experiment ",
                                  "carried out in the orginal CIBERSORT manuscript (Newman et al. Nature Methods 2015). ",
                                  "Here the authors used a collection of peripheral blood mononuclear cells (PBMCs) isolated ",
                                  "from 20 adults receiving influenza immunization. Performance of methods was assessed by ",
                                  "generating estimates for each of the methods from RNA and correlated against flow ",
                                  "cytometry fractions. Due to methods estimating different cell types, where appropriate ",
                                  "cell sub-types are combined together in order to match flow cytometry categories. ",
                                  "While having the benefits of giving an estimate of cell type specific performance ",
                                  "this form of benchmark suffers from <b> using immune cells from peripheral blood of ",
                                  "healthly individuals </b> which doesn't reflect the complexity of the tumour microenvironment. </font>")),
                      tags$br(),
                      tags$br(),
                      withSpinner(
                        plotlyOutput("cibersortBoxplots",
                                     height = 500,
                                     width = "100%"),
                        type = 6
                        )
                      ),
             tabPanel("HGSOC Dataset",
                      HTML("<h2> High Grade Serous Ovarian Cancer Benchmark </h2>" ),
                      tags$hr(),
                      HTML(paste0("<font color=\"#101011\"> This benchmark replicates a validation experiment ",
                                  "carried out in an independent study involving Consensus<sup>TME</sup> ",
                                  "(Jimenez-Sanchez et al. 2018, Under Revision, ",
                                  "<a href='https://www.biorxiv.org/content/10.1101/441428v2.full'>",
                                  " bioRxiv Preprint available </a>). Tumours from patients with high-grade ",
                                  " serous ovarian cancer (HGSOC) were used. Methods were benchmarked by ",
                                  " correlating bulk tumour mRNA-based immune estimates against IHC counts for ",
                                  "CD4<sup>+</sup> T Cells, CD8<sup>+</sup> T Cells and T Regulatory Cells. ",
                                  "<b> N.B. </b> This dataset was used in the orginal development of Consensus",
                                  "<sup>TME</sup></font>")),
                      tags$br(),
                      tags$br(),
                      HTML("<center><img src='HGSOC_Graphical_Abstract.png' width='70%' height='auto' ></center>"),
                      tags$br(),
                      tags$br(),
                      withSpinner(
                        plotlyOutput("hgsocBar"),
                        type = 6
                      )
             )
             ),
  #### Results Overview ####
  tabPanel("Results Overview",
           HTML("<h2> Results Overview </h2>" ),
           tags$hr(),
           HTML(paste0("<font color=\"#101011\"> Here the overall performance of methods across all benchmarks ",
                       "can be visualised. Benchmarking experiments fall into three broad categories: <b>1) TCGA</b> - ",
                       "Benchmarks carried out using orthogonal inferences from TCGA data. <b>2) PBMCs</b> - Benchmarks ",
                       "carried out in which authours used peripheral blood mononuclear cells (PBMCs) derived from ",
                       "circulation. <b>3) Bulk Tumour</b> - Benchmarks using samples derived from the setting of a bulk ",
                       "tumour. The mean rank of each method is also plotted </font>")),
           tags$br(),
           tags$br(),
           withSpinner(
             plotlyOutput("overviewLinePlot",
                          height = 500,
                          width = "100%"),
             type = 6
           )
  ),
  #### About ####
  tabPanel("About",
           wellPanel(
             HTML("<h3> Overview </h3>"),
             tags$br(),
             fluidRow(
               column(
                 HTML(paste0("This online portal serves as a companion to the the manuscript:<br><br>",
                             "<a href='https://www.biorxiv.org/content/10.1101/437533v2.full'>",
                             "\"Comprehensive Benchmarking and Integration of Tumour Microenvironment ",
                             "Cell Estimation Methods </a><br><br>",
                             "This serves a dual purpose of allowing interactive exploration of large, ",
                             "multi-dimensional data but also allows benchmarking results to be evolvable; ",
                             "as both new methods and new benchmarking datasets become available the portal ",
                             "can be updated to ensure orginal benchmarks are more than a snapshot in time of ",
                             "method performance."
                             )),
                 width = 11,
                 offset = 1
                 )
             ),
             tags$hr(),
             tags$br(),
             HTML("<h3> Consensus<sup>TME</sup> Approach </h3>"),
             tags$br(),
             tags$br(),
             fluidRow(
               column(
                 HTML(paste0("The approach taken by Consensus<sup>TME</sup> is fully described in ",
                             "<a href='https://www.biorxiv.org/content/10.1101/437533v2.full'>",
                             "the manuscript</a>. However, to summerise briefly, the approach consists of ",
                             "X main steps: <br><br><b>1)</b> Multiple sources are used for initial compiling of ",
                             "gene sets. Either from pre-existing carefully curated gene sets (e.g ",
                             "<a href='https://jitc.biomedcentral.com/articles/10.1186/s40425-017-0215-8'>",
                             "Danaher gene sets </a>) or through analysis of signature matricies leveraged ",
                             "by other methods (e.g LM22 matrix used by ",
                             "<a href='https://www.nature.com/articles/nmeth.3337'>",
                             "CIBERSORT </a>).<br><br><b>2)</b> Defining cell types for which gene sets can be ",
                             "derived from at least two methods.<br><br><b>3)</b> Create a unique union of genes from ",
                             "multiple sources. <br><br><b>4)</b> Using an approach, orginally used by the ",
                             "<a href='https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1028-7'>",
                             "TIMER </a> tool, we ensure that only genes whose expression has a negative correlation ",
                             "with tumour purity are included. This is done on a cancer by cancer basis for each of ",
                             "the cancer types available from TCGA and increases the confidence that if the expression ",
                             "of a gene is increasing it is due to the presence of an immune cell instead of spurious ",
                             "up-regulation by cancer cells.<br><br><b>5)</b> Use the gene sets with a statistical ",
                             "framework to produce normalised enrichment scores (NESs). Currently single sample gene set ",
                             "enrichment analysis (ssGSEA) is employed and benchmarked. Future versions of ",
                             "Consensus<sup>TME</sup> may use different approaches as they become available. ",
                             "<br><br><b>6)</b> NES output can be used to identify differences between immune cell ",
                             "subtype abundance between patients.<br><br> Consensus<sup>TME</sup> is available as a ",
                             "GitHub downloadable package: <a href='https://github.com/cansysbio/ConsensusTME'>",
                             "https://github.com/cansysbio/ConsensusTME</a>"
                             )
                      ),
                 width = 6, offset = 1
                 ),
               column(
                 HTML("<center><img src='ConsensusTME_Approach.png' width='50%' height='auto' ></center>"),
                 width = 5
               )
               ),
             tags$hr(),
             HTML("<h3>Contact</h3>"),
             tags$br(),
             fluidRow(
               column(
                 HTML(paste0("Alejandro Jim&eacute;nez-S&aacute;nchez<sup>&#8224;</sup>: Alejandro.JimenezSanchez@cruk.cam.ac.uk<br><br>",
                           "Oliver Cast: Oliver.Cast@cruk.cam.ac.uk<br><br>",
                           "Martin Miller<sup>&#8224;</sup>: Martin.Miller@cruk.cam.ac.uk <br>",
                           "<br>&#8224; = Corresponding Author <br><br>",
                           "For suggestions & queries regarding the portal:<br><br> devonvolutionbenchmarks@gmail.com")
                      ),
                 width = 7),
               column(
                 HTML("<center><img src='CI_Logo.png' width='70%' height='auto' ></center>"),
                 width = 5
                 )
               )
             )
           ),
  
  #### Help ####
  tabPanel("FAQ")
  #### Method Comparison ####
  #tabPanel("Method Comparison")
)

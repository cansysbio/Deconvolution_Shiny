library(shiny)
library(bslib)
library(shinycssloaders)
library(plotly)
library(DT)
library(crosstalk)
library(heatmaply)
library(RColorBrewer)
library(ggbeeswarm)
library(data.table)
library(markdown)

# UI for Application for Benchmarking Gene Sets
navbarPage(HTML("Consensus<sup>TME</sup>"),
           collapsible = TRUE,
           theme = bs_theme(fg = "rgb(14, 54, 107)",
                            font_scale = NULL,
                            `enable-rounded` = TRUE,
                            bootswatch = "lux",
                            bg = "rgb(255, 255, 255)"),
           windowTitle = "ConsensusTME",
           header = tags$head(
             includeCSS("./www/DT.css")
           ),
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
                                  "<a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6006237/' target='_blank'>Finotello & Trajanoski (2018) </a>.")),
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
           #### Cell type estimation ####
           tabPanel("Cell Type Estimation", icon = icon("fas fa-chart-simple"),
                    HTML("<h2> Running Consensus<sup>TME</sup> </h2>" ),
                    sidebarPanel(
                      HTML("<h3> Within the R environemnt </h3>" ),
                      tags$hr(),
                      HTML(paste0("Currently Consensus<sup>TME</sup> can only be run by installing the R ",
                                  "package. Instructions for installation are on the github <a href='https://github.com/cansysbio/ConsensusTME' target='_blank'>",
                                  "https://github.com/cansysbio/ConsensusTME</a> or below:")),
                      tags$br(),
                      tags$br(),
                      includeMarkdown("www/cons_git_README.md"),
                      width = 12
                    )
           ),
           #### Benchmarking Experiments ####  
           navbarMenu("Benchmarking Experiments", icon = icon("fas fa-chart-line"),
                      "TCGA Benchmarks:",
                      #### Tumour Purity Benchmark ####
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
                      ),
                      "----",
                      #### Cell-Type Specific Benchmarking ####
                      
                      "Cell-Type Specific Benchmarks:",
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
                                           "(Jimenez-Sanchez et al. 2020, Nature Genetics, ",
                                           "<a href='https://rdcu.be/dcsJz' target='_blank'>",
                                           " free access article </a>). Tumours from patients with high-grade ",
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
                      ),
                      "----",
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
                      
           ),
           
           #### Literature Resources ####
           navbarMenu("Literature Resources", icon = icon("fas fa-book"),
                      tabPanel("Current Approaches",
                               HTML("<h2> Overview of approaches & tools for cell type estimation </h2>"),
                               tags$p("To suggest additional tools to be added to this list please contact ", 
                                      tags$a(href = "mailto:Oliver.Cast@cruk.cam.ac.uk", "Oliver.Cast@cruk.cam.ac.uk")),
                               sidebarPanel(
                                 DTOutput("currAppTab"),
                                 width = 12
                               )
                      ),
                      "----",
                      tabPanel("Benchmarking Datasets",
                               HTML("<h2> Datasets for benchmarking & validation of cell estimation tools </h2>"),
                               tags$p("To suggest additional datasets to be added to this list please contact ", 
                                      tags$a(href = "mailto:Oliver.Cast@cruk.cam.ac.uk", "Oliver.Cast@cruk.cam.ac.uk")),
                               sidebarPanel(
                                 DTOutput("benchDsTab"),
                                 width = 12
                               )
                      ),
                      "----",
                      tabPanel("Review Articles",
                               HTML("<h2> Cell type estimation review and benchmarking studies </h2>"),
                               tags$p("To suggest additional articles to be added to this list please contact ", 
                                      tags$a(href = "mailto:Oliver.Cast@cruk.cam.ac.uk", "Oliver.Cast@cruk.cam.ac.uk")),
                               sidebarPanel(
                                 DTOutput("reviewTab"),
                                 width = 12
                               )
                      )
           ),
           
           #### About ####
           tabPanel("About",
                    wellPanel(
                      HTML("<h3> Overview </h3>"),
                      tags$br(),
                      fluidRow(
                        column(
                          HTML(paste0("This online portal serves as a companion to the the manuscript published in Cancer Research:<br><br>",
                                      "<a href='https://aacrjournals.org/cancerres/article/79/24/6238/639705/Comprehensive-Benchmarking-and-Integration-of' target='_blank'>",
                                      "\"Comprehensive Benchmarking and Integration of Tumour Microenvironment ",
                                      "Cell Estimation Methods\" </a><br><br>",
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
                                      "<a href='https://aacrjournals.org/cancerres/article/79/24/6238/639705/Comprehensive-Benchmarking-and-Integration-of' target='_blank'>",
                                      "the manuscript</a>. However, to summerise briefly, the approach consists of ",
                                      "6 main steps: <br><br><b>1)</b> Multiple sources are used for initial compiling of ",
                                      "gene sets. Either from pre-existing carefully curated gene sets (e.g ",
                                      "<a href='https://jitc.biomedcentral.com/articles/10.1186/s40425-017-0215-8' target='_blank'>",
                                      "Danaher gene sets </a>) or through analysis of signature matricies leveraged ",
                                      "by other methods (e.g LM22 matrix used by ",
                                      "<a href='https://www.nature.com/articles/nmeth.3337' target='_blank'>",
                                      "CIBERSORT </a>).<br><br><b>2)</b> Defining cell types for which gene sets can be ",
                                      "derived from at least two methods.<br><br><b>3)</b> Create a unique union of genes from ",
                                      "multiple sources. <br><br><b>4)</b> Using an approach, orginally used by the ",
                                      "<a href='https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1028-7' target='_blank'>",
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
                                      "GitHub downloadable R package: <a href='https://github.com/cansysbio/ConsensusTME' target='_blank'>",
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
                          HTML(paste0("Oliver Cast: Oliver.Cast@cruk.cam.ac.uk<br><br>",
                                      "Alejandro Jim&eacute;nez-S&aacute;nchez: ajs.scientia@gmail.com<br><br>",
                                      "Martin Miller<sup>&#8224;</sup><br>",
                                      "For suggestions & queries regarding the portal:<br><br> Oliver.Cast@cruk.cam.ac.uk")
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
)


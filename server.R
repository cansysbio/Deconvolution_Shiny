source('./useful_functions.R')
# cancerTypes <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC",
#                         "KICH", "KIRC", "KIRP","LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD",
#                         "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM",
#                         "UCEC", "UCS", "UVM")
cancerTypes <- c("UCEC", "UCS", "UVM")

cancerList <- as.list(cancerTypes)

estProgressBox <- 0

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

heatPrevButton <- tags$div(actionButton("heatPrevPlot",
                                         HTML('<div class="col-sm-4"><i class="fa fa-arrow-alt-circle-left"></i></div>')))
heatNextButton <- div(actionButton("heatNextPlot",
                             HTML('<div class="col-sm-4"><i class="fa fa-arrow-alt-circle-right"></i></div>')))
corrPrevButton <- tags$div(actionButton("corrPrevPlot",
                                        HTML('<div class="col-sm-4"><i class="fa fa-arrow-alt-circle-left"></i></div>')))
corrNextButton <- div(actionButton("corrNextPlot",
                                   HTML('<div class="col-sm-4"><i class="fa fa-arrow-alt-circle-right"></i></div>')))

function(input, output, session) {
  
  #### Load Required Data ####
  
  load("./data/Tumour_Purity_Corrs.RData")
  load("./data/TCGA_Annotation.RData")
  
  #### Input Gene Signatures ####

  # Check File Type
  inType <- eventReactive(input$geneFile, {
    if(!is.null(input$geneFile)){
      genesetFile <- input$geneFile
      extension <- strsplit(genesetFile$datapath, split = ".", fixed = TRUE)[[1]]
      if(extension[length(extension)] == "rds"){
        return("fileList")
      } else if(extension[length(extension)] %in% c("csv","txt")){
        return("fileSingle")
      } else {
        return(NULL)
      }
    } else {
      return(NULL)
    }
  })

  observeEvent(input$geneFile, {
    if(!is.null(input$geneFile)){
      geneset.file <- input$geneFile
      extension <- strsplit(geneset.file$datapath, split = ".", fixed = TRUE)[[1]]
      if(!extension[length(extension)] %in% c("csv","txt","rds")){
        sendSweetAlert(
          session = session,
          title = "Input File",
          text = "Unexpected File Extension - Please Use: \".txt\", \".csv\" or \".rds\"",
          type = "warning"
        )
      }
    }
  })

  # Update Radio Buttons and Add Tool Tip After Loading File

  observeEvent(input$geneFile, {
    updateRadioButtons(session, "fileType",
                       "Gene Set Type",
                       choices = c("Single Set of Gene Signatures - Dataframe" = "fileSingle",
                                   "Cancer Specific Signatures List - List (.rds)" = "fileList"
                       ),
                       selected = inType()
    )
    addTooltip(session,
             id = "fileType",
             "Upload New File \n to Change Gene Set Type"
             )


    delay(100, disable("fileType"))
  })


  ## Render UI For Input Type ##
  observeEvent(input$fileType,{
    if(is.null(input$fileType) | input$fileType == "fileSingle"){
      output$inSettings <- renderUI({
        radioButtons("sep", "Separator",
                     choices = c("Tab (.txt)" = "\t",
                                 "Comma (.csv)" = ",",
                                 "Semicolon (.csv)" = ";"
                     ),
                     selected = "\t",
                     inline = TRUE)
      })
    } else if (input$fileType == "fileList") {
      output$inSettings <- renderUI({
        selectInput(
          "listPrev",
          HTML("<h4> Select Cancer Type to Preview: </h4>"),
          choices = cancerList)
        })
    }
  })


  ## Keep Track of Active Input ##
  currentSig <- reactiveValues(reactInd = 0)

  observeEvent(input$loadCons, {
    currentSig$reactInd <- "sigConsensus"
  })

  observeEvent(input$geneFile, {
    currentSig$reactInd <- "sigInput"
  })

  ## Selected Correct Cancer Type Dependant on Dropdown Selection ##
  # Consensus Input
  consCanc <- reactive({
    input$consCanIn
  })
  # New Signature Input
  newSigCanc <- reactive({
    if(!is.null(input$listPrev)){
      return(input$listPrev)
    } else {
      'ACC'
    }
  })

  ## Load Gene Sets: Either from Consensus or Uploaded File ##

  # Single Signature Set Object
  newSigDf <- reactive({
    if(inType() == "fileSingle"){
      genesetFile <- input$geneFile
      if(!is.null(genesetFile)) {
        geneset <- read.table(genesetFile$datapath, header = TRUE, sep = input$sep)
        return(geneset)
      } else {
        return(NULL)
      }
    } else {
      return(NULL)
    }
  })

  # Multiple Signatures List
  newSigList <- reactive({
    if(inType() == "fileList"){
      genesetFile <- input$geneFile
      if(!is.null(genesetFile)){
        cancerGeneList <-readRDS(genesetFile$datapath)
        if(!is.list(cancerGeneList)){
          sendSweetAlert(
            session = session,
            title = "Gene List Error",
            text = "Uploaded .rds File Must Contain Gene Sets In A List Format (See Help Section For Guidance)",
            type = "warning"
          )
          return(NULL)
        }
        if(FALSE %in% unique(rapply(cancerGeneList, is.character))){
          sendSweetAlert(
            session = session,
            title = "Gene List Error",
            text = "Gene Lists Contains Objects of Type Other Than Character (See Help Section For Guidance)",
            type = "warning"
          )
          return(NULL)
        }
        if(FALSE %in% (names(cancerGeneList) == cancerList)){
          unused <- names(cancerGeneList[!names(cancerGeneList) == cancerList])
          missingSets <- as.character(cancerList)[!names(cancerGeneList) == cancerList]
          if(unused > 0 | !is.na(unused)){
            cancerGeneList[[unused]] <- NULL
            sendSweetAlert(
              session = session,
              html = TRUE,
              title = "Gene List Error",
              text = tags$span(
                tags$body(
                        "The Following Gene Sets Could Not Be Matched To One Of The Available TCGA Cancer Types So Will Removed:",
                        tags$br(),
                        tags$b(paste0('"', unused, '"')),
                        tags$br(),
                        "(See Help Section For Guidance)"
                        )
                ),
              type = "warning"
            )
          }
          if(missingSets > 0 & !is.na(missingSets)) {
            sendSweetAlert(
              session = session,
              title = "Gene List Error",
              text = sprintf("Gene List Names Need to Match Cancer Types, /n %s /n Could Not Be Found In Names of Uploaded Gene Sets Lists (See Help Section For Guidance)",
                             missingSets),
              type = "warning"
            )
          }
        }
        return(cancerGeneList)
      } else {
        return(NULL)
      }
    } else {
      return(NULL)
    }
  })


  # Gene Set For Previewing

  previewSet <- reactive({
    if(currentSig$reactInd == "sigConsensus"){
      consensusCancer <- readRDS('./data/ConsensusTME_GeneSets.rds')
      genelist <- ldply(consensusCancer[[consCanc()]], rbind)
      geneset <- t(genelist)
      colnames(geneset) <- geneset[1, ]
      geneset<- geneset[-1, ]
      geneset[is.na(geneset)] <- ''
      geneset <- as.data.frame(geneset)
      return(geneset)
    } else if(currentSig$reactInd == "sigInput") {
      if(input$fileType == "fileSingle"){
        geneset <- newSigDf()
        return(geneset)
      } else if(input$fileType == "fileList"){
        if(!is.null(newSigList())){
          cancerGeneList <- newSigList()
          genelist <- ldply(cancerGeneList[[newSigCanc()]], rbind)
          geneset <- t(genelist)
          colnames(geneset) <- geneset[1, ]
          geneset<- geneset[-1, ]
          geneset[is.na(geneset)] <- ''
          geneset <- as.data.frame(geneset)
          return(geneset)
        }
      }
    }
  })

  ## Render Previews

  # Render Dataframe
  output$tableHead <- renderTable({
    return(previewSet()[1:10, 1:5])
    })
  output$tableAll <- renderTable({
    return(previewSet())
  })



  ## Render Plotly Barchart for Gene Numbers ##
    output$sigsBar <- renderPlotly({
        geneSigs <- previewSet()
        sigChar <- as.data.frame(apply(geneSigs,2,
                                        FUN = function(x) length(x[nchar(x) > 0])))
        sigChar <- cbind(row.names(sigChar), sigChar)
        colnames(sigChar) <- c('Signature','Number_Of_Genes')
        row.names(sigChar) <- NULL

        nsigs <- nrow(sigChar)
        pal <- list()
        pal$n1 <- ifelse(nsigs < 9, 9 - nsigs, 9)
        pal$n2 <- ifelse(nsigs - 9 < 8, nsigs - 9, 8)
        pal$n3 <- ifelse(nsigs - 17 < 12, nsigs - 17, 12)

        pal <- lapply(pal, FUN = function (x) ifelse(x < 3, 3, x))
        sigCols = c(brewer.pal(name = 'Set1', n = pal$n1),
                     brewer.pal(name = 'Set2', n = pal$n2),
                     brewer.pal(name = 'Set3', n = pal$n3))
        sigCols <- sigCols[1:nsigs]

        sigPlot <- ggplot(data = sigChar,
                           aes(x = Signature,
                               y = Number_Of_Genes,
                               fill = Signature,
                               text = paste('Gene Signature:', gsub('_',' ',Signature),
                                                '<br>Number of Genes in Signature:', Number_Of_Genes))) +
                       geom_bar(stat = "identity") +
                       theme_minimal() +
                       theme(axis.text.x = element_text(angle = 45, vjust = 1.3, hjust = 1),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank()
                             ) +
                       scale_fill_manual(values = sigCols)
        intSigPlt <- ggplotly(sigPlot,
                                tooltip = "text") %>%
                       layout(yaxis = list(title = "Number of Genes in Signature"),
                              showlegend = FALSE) %>%
                       config(displayModeBar = F) %>%
        onRender("function(el, x) {
                       Plotly.d3.select('.cursor-crosshair').style('cursor', 'default')
            }"
        )
        intSigPlt
    })

  ## Render Overlap Matrix Of Genesets ##
  output$gsOverPltly <- renderPlotly({
      geneSigs <- previewSet()
      get_lower_tri<-function(mat){
        mat[upper.tri(mat, diag = F)] <- NA
        return(mat)
      }

      geneList <- convertRowsToList(t(geneSigs))
      overlapM <- as.data.frame(computeGeneSetsOverlap(geneList, unique(unlist(geneList))))
      overlapTri <- get_lower_tri(overlapM)
      overlapTri <- cbind(row.names(overlapTri), overlapTri)
      colnames(overlapTri)[1] <- 'sigs'

      overlapMelt <- reshape2::melt(overlapTri, na.rm = T)
      colnames(overlapMelt) <- c('sig1','sig2','overlap')
      overlapMelt$sig1 <- factor(overlapMelt$sig1, levels = rev(levels(overlapMelt$sig2)))

      corrPlt <- ggplot(data = overlapMelt,
                        aes(x = sig1,
                            y = sig2,
                            fill = overlap,
                            text = sprintf('Overlap of the %s & %s\n Gene Sets is: %s%%',
                                           gsub('_',' ', sig1),
                                           gsub('_',' ', sig2),
                                           format(round(overlap*100, 2), nsmall = 2)))) +
        geom_tile() +
        scale_fill_distiller(palette = "YlOrRd", direction = 1) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1.3, hjust = 1),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank())

      corrPltly <- ggplotly(corrPlt,
                             tooltip = "text") %>%
        config(displayModeBar = F)
  })



      output$previews <- renderUI({
        if(!is.null(previewSet())){
        tabsetPanel(
          tabPanel("Gene Numbers", plotlyOutput("sigsBar")),
          tabPanel("Overlap", plotlyOutput("gsOverPltly")),
          tabPanel("Table (Head)", tableOutput("tableHead")),
          tabPanel("Full Table", tableOutput("tableAll"))
        )
        } else {
          NULL
        }
      })

  #### Generate Estimates ####


      ## Create Reactive That Contains The Name of All signatures:
      ## For single geneset use colnames for list use unique union of all colnames

      output$isSelector <- renderUI({
        if(inType() == "fileSingle"){
          pickerInput(
            inputId = "immuneScorePicker",
            label = "Select All Immune Cell Signatures - These Will Be Used To Generate An Immune Score",
            choices = setNames(as.list(colnames(newSigDf())),gsub('_',' ',colnames(newSigDf()))),
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            ),
            multiple = TRUE
          )

        }
      })

 ## Ensure gene set name is valid and gene set is present
 validGsName <- eventReactive(input$runEstimation,{
    gsName <- input$geneSetName
    if(gsName == ''){
      sendSweetAlert(
        session = session,
        title = "Gene Set Name",
        text = "Gene set name is currently blank please enter a name",
        type = "warning"
      )
      return(FALSE)
      } else if(!is.character(gsName)){
        sendSweetAlert(
          session = session,
          title = "Gene Set Name",
          text = "Gene set name needs to be characters",
          type = "warning"
        )
        return(FALSE)
      } else if(grepl('[^A-z0-9_]', gsName)){
        sendSweetAlert(
          session = session,
          title = "Gene Set Name",
          text = "Gene set name cannot contain special characters please enter a new name",
          type = "warning"
        )
        return(FALSE)
      } else if(nchar(gsName) > 50){
        sendSweetAlert(
          session = session,
          title = "Gene Set Name",
          text = "Gene set name length limited to 50 characters please shorten",
          type = "warning"
        )
        return(FALSE)
      } else if(is.null(input$geneFile)){
        sendSweetAlert(
          session = session,
          title = "Gene Set",
          html = TRUE,
          text = tags$body(
            "A new gene set has not been uploaded",
            tags$br(),
            tags$br(),
            "(See Help Section For Guidance)"
          ),
          type = "warning"
        )
        return(FALSE)
        } else {
        return(TRUE)
        }
  })

 ## Generate Estimates

  ## Parallel Estimate Generation ##
 shinyjs::hide(id = "progressPanel")
 shinyjs::hide(id = "downloadPanel")

 # Status File

 status_file <- tempfile()

 # Update Status File
 set_status <- function(msg){
   write(msg, status_file, append = TRUE)
 }

 fire_ready <- function(){
   write("Ready", status_file, append = FALSE)
 }

 # Update Status File With Completed Cancers
 fire_complete <- function(canc){
   set_status(canc)
 }

 # Create Status File
 fire_ready()

 # Delete file at end of session
 onStop(function(){
   print(status_file)
   if(file.exists(status_file))
     unlink(status_file)
 })

 nclicks <- reactiveVal(0)

 # Create Reactive Output

 allEst <- reactiveValues()
 estimateProgress <- reactiveVal(FALSE)

 observeEvent(input$runEstimation,{

   # Don't do anything if analysis is already being run
   if(nclicks() != 0){
     showNotification("Already running analysis")
     return(NULL)
   }

   # Create reactive output
   shinyjs::disable(id = "runEstimation")

   estimateProgress(TRUE)

   # Assign all reactive values for futures

   gsBool <- validGsName()
   inTyp <- inType()
   newDf  <- newSigDf()
   isCell <- input$immuneScorePicker
   statFrame <- input$statMethod

   if(validGsName()){
     shinyjs::hide(id = "downloadPanel")
     estProgressBox <- shiny::Progress$new(session = session)
     estProgressBox$set(message = "Generating Estimates -", detail = "This May Take Several Minutes Please Explore Current Benchmarking While You Wait")

     futs <- promise_all(.list = lapply(cancerTypes, function(canc){
       future({
         lazyLoad("./data/All_RNA")
         if(inTyp == "fileSingle"){
           geneset <- t(newDf)
           genelist <- setNames(split(geneset, seq(nrow(geneset))), rownames(geneset))
           if(!is.null(isCell)){
             genelist[['Immune_Score']] <- unique(unlist(genelist[isCell]))
           }
           genelist <- removeBlanks(genelist)
           tmpM <- get(canc)
           if (statFrame == "statSsgsea") {
             estimates <- gsva(tmpM,
                               genelist,
                               method='ssgsea',
                               min.sz=0,
                               max.sz=Inf,
                               ssgsea.norm=T)
             cancEst <- as.data.frame(estimates)
           } else if (statFrame == "statSingscore") {
             rankRna <- rankGenes(tmpM)
             genelist <- genelist[sapply(genelist, function(x){
               length(x) > 0
             })]
             cancerGeneSet <- lapply(names(genelist), function(cellType){
               GeneSet(genelist[[cellType]], setName = cellType)
             })
             cancerGeneCol <- GeneSetCollection(cancerGeneSet)
             scores <- multiScore(rankData = rankRna, upSetColc = cancerGeneCol)
             cancEst <- as.data.frame(scores$Scores)
           }
           fire_complete(canc)
           rm(tmpM)
           gc()
           return(cancEst)
         }
       }) %...>% {
           estProgressBox$inc(1/(length(cancerTypes)), message = "Generating Estimates:  ", detail = 'This may take several minutes - See Estimation Panel For More Detail')
           allEst[[canc]] <<- .
         }
     }
     )
     )
   }

   futs <- finally(futs,
                     function(){
                       fire_ready()
                       nclicks(0)
                       sendSweetAlert(
                                 session = session,
                                 title = "Estimates Generated",
                                 html = TRUE,
                                 text = tags$body(
                                   "Estimates generated for TCGA and benchmarking datasets",
                                   tags$br(),
                                   tags$br(),
                                   tags$b("To Avoid Repeating Long Generation Process For This Geneset Please Download Results")
                                 ),
                                 type = "success",
                                 closeOnClickOutside = FALSE
                               )
                       shinyjs::hide(id = "progressPanel")
                       shinyjs::enable(id = "runEstimation")
                       estimateProgress(FALSE)
                       estProgressBox$close()
                       shinyjs::show(id = "downloadPanel")

                     })

   # Return Something To Keep Shiny Responsive
   NULL

 })

 # Observe tmp file for progress bar only during running

 autoInvalidate <- reactiveTimer(2000)

 datasetsProgress <- reactiveValues()

 # Progress update via reactive timer based observer reading tmp file
 observe({
   req(estimateProgress())
   shinyjs::show(id = "progressPanel")

   if(estimateProgress()){
     autoInvalidate()

     tmpContents <- scan(status_file, what = "character",sep="\n")

     datasetsProgress$toDo <- cancerTypes[!cancerTypes %in% tmpContents]

     datasetsProgress$finished <- cancerTypes[cancerTypes %in% tmpContents]

     percentComplete <- length(datasetsProgress$finished)/length(cancerTypes)

     updateProgressBar(session = session, id = "estProgressBar", value = length(datasetsProgress$finished), total = length(cancerTypes))

   }
 })

 # Render table of datasets to do and completed

 output$progressToDo <- renderUI({

   HTML(paste0("<h4> Generation in Progress For Datasets: <h4/>", paste(datasetsProgress$toDo, sep = "<br/>", collapse = " ")))
 })

 output$progressFinished<- renderUI({
   HTML(paste0("<h4> Finished: <h4/>", paste(datasetsProgress$finished, sep = "<br/>", collapse = " ")))
 })


  ## Loaded Estimates
  loadEst <- eventReactive(input$estFile, {
    estLoadFile <- input$estFile
    loadedEst <- readRDS(estLoadFile$datapath)
    return(loadedEst)
  })

  ## Keep Track of Last Input ##
  currentEst <- reactiveValues(reactInd = 0)

  observeEvent(input$runEstimation, {
    if(validGsName()) {
      currentEst$reactInd <- "estGenerated"
    }
  })

  observeEvent(input$estFile, {
    currentEst$reactInd <- "estLoaded"
  })

  ## Assign Most Recent Input To Be Activate Estimate Set

  activeEst <- reactive({
    if(currentEst$reactInd == "estGenerated") {
      return(isolate(reactiveValuesToList(allEst)))
    } else if (currentEst$reactInd == "estLoaded")
      return(loadEst())
  })

  ## Observer to enable/disable generate estimates panel if geneset hasn't been uploaded
  observe(
    if(is.null(input$geneFile)){
      disable("estSidebar")
      addTooltip(session,
                 id = "estSidebar",
                 "Upload Gene Set To Enable Panel"
      )
    } else {
      enable("estSidebar")
      addTooltip(session,
                 id = "estSidebar",
                 ""
      )
    }
  )

  ## Txt Download Handler

  output$estTxtDownload <- downloadHandler(
    filename = function() {
      paste0(input$geneSetName, "_Estimates_", format(Sys.Date(), format = "%d%b%y"), ".zip")
    },
    content = function(file) {
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      files <- NULL;

      estimates <- activeEst()
      for(i in 1:length(estimates)){
        tmpDf <- estimates[[i]]
        fileName <- paste0(input$geneSetName,
                           "_",
                           names(estimates)[i],
                           "_",
                           format(Sys.Date(), format = "%d%b%y"),
                           ".txt")
        write.table(tmpDf,
                    fileName,
                    sep = "\t",
                    col.names = NA)
        files <- c(fileName, files)
      }
      zip(file, files)
    }
  )

  ## rds Download Handler

  output$estRdsDownload <- downloadHandler(
    filename = function() {
      paste0(input$geneSetName, "_Estimates_", format(Sys.Date(), format = "%d%b%y"), ".rds")
    },
    content = function(file) {
      write_rds(activeEst(), file, compress = "gz", compression = 9L)
    }
  )

  #### Explore Estimates ####

  ## Panel To Explore Generated Estimates

  output$estPreviews <- renderUI({
    if (estimateProgress()) {
      sidebarPanel(
        HTML("<center> <h3> Generating Estimates For Benchmarking Datasets </h3> </center>"),
        offset = 1,
        width = 11
      )
    } else if (!is.null(activeEst())){
      fluidPage(
        tabsetPanel(
          tabPanel("Heatmap",
                   fluidRow(
                     dropdown(
                       tags$h3("Heatmap Options"),
                       tags$hr(),
                       pickerInput(
                         inputId = "nesPrevSelection",
                         label = "Select Dataset To Preview:",
                         choices = list(
                           TCGA = cancerTypes,
                           Independent = c("MCP-Counter", "TIMER", "CIBERSORT")
                         )
                       ),
                       checkboxGroupButtons(
                         inputId = "nesCluster",
                         label = "Cluster By:",
                         choiceNames = c("Row", "Column"),
                         choiceValues = c("row", "column"),
                         individual = TRUE,
                         checkIcon = list(
                           yes = icon("ok", lib = "glyphicon"),
                           no = icon("remove", lib = "glyphicon")
                         ),
                         selected = c("row", "column")
                       ),
                       style = "pill",
                       status = "success",
                       icon = icon("gear"),
                       size = "lg"
                       ),
                     uiOutput("heatNextPrev")
                     ),
                   withSpinner(plotlyOutput("nesPlotly"),
                                          type = 6)
                   ),
          tabPanel("Correlations",
                   fluidRow(
                     dropdown(
                       tags$h3("Correlation Plot Options"),
                       tags$hr(),
                       pickerInput(
                         inputId = "nesCorrSelection",
                         label = "Select Dataset To Preview:",
                         choices = list(
                           TCGA = cancerTypes,
                           Independent = c("MCP-Counter", "TIMER", "CIBERSORT")
                         )
                       ),
                       radioGroupButtons(
                         inputId = "corrTestSelection",
                         label = "Select Correlation Test",
                         choiceNames = c("Pearson", "Spearman", "Kendall"),
                         choiceValues = c("pearson", "spearman", "kendall"),
                         selected = "kendall",
                         individual = TRUE
                       ),
                       style = "pill",
                       status = "success",
                       icon = icon("gear"),
                       size = "lg"
                     ),
                     uiOutput("corrNextPrev")
                   ),
                   withSpinner(plotlyOutput("estCorrPltly"),
                               type = 6)
          ),
          tabPanel("Normalised Enrichment Scores (NES) (Head)",
                   do.call(tabsetPanel,
                           lapply(cancerTypes, function(cancer) {
                             tabPanel(title = cancer, withSpinner(tableOutput(paste0("estTableHead", cancer)),
                                                                  type = 6))
                             })
                           )
                   ),
          tabPanel("NES - Full Table",
                   do.call(tabsetPanel,
                           lapply(cancerTypes, function(cancer) {
                             tabPanel(title = cancer, withSpinner(tableOutput(paste0("estTableAll", cancer)),
                                                                  type = 6))
                           })
                   )
          )
        )
      )
    } else {
      sidebarPanel(
        HTML("<center> <h3> <font color=\"#8f9399\">Generate or Load New Estimates To Preview Estimates </font> </h3> </center>"),
        width = 12
      )
    }
  })


  ## Render Table Output
  lapply(cancerTypes, function(cancer) {
  output[[paste0("estTableHead", cancer)]] <- renderTable({
    estimates <- activeEst()
    if(is.null(input$estPrevDataset)){
      return(estimates[[cancer]][1:nrow(estimates[[1]]),1:5])
    } else {
      return(estimates[[input$estPrevDataset]][1:nrow(estimates[[1]]),1:5])
    }
  },rownames = TRUE)})

  lapply(cancerTypes, function(cancer) {
    output[[paste0("estTableAll", cancer)]] <- renderTable({
      estimates <- activeEst()
      if(is.null(input$estPrevDataset)){
        return(estimates[[cancer]])
      } else {
        return(estimates[[input$estPrevDataset]])
      }
    },rownames = TRUE)})

  ## Render Interactive Heatmap With Options

  output$nesPlotly <- renderPlotly({
    estimates <- isolate(activeEst())
    mat <- estimates[[input$nesPrevSelection]]
    if (length(input$nesCluster) == 0){
       dend_clust <- "none"
     } else if (length(input$nesCluster) == 1) {
       dend_clust <- input$nesCluster
     } else {
       dend_clust <- "both"
     }
    heatmaply(mat,
              colors = PRGn,
              dendrogram = dend_clust,
              branches_lwd = 0.3,
              grid_gap = 0.5,
              xlab = "Tumour Samples \n <sup> (hover for details) </sup>",
              main = input$nesPrevSelection,
              showticklabels = c(F,T),
              label_names = c("Immune Cell", "Sample ID", "NES"),
              key.title='Normalised Enrichment \n Score (NES)') %>%
      layout(xaxis=list(fixedrange=TRUE)) %>%
      layout(xaxis2=list(fixedrange=TRUE)) %>%
      layout(yaxis=list(fixedrange=TRUE)) %>%
      layout(yaxis2=list(fixedrange=TRUE)) %>%
      layout(height=700) %>%
      config(displayModeBar = F) %>%
      layout(paper_bgcolor='transparent')
    })
  
  ## Render Correlation Matrix
  
  output$estCorrPltly <- renderPlotly({
    estimates <- isolate(activeEst())
    mat <- t(estimates[[input$nesCorrSelection]])
    
    get_lower_tri<-function(mat){
      mat[upper.tri(mat, diag = F)] <- NA
      return(mat)
    }
    
    corrMat <- cor(mat, method = input$corrTestSelection)
    corrMat <- round(corrMat, 2)
    corrTri <- get_lower_tri(corrMat)
    corrMelt <- reshape2::melt(corrTri, na.rm = T)
    colnames(corrMelt) <- c("sig1", "sig2", "correlation")
    corrMelt$sig1 <- factor(corrMelt$sig1, levels = rev(levels(corrMelt$sig2)))
    
    limit <- max(abs(corrMelt$correlation)) * c(-1, 1)
    
    corrPlt <- ggplot(data = corrMelt,
                      aes(x = sig1,
                          y = sig2,
                          fill = correlation,
                          text = sprintf("Correlation of %s & %s\n NESs is: %s",
                                         gsub("_"," ", sig1),
                                         gsub("_"," ", sig2),
                                         correlation))) +
      geom_tile() +
      scale_fill_distiller(type = "div", palette = "RdBu", limit = limit) +
      theme_minimal() +
      ggtitle(input$nesCorrSelection) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1.3, hjust = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(hjust = 0.5))
    
    ggplotly(corrPlt,
             tooltip = "text") %>%
      config(displayModeBar = F) %>%
      layout(paper_bgcolor='transparent')
    
  })
  
  ## Heatmaps Next & Previous Buttons
  
  output$heatNextPrev <- renderUI({
    plotList <- cancerTypes
    nbTab <- length(plotList)
    if (match(input$nesPrevSelection, cancerTypes) == nbTab)
      column(1, offset = 1, heatPrevButton)
    else if (match(input$nesPrevSelection, cancerTypes) == 1)
      column(1, offset = 10, heatNextButton)
    else
      div(column(1, offset = 1, heatPrevButton), column(1, offset=8, heatNextButton))
  })
  
  output$corrNextPrev <- renderUI({
    plotList <- cancerTypes
    nbTab <- length(plotList)
    if (match(input$nesCorrSelection, cancerTypes) == nbTab)
      column(1, offset = 1, corrPrevButton)
    else if (match(input$nesCorrSelection, cancerTypes) == 1)
      column(1, offset = 10, corrNextButton)
    else
      div(column(1, offset = 1, corrPrevButton), column(1, offset=8, corrNextButton))
  })
  
  observeEvent(input$heatPrevPlot,
               {
                 plotList <- cancerTypes
                 currentPlot <- match(input$nesPrevSelection, cancerTypes)
                 updatePickerInput(session,"nesPrevSelection", selected = plotList[currentPlot - 1])
               })
  
  observeEvent(input$heatNextPlot,
               {
                 plotList <- cancerTypes
                 currentPlot <- match(input$nesPrevSelection, cancerTypes)
                 updatePickerInput(session,"nesPrevSelection", selected = plotList[currentPlot + 1])
               })
  
  observeEvent(input$corrPrevPlot,
               {
                 plotList <- cancerTypes
                 currentPlot <- match(input$nesCorrSelection, cancerTypes)
                 updatePickerInput(session,"nesCorrSelection", selected = plotList[currentPlot - 1])
               })
  
  observeEvent(input$corrNextPlot,
               {
                 plotList <- cancerTypes
                 currentPlot <- match(input$nesCorrSelection, cancerTypes)
                 updatePickerInput(session,"nesCorrSelection", selected = plotList[currentPlot + 1])
               })

  #### Tumour Purity Benchmark ####
  
  load("./data/Tumour_Purity_Corrs.RData")
  load("./data/TCGA_Annotation.RData")
  leukSideCols <- side_cols[!row.names(side_cols) %in% c("THYM", "DLBC"), ]
  ## Plot Heatmap 
  
  output$purCorrPlotly <- renderPlotly({
    purWideTau <- wideTau
    purWideTau <- as.data.frame(sapply(purWideTau, as.numeric))
    row.names(purWideTau) <- row.names(wideQval)
    row.names(purWideTau) <- gsub("ConsensusTME", "Consensus<sup>TME</sup>", row.names(purWideTau))
    heatmaply(as.matrix(purWideTau), 
              colors = PRGn(100)[1:50],
              method = "plotly",
              dendrogram = "none",
              grid_gap = 2,
              col_side_colors = side_cols,
              col_side_palette = side_cols_comb,
              subplot_heights = c(0.2, 0.8),
              key.title = "Immune Score vs\n Tumour Purity\n (Kendalls &#964;)",
              label_names = c("Method", "Cancer", "Correlation"),
              ylab = "Estimation Method",
              xlab = "TCGA Cancer Type",
              heatmap_layers = theme(axis.text.x = element_text(angle = 90),
                                     axis.line = element_blank())) %>%
      layout(showlegend = FALSE) %>%
      config(displayModeBar = F)
  })
  
  ## Plot Purity Correlations Boxplot
  output$purCorrBoxPlotly <- renderPlotly({
    purLongCorrs <- longCorrs
    purLongCorrs <- purLongCorrs %>% dplyr::distinct()
    purLongCorrs$Method <- gsub("_Immune_Score", "", purLongCorrs$Method)
    purLongCorrs$Method <- gsub("MCP", "MCP-Counter", purLongCorrs$Method)
    purLongCorrs$Method <- gsub("ConsensusTME", "Consensus<sup>TME</sup>", purLongCorrs$Method)
    
    methodCols <- c("#F6EB16", "#6FCCDD", "#3953A4" , "#087022",
                     "#231F20", "#A65628", "#B9529F", "#ED2224", "#F7931D")
    
    names(methodCols) <- unique(purLongCorrs$Method)
    
    qualColPals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    colVector = unlist(mapply(brewer.pal, qualColPals$maxcolors, rownames(qualColPals)))
    cancCols <- rep(colVector, 9)
    
    names(cancCols) <- unique(purLongCorrs$Cancer)
    
    plot <- ggplot(data = purLongCorrs,
                   aes(x = reorder(Method, Tau, FUN = median),
                       y = Tau, labels = Method)) +
      geom_boxplot(fill = methodCols[levels(reorder(purLongCorrs$Method,
                                                    purLongCorrs$Tau,
                                                    FUN = median))],
                   outlier.shape = NA) +
      geom_beeswarm(aes(fill = Cancer),
                    shape = 21,
                    alpha = 0.65,
                    cex = 1.75,
                    priority = "ascending"
      ) +
      theme_minimal() +
      xlab("Estimation Method") +
      ylab("Correlation Between RNA & IHC Estimates <br> <sup>(Kendall's &#964; )</sup>")
      theme(axis.text.x = element_text(angle = 45, vjust = 1.3, hjust = 1),
            panel.grid.minor = element_blank()
      ) 
    
    p <- ggplotly(plot, tooltip = c("Tau", "Method", "Cancer"))
    
    p$x$data[1:9] <- lapply(p$x$data[1:9], FUN = function(x){
      x$marker = list(opacity = 0)
      return(x)
    })
    
    p %>%
      config(displayModeBar = F)
  })
  
  
  #### Leukocyte Methylation Benchmark ####
  
  ## Load Data
  
  load("./data/Methylation_Fit_Metrics.RData")
  
  colnames(longFitMetric) <- gsub("_Z", " Z-Score", colnames(longFitMetric))
  
  colnames(longFitMetric) <- gsub("_R_Squared", " R<sup>2</sup>", colnames(longFitMetric))
  
  ## Fit Metric Heatmaps
  
  ## Adjusted R-Squared
  
  output$leukR2FitHeatmap <- renderPlotly({
    
    leukWideR2 <- wideRSquared
    
    leukWideR2 <- as.data.frame(sapply(wideRSquared, as.numeric))
    
    row.names(leukWideR2) <- row.names(wideRSquared)
    
    leukSideCols <- side_cols[!row.names(side_cols) %in% c("THYM", "DLBC"), ]
    
    row.names(leukWideR2) <- gsub("consensusTME_v2", "Consensus<sup>TME</sup>", row.names(leukWideR2))
    row.names(leukWideR2) <- gsub("MCPcounter", "MCP-Counter", row.names(leukWideR2))
    
    method_cols <- c("#F6EB16", "#6FCCDD", "#3953A4" , "#087022",
                     "#231F20", "#A65628", "#B9529F", "#ED2224", "#F7931D")
    
    names(method_cols) <- unique(longFitMetric$Method)
    
    heatmaply(as.matrix(leukWideR2), 
              colors = BrBG(100)[50:100],
              method = "plotly",
              dendrogram = "none",
              grid_gap = 2,
              col_side_colors = leukSideCols,
              col_side_palette = side_cols_comb,
              subplot_heights = c(0.2, 0.8),
              key.title = "Leukocyte Fraction ~ \nRNA Cell-Type Estimates \n(Adjusted R<sup>2</sup>)",
              label_names = c("Method", "Cancer", "R<sup>2</sup>"),
              ylab = "Estimation Method",
              xlab = "TCGA Cancer Type",
              heatmap_layers = theme(axis.text.x = element_text(angle = 90),
                                     axis.line = element_blank())) %>%
      layout(showlegend = FALSE) %>%
      config(displayModeBar = F)

  })
  
  ## AIC
  
  output$leukAicFitHeatmap <- renderPlotly({
    leukWideAic <- wideAIC
    
    leukWideAic <- as.data.frame(sapply(leukWideAic, as.numeric))
    
    row.names(leukWideAic) <- row.names(wideAIC)
    
    leukSideCols <- side_cols[!row.names(side_cols) %in% c("THYM", "DLBC"), ]
    
    row.names(leukWideAic) <- gsub("consensusTME_v2", "Consensus<sup>TME</sup>", row.names(leukWideAic))
    row.names(leukWideAic) <- gsub("MCPcounter", "MCP-Counter", row.names(leukWideAic))
    
    heatmaply(as.matrix(leukWideAic), 
              colors = RdBu(100),
              method = "plotly",
              dendrogram = "none",
              grid_gap = 2,
              col_side_colors = leukSideCols,
              col_side_palette = side_cols_comb,
              subplot_heights = c(0.2, 0.8),
              key.title = "Leukocyte Fraction ~ \nRNA Cell-Type Estimates \n(<b>A</b>kaike <b>I</b>nformation <b>C</b>riterion)",
              label_names = c("Method", "Cancer", "AIC"),
              ylab = "Estimation Method",
              xlab = "TCGA Cancer Type",
              heatmap_layers = theme(axis.text.x = element_text(angle = 90),
                                     axis.line = element_blank())) %>%
      layout(showlegend = FALSE) %>%
      config(displayModeBar = F)
    
  })
  
  ## BIC
  
  output$leukBicFitHeatmap <- renderPlotly({
    leukWideBic <- wideBIC
    
    leukWideBic <- as.data.frame(sapply(leukWideBic, as.numeric))
    
    row.names(leukWideBic) <- row.names(wideBIC)
    
    leukSideCols <- side_cols[!row.names(side_cols) %in% c("THYM", "DLBC"), ]
    
    row.names(leukWideBic) <- gsub("consensusTME_v2", "Consensus<sup>TME</sup>", row.names(leukWideBic))
    
    row.names(leukWideBic) <- gsub("MCPcounter", "MCP-Counter", row.names(leukWideBic))
    
    heatmaply(as.matrix(leukWideBic), 
              colors = RdBu(100),
              method = "plotly",
              dendrogram = "none",
              grid_gap = 2,
              col_side_colors = leukSideCols,
              col_side_palette = side_cols_comb,
              subplot_heights = c(0.2, 0.8),
              key.title = "Leukocyte Fraction ~ \nRNA Cell-Type Estimates \n(<b>B</b>ayesian <b>I</b>nformation <b>C</b>riterion)",
              label_names = c("Method", "Cancer", "BIC"),
              ylab = "Estimation Method",
              xlab = "TCGA Cancer Type",
              heatmap_layers = theme(axis.text.x = element_text(angle = 90),
                                     axis.line = element_blank())) %>%
      layout(showlegend = FALSE) %>%
      config(displayModeBar = F)
    
  })
  
  ## Fit Metrics Boxplots
  
  ## Adjusted R-Squared
  output$leukR2FitBoxplot <- renderPlotly({
    
    leukLong <- longFitMetric %>% dplyr::distinct()
    
    leukLong$Method <- gsub("MCPcounter", "MCP-Counter", leukLong$Method)
    
    leukLong$Method <- gsub("consensusTME_v2", "Consensus<sup>TME</sup>", leukLong$Method)
    
    leukMethodCols <- c("#ED2224", "#6FCCDD", "#F6EB16" , "#231F20",
                        "#087022", "#F7931D", "#B9529F", "#3953A4")
    
    names(leukMethodCols) <- unique(leukLong$Method)
    
    qualColPals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    colVector = unlist(mapply(brewer.pal, qualColPals$maxcolors, rownames(qualColPals)))
    cancCols <- rep(colVector, 8)
    
    leukR2plot <- ggplot(data = leukLong,
                         aes(x = reorder(Method, `Adjusted R<sup>2</sup>`, FUN = function(x) -(median(x, na.rm = TRUE))),
                             y = `Adjusted R<sup>2</sup>`, labels = Method)) +
      geom_boxplot(fill = leukMethodCols[levels(reorder(leukLong$Method,
                                                        leukLong$`Adjusted R<sup>2</sup>`,
                                                        FUN = function(x) -(median(x, na.rm = TRUE))))],
                   outlier.shape = NA) +
      geom_beeswarm(aes(fill = Cancer),
                    shape = 21,
                    alpha = 0.65,
                    cex = 1.75,
                    priority = "ascending"
      ) +
      theme_minimal() +
      xlab("Estimation Method") +
      theme(axis.text.x = element_text(angle = 45, vjust = 1.3, hjust = 1),
            panel.grid.minor = element_blank()
      ) 
    
    pR2 <- ggplotly(leukR2plot, tooltip = c("Adjusted R<sup>2</sup>", "Method", "Cancer"))
    
    pR2$x$data[1:8] <- lapply(pR2$x$data[1:8], FUN = function(x){
      x$marker = list(opacity = 0)
      return(x)
    })
    
    pR2 %>%
      config(displayModeBar = F)
  })
  
  ## AIC
  output$leukAicFitBoxplot <- renderPlotly({
    
    leukLong <- longFitMetric %>% dplyr::distinct()
    
    leukLong$Method <- gsub("MCPcounter", "MCP-Counter", leukLong$Method)
    
    leukLong$Method <- gsub("consensusTME_v2", "Consensus<sup>TME</sup>", leukLong$Method)
    
    leukMethodCols <- c("#ED2224", "#6FCCDD", "#F6EB16" , "#231F20",
                        "#087022", "#F7931D", "#B9529F", "#3953A4")
    
    names(leukMethodCols) <- unique(leukLong$Method)
    
    qualColPals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    colVector = unlist(mapply(brewer.pal, qualColPals$maxcolors, rownames(qualColPals)))
    cancCols <- rep(colVector, 8)
    
    leukAicplot <- ggplot(data = leukLong,
                         aes(x = reorder(Method, `AIC Z Score`, FUN = function(x) (median(x, na.rm = TRUE))),
                             y = `AIC Z Score`, labels = Method)) +
      geom_boxplot(fill = leukMethodCols[levels(reorder(leukLong$Method,
                                         leukLong$`AIC Z Score`,
                                         FUN = function(x) (median(x, na.rm = TRUE))))],
                   outlier.shape = NA) +
      geom_beeswarm(aes(fill = Cancer),
                    shape = 21,
                    alpha = 0.65,
                    cex = 1.75,
                    priority = "ascending"
      ) +
      theme_minimal() +
      xlab("Estimation Method") +
      theme(axis.text.x = element_text(angle = 45, vjust = 1.3, hjust = 1),
            panel.grid.minor = element_blank()
      ) 
    
    pAic <- ggplotly(leukAicplot, tooltip = c("AIC Z Score", "Method", "Cancer"))
    
    pAic$x$data[1:8] <- lapply(pAic$x$data[1:8], FUN = function(x){
      x$marker = list(opacity = 0)
      return(x)
    })
    
    pAic %>%
      config(displayModeBar = F)
  })
  
  ## BIC
  
  output$leukBicFitBoxplot <- renderPlotly({
    
    leukLong <- longFitMetric %>% dplyr::distinct()
    
    leukLong$Method <- gsub("MCPcounter", "MCP-Counter", leukLong$Method)
    
    leukLong$Method <- gsub("consensusTME_v2", "Consensus<sup>TME</sup>", leukLong$Method)
    
    
    leukMethodCols <- c("#ED2224", "#6FCCDD", "#F6EB16" , "#231F20",
                        "#087022", "#F7931D", "#B9529F", "#3953A4")
    
    names(leukMethodCols) <- unique(leukLong$Method)
    
    qualColPals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    colVector = unlist(mapply(brewer.pal, qualColPals$maxcolors, rownames(qualColPals)))
    cancCols <- rep(colVector, 8)
    
    leukBicplot <- ggplot(data = leukLong,
                          aes(x = reorder(Method, `BIC Z Score`, FUN = function(x) (median(x, na.rm = TRUE))),
                              y = `BIC Z Score`, labels = Method)) +
      geom_boxplot(fill = leukMethodCols[levels(reorder(leukLong$Method,
                                                        leukLong$`BIC Z Score`,
                                                        FUN = function(x) (median(x, na.rm = TRUE))))],
                   outlier.shape = NA) +
      geom_beeswarm(aes(fill = Cancer),
                    shape = 21,
                    alpha = 0.65,
                    cex = 1.75,
                    priority = "ascending"
      ) +
      theme_minimal() +
      xlab("Estimation Method") +
      theme(axis.text.x = element_text(angle = 45, vjust = 1.3, hjust = 1),
            panel.grid.minor = element_blank()
      ) 
    
    pBic <- ggplotly(leukBicplot, tooltip = c("BIC Z Score", "Method", "Cancer"))
    
    pBic$x$data[1:8] <- lapply(pBic$x$data[1:8], FUN = function(x){
      x$marker = list(opacity = 0)
      return(x)
    })
    
    pBic %>%
      config(displayModeBar = F)
  })

  #### Image Analysis Benchmark ####
  load("./data/Image_Analysis_Fit_Metrics.RData")
  
  colnames(imageLongFitMetric) <- gsub("_Z", " Z Score",
                                       colnames(imageLongFitMetric))
  
  colnames(imageLongFitMetric) <- gsub("_R_Squared", " R<sup>2</sup>",
                                       colnames(imageLongFitMetric))
  
  ## Fit Metrics Heatmap
  
  ## Adjusted R-Squared
  
  output$imageR2FitHeatmap <- renderPlotly({
    
    imageWideR2 <- as.data.frame(sapply(imageWideRSquared, as.numeric))
    
    row.names(imageWideR2) <- row.names(imageWideRSquared)
    
    imageSideCols <- side_cols[row.names(side_cols) %in% colnames(imageWideR2), ]
    
    row.names(imageWideR2) <- gsub("consensusTME_v2", "Consensus<sup>TME</sup>", row.names(imageWideR2))
    row.names(imageWideR2) <- gsub("MCPcounter", "MCP-Counter", row.names(imageWideR2))
    
    method_cols <- c("#F6EB16", "#6FCCDD", "#3953A4" , "#087022",
                     "#231F20", "#A65628", "#B9529F", "#ED2224", "#F7931D")
    
    names(method_cols) <- unique(imageLongFitMetric$Method)
    
    heatmaply(as.matrix(imageWideR2), 
              colors = BrBG(100)[50:100],
              method = "plotly",
              dendrogram = "none",
              grid_gap = 2,
              col_side_colors = imageSideCols,
              col_side_palette = side_cols_comb,
              subplot_heights = c(0.2, 0.8),
              key.title = "Lymphocyte Score ~ \nRNA Cell-Type Estimates \n(Adjusted R<sup>2</sup>)",
              label_names = c("Method", "Cancer", "R<sup>2</sup>"),
              ylab = "Estimation Method",
              xlab = "TCGA Cancer Type",
              heatmap_layers = theme(axis.text.x = element_text(angle = 90),
                                     axis.line = element_blank())) %>%
      layout(showlegend = FALSE) %>%
      config(displayModeBar = F)
    
  })
  
  ## AIC
  
  output$imageAicFitHeatmap <- renderPlotly({
    
    imageAic <- as.data.frame(sapply(imageWideAIC, as.numeric))
    
    row.names(imageAic) <- row.names(imageWideAIC)
    
    imageSideCols <- side_cols[row.names(side_cols) %in% colnames(imageAic), ]
    
    row.names(imageAic) <- gsub("consensusTME_v2", "Consensus<sup>TME</sup>", row.names(imageAic))
    row.names(imageAic) <- gsub("MCPcounter", "MCP-Counter", row.names(imageAic))
    
    method_cols <- c("#F6EB16", "#6FCCDD", "#3953A4" , "#087022",
                     "#231F20", "#A65628", "#B9529F", "#ED2224", "#F7931D")
    
    names(method_cols) <- unique(imageLongFitMetric$Method)
    
    heatmaply(as.matrix(imageAic), 
              colors = RdBu(100),
              method = "plotly",
              dendrogram = "none",
              grid_gap = 2,
              col_side_colors = imageSideCols,
              col_side_palette = side_cols_comb,
              subplot_heights = c(0.2, 0.8),
              key.title = "Leukocyte Fraction ~ \nRNA Cell-Type Estimates \n(<b>A</b>kaike <b>I</b>nformation <b>C</b>riterion)",
              label_names = c("Method", "Cancer", "AIC"),
              ylab = "Estimation Method",
              xlab = "TCGA Cancer Type",
              heatmap_layers = theme(axis.text.x = element_text(angle = 90),
                                     axis.line = element_blank())) %>%
      layout(showlegend = FALSE) %>%
      config(displayModeBar = F)
    
  })
  
  output$imageBicFitHeatmap <- renderPlotly({
    
    imageBic <- as.data.frame(sapply(imageWideBIC, as.numeric))
    
    row.names(imageBic) <- row.names(imageWideBIC)
    
    imageSideCols <- side_cols[row.names(side_cols) %in% colnames(imageBic), ]
    
    row.names(imageBic) <- gsub("consensusTME_v2", "Consensus<sup>TME</sup>", row.names(imageBic))
    row.names(imageBic) <- gsub("MCPcounter", "MCP-Counter", row.names(imageBic))
    
    method_cols <- c("#F6EB16", "#6FCCDD", "#3953A4" , "#087022",
                     "#231F20", "#A65628", "#B9529F", "#ED2224", "#F7931D")
    
    names(method_cols) <- unique(imageLongFitMetric$Method)
    
    heatmaply(as.matrix(imageBic), 
              colors = RdBu(100),
              method = "plotly",
              dendrogram = "none",
              grid_gap = 2,
              col_side_colors = imageSideCols,
              col_side_palette = side_cols_comb,
              subplot_heights = c(0.2, 0.8),
              key.title = "Leukocyte Fraction ~ \nRNA Cell-Type Estimates \n(<b>B</b>ayesian <b>I</b>nformation <b>C</b>riterion)",
              label_names = c("Method", "Cancer", "BIC"),
              ylab = "Estimation Method",
              xlab = "TCGA Cancer Type",
              heatmap_layers = theme(axis.text.x = element_text(angle = 90),
                                     axis.line = element_blank())) %>%
      layout(showlegend = FALSE) %>%
      config(displayModeBar = F)
    
  })
  
  ## Summary Boxplots
  
  ## Adjusted R-Squared
  
  output$imageR2FitBoxplot <- renderPlotly({
    
    imLong <- imageLongFitMetric %>% dplyr::distinct()
    
    imLong$Method <- gsub("MCPcounter", "MCP-Counter", imLong$Method)
    
    imLong$Method <- gsub("consensusTME_v2", "Consensus<sup>TME</sup>", imLong$Method)
    
    imMethodCols <- c("#ED2224", "#6FCCDD", "#F6EB16" , "#231F20",
                      "#087022", "#F7931D", "#B9529F", "#3953A4")
    
    names(imMethodCols) <- unique(imLong$Method)
    
    imMethodColScale <- scale_fill_manual(name = "Method ", values = imMethodCols)
    
    cancCols <- rep(c("#313131", "#776f7a", "#0436bf", "#3ee5ea", "#04601d",
                      "#0ce747", "#c0cc00", "#e5f656", "#f9a53d", "#760606",
                      "#f94e4e","#85037f","#bd3eea"), 8)
    
    imR2plot <- ggplot(data = imLong,
                       aes(x = reorder(Method, `Adjusted R<sup>2</sup>`, FUN = function(x) -(median(x, na.rm = TRUE))),
                           y = `Adjusted R<sup>2</sup>`, labels = Method)) +
      geom_boxplot(fill = imMethodCols[levels(reorder(imLong$Method, imLong$`Adjusted R<sup>2</sup>`, FUN = function(x) -(median(x, na.rm = TRUE))))],
                   outlier.shape = NA) +
      geom_beeswarm(aes(fill = Cancer),
                    shape = 21,
                    alpha = 0.65,
                    cex = 1.75,
                    priority = "ascending"
      ) +
      
      theme_minimal() +
      xlab("Estimation Method") +
      theme(axis.text.x = element_text(angle = 45, vjust = 1.3, hjust = 1),
            panel.grid.minor = element_blank()
      ) 
    
    pR2 <- ggplotly(imR2plot, tooltip = c("Adjusted R<sup>2</sup>", "Method", "Cancer"))
    
    pR2$x$data[1:8] <- lapply(pR2$x$data[1:8], FUN = function(x){
      x$marker = list(opacity = 0)
      return(x)
    })
    
    pR2 %>%
      config(displayModeBar = F)
    
  })
  
  ## AIC
  
  output$imageAicFitBoxplot <- renderPlotly({
    
    imLong <- imageLongFitMetric %>% dplyr::distinct()
    
    imLong$Method <- gsub("MCPcounter", "MCP-Counter", imLong$Method)
    
    imLong$Method <- gsub("consensusTME_v2", "Consensus<sup>TME</sup>", imLong$Method)
    
    imMethodCols <- c("#ED2224", "#6FCCDD", "#F6EB16" , "#231F20",
                      "#087022", "#F7931D", "#B9529F", "#3953A4")
    
    names(imMethodCols) <- unique(imLong$Method)
    
    imMethodColScale <- scale_fill_manual(name = "Method ", values = imMethodCols)
    
    cancCols <- rep(c("#313131", "#776f7a", "#0436bf", "#3ee5ea", "#04601d",
                      "#0ce747", "#c0cc00", "#e5f656", "#f9a53d", "#760606",
                      "#f94e4e","#85037f","#bd3eea"), 8)
    
    imAicplot <- ggplot(data = imLong,
                       aes(x = reorder(Method, `AIC Z Score`, FUN = function(x) (median(x, na.rm = TRUE))),
                           y = `AIC Z Score`, labels = Method)) +
      geom_boxplot(fill = imMethodCols[levels(reorder(imLong$Method, imLong$`AIC Z Score`,
                                                      FUN = function(x) (median(x, na.rm = TRUE))))],
                   outlier.shape = NA) +
      geom_beeswarm(aes(fill = Cancer),
                    shape = 21,
                    alpha = 0.65,
                    cex = 1.75,
                    priority = "ascending"
      ) +
      theme_minimal() +
      xlab("Estimation Method") +
      theme(axis.text.x = element_text(angle = 45, vjust = 1.3, hjust = 1),
            panel.grid.minor = element_blank()
      ) 
    
    pAic <- ggplotly(imAicplot, tooltip = c("AIC Z Score", "Method", "Cancer"))
    
    pAic$x$data[1:8] <- lapply(pAic$x$data[1:8], FUN = function(x){
      x$marker = list(opacity = 0)
      return(x)
    })
    
    pAic %>%
      config(displayModeBar = F)
    
  })
  
  ## BIC
  
  output$imageBicFitBoxplot <- renderPlotly({
    
    imLong <- imageLongFitMetric %>% dplyr::distinct()
    
    imLong$Method <- gsub("MCPcounter", "MCP-Counter", imLong$Method)
    
    imLong$Method <- gsub("consensusTME_v2", "Consensus<sup>TME</sup>", imLong$Method)
    
    imMethodCols <- c("#ED2224", "#6FCCDD", "#F6EB16" , "#231F20",
                      "#087022", "#F7931D", "#B9529F", "#3953A4")
    
    names(imMethodCols) <- unique(imLong$Method)
    
    imMethodColScale <- scale_fill_manual(name = "Method ", values = imMethodCols)
    
    cancCols <- rep(c("#313131", "#776f7a", "#0436bf", "#3ee5ea", "#04601d",
                      "#0ce747", "#c0cc00", "#e5f656", "#f9a53d", "#760606",
                      "#f94e4e","#85037f","#bd3eea"), 8)
    
    imBicplot <- ggplot(data = imLong,
                        aes(x = reorder(Method, `BIC Z Score`, FUN = function(x) (median(x, na.rm = TRUE))),
                            y = `BIC Z Score`, labels = Method)) +
      geom_boxplot(fill = imMethodCols[levels(reorder(imLong$Method, imLong$`BIC Z Score`, FUN = median))],
                   outlier.shape = NA) +
      geom_beeswarm(aes(fill = Cancer),
                    shape = 21,
                    alpha = 0.65,
                    cex = 1.75,
                    priority = "ascending"
      ) +
      
      theme_minimal() +
      xlab("Estimation Method") +
      theme(axis.text.x = element_text(angle = 45, vjust = 1.3, hjust = 1),
            panel.grid.minor = element_blank()
      ) 
    
    pBic <- ggplotly(imBicplot, tooltip = c("BIC Z Score", "Method", "Cancer"))
    
    pBic$x$data[1:8] <- lapply(pBic$x$data[1:8], FUN = function(x){
      x$marker = list(opacity = 0)
      return(x)
    })
    
    pBic %>%
      config(displayModeBar = F)
    
  })
  
  
  #### MCP-Counter Benchmark ####
  
  load("./data/MCP_Benchmark_Data.RData")
  
  ## Scatter Plot
  output$mcpScatter <- renderPlotly({
    
    imMethodCols <- c("#3953A4", "#087022", "#B9529F", "#231F20",
                      "#F7931D", "#ED2224", "#F6EB16", "#6FCCDD")
    
    names(imMethodCols) <- unique(mcpBenchVals$Method)
    
    imMethodColScale <- scale_fill_manual(name = "Method ", values = imMethodCols)
    
    mcpScatter <- ggplot(mcpBenchVals, aes(RNA_Estimate,
                                           IHC_Estimate,
                                           fill = Method,
                                           text = paste0("Cell Type: ",
                                                         gsub("_"," ",Cell_Type),
                                                         "<br>Method: ", Method,
                                                         "<br> Correlation: ", format(round(Tau, 2), nsmall = 2),
                                                         "<br> Adjusted P-Value: ", formatC(signif(Q_Value,digits = 2), digits = 1, format = "e"))
    )
    ) +
      theme_bw() +
      imMethodColScale +
      theme(panel.spacing.y=unit(1, "lines"),
            plot.margin=unit(c(1,1,1.5,1.2),"cm"),
            legend.title=element_blank()) +
      geom_smooth(method = 'lm',
                  fill = "black",
                  se = FALSE,
                  linetype = 'dashed',
                  colour = 'black',
                  size = 0.5 ) +
      geom_point(shape = 21,
                 colour = 'black',
                 size = 2,
                 alpha = 0.8) +
      scale_x_continuous(breaks = c(0, 0.5, 1)) +
      xlab(" <br>RNA Estimates<br><sup>(Scaled)</sup>") +
      ylab("IHC Estimates (Cells/mm<sup>2</sup>)") +
      facet_grid(Cell_Type ~ Method_F , scale = 'free')
    
    ggplotly(mcpScatter,
             tooltip = "text") %>%
      layout(
        legend = list(
          orientation = "h",
          y = -0.25,
          x = 0.5
        )
      ) %>%
      config(displayModeBar = F)
    
  })
  
  ## Boxplots
  
  output$mcpBar <- renderPlotly({
    
    imMethodCols <- c("#3953A4", "#087022", "#B9529F" , "#231F20",
                      "#F7931D", "#ED2224", "#F6EB16", "#6FCCDD")
    
    names(imMethodCols) <- unique(mcpBenchVals$Method)
    
    imMethodColScale <- scale_fill_manual(name = "Method ", values = imMethodCols)
    
    mcpCorrs <- as.data.table(mcpCorrs)
    
    mcpCorrs[, ord := sprintf("%02i", frank(mcpCorrs, -Tau, ties.method = "first"))]
    
    mcpBar <- ggplot(data = mcpCorrs,
                     aes(x = ord,
                         y = Tau,
                         fill = Method,
                         text = paste0("Cell Type: ",
                                       gsub("_"," ",Cell_Type),
                                       "<br>Method: ", Method,
                                       "<br> Correlation: ", format(round(Tau, 2), nsmall = 2),
                                       "<br> Adjusted P-Value: ", formatC(signif(q_val, digits = 2), digits = 1, format = "e")
                                       )
                     )
    ) +
      geom_col(colour = "#333333",
               size = 0.25,
               width = 0.75) + 
      facet_grid(~ Cell_Type,
                 scales = "free_x",
                 drop = TRUE) +
      scale_x_discrete(labels = mcpCorrs[, setNames(as.character(Method), ord)]) + 
      theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=.5),
            panel.spacing.x=unit(2, "lines"),
            panel.background = element_rect(fill = "white", colour = "grey50"),
            plot.margin=unit(c(1,1,1.5,1.2),"cm")
            ) +
      xlab(" <br> <br> <br> Immune Estimation Methods (<i>Sorted</i>)") +
      ylab("Correlation Between RNA & IHC Estimates <br> <sup>(Kendall's &#964; )</sup>") +
      imMethodColScale
    
    ggplotly(mcpBar,
             tooltip = "text") %>%
      config(displayModeBar = F)
    
  })
  
  #### TIMER Boxplots ####
  
  timerLong <- readRDS("./data/TIMER_Benchmark_Data.rds")
  
  output$timerBoxplots <- renderPlotly({
    imMethodCols <- c("#3953A4", "#F6EB16", "#ED2224" , "#087022",
                      "#F7931D", "#6FCCDD", "#B9529F")
    
    names(imMethodCols) <- unique(timerLong$Method)
    
    imMethodColScale <- scale_fill_manual(name = "Method ", values = imMethodCols)
    
    timerPlot <- ggplot(data = timerLong,
                        aes (x = Pathologist_Category,
                             y = RNA_Z_Scored,
                             fill = Method,
                             text = paste0(Method,
                                           "<br><br>Tukey HSD Post-Hoc <i>p</i> values: ",
                                           "<br><br>Low vs Medium: ", `Low-Medium`,
                                           "<br>Medium vs High: ", `Medium-High`,
                                           "<br>Low vs High: ", `Low-High`
                             )
                        )
    ) +
      geom_boxplot(outlier.shape = NA,
                   outlier.size = NA,
                   outlier.stroke = NA,
                   notch = TRUE) +
      geom_quasirandom(varwidth = TRUE,
                       alpha = 0.5,
                       size = 0.85,
                       stroke = 0.25) +
      imMethodColScale +
      theme_minimal() +
      theme(panel.spacing.y=unit(1.2, "lines"),
            plot.margin=unit(c(0, 0, 0.5,1.2),"cm")) +
      ylab("RNA-Based Neutrophil Estimate<br><sup>(Z-Score)</sup>") +
      xlab("Pathologist-Based Neutrophil Quantity Category") +
      facet_wrap(Method~., ncol = 4,
                 drop = TRUE,
                 scales = "free_y")
    
    timPlotly <- ggplotly(timerPlot,
                          tooltip = "text")
    
    timPlotly$x$data[1:7] <- lapply(timPlotly$x$data[1:7], FUN = function(x){
      x$marker = list(opacity = 0)
      return(x)
    })
    
    
    
    timPlotly %>%
      layout(showlegend = F) %>%
      config(displayModeBar = F)
  })
  
  
  
  #### CIBERSORT Boxplots ####
  
  cibersortCorrs <- readRDS("./data/CIBERSORT_Benchmark_Correlations.rds")
  
  output$cibersortBoxplots <- renderPlotly({
    imMethodCols <- c("#F6EB16", "#6FCCDD", "#3953A4" , "#087022",
                      "#231F20", "#B9529F", "#ED2224", "#F7931D")
    
    names(imMethodCols) <- unique(cibersortCorrs$Method)
    
    imMethodColScale <- scale_fill_manual(name = "Method ", values = imMethodCols)
    
    cellCols = c("#800000", "#e6194b", "#f58231", "#914800", "#808000",
                 "#31b743", "#422100", "#004401", "#46f0f0", "#0082c8",
                 "#fabebe", "#e8e8e8")
    
    names(cellCols) <- c("B cells", "B cells memory", "B cells naive",
                         "Monocytes", "NK cells", "T cells CD4<sup>+</sup>",
                         "T cells CD4<sup>+</sup> memory resting",
                         "T cells CD4<sup>+</sup> memory activated", "T cells CD4<sup>+</sup> naive",
                         "T cells CD8<sup>+</sup>", "T cells gd", "T cells CD4<sup>+</sup> memory")
    
    imCellColScale <- scale_fill_manual(name = "Cell Type ", values = cellCols)
    
    imMethodCols<- imMethodCols[levels(reorder(cibersortCorrs$Method,
                                               cibersortCorrs$Tau,
                                               FUN = function(x) -(median(x, na.rm = TRUE))))]
    
    cibBox <- ggplot(data = cibersortCorrs,
                     aes(x = reorder(Method, Tau,
                                     FUN = function(x) -(median(x, na.rm = TRUE))),
                         y = Tau,
                         text = sprintf("Method: %s <br>Cell Type: %s <br>Correlation: %s <br>Adjusted P-Value: %s",
                                        Method,
                                        cibersortCorrs$`Cell Type`,
                                        format(round(Tau, 2), nsmall = 2),
                                        signif(q_val, digits = 2)
                                        )
                         )
                     ) +
      geom_beeswarm(aes(fill = `Cell Type`),
                    cex = 3,
                    alpha = 0.85,
                    size = 2.25
                    ) + 
      geom_boxplot(aes(x = reorder(Method, Tau,
                                   FUN = function(x) -(median(x, na.rm = TRUE))),
                       y = Tau),
                   fill = imMethodCols[levels(reorder(cibersortCorrs$Method,
                                                      cibersortCorrs$Tau,
                                                      FUN = function(x) -(median(x, na.rm = TRUE))))],
                   outlier.shape = NA,
                   alpha = 0.55,
                   coef = 3,
                   inherit.aes = FALSE) +
      theme_minimal() +
      xlab("Immune Estimation Method (<i>Sorted</i>)") +
      ylab("Correlation Between RNA & Flow Cytometry Estimates <br> <sup>(Kendall's &#964; )</sup>") +
      imCellColScale
    
    ggplotly(cibBox,
             tooltip = "text") %>%
      config(displayModeBar = F)
  })
  
  #### xCell Benchmarks ####
  
  load("./data/xCell_Benchmarking_Correlations.RData")
  
  output$xCellBoxplots311 <- renderPlotly({
    imMethodCols <- c("#F6EB16", "#6FCCDD", "#3953A4" , "#087022",
                      "#231F20", "#B9529F", "#ED2224", "#F7931D")
    
    names(imMethodCols) <- unique(xCellCorrs311$Method)
    
    imMethodColScale <- scale_fill_manual(name = "Method ", values = imMethodCols)
    
    cellCols = c("#800000", "#e6194b", "#f58231", "#914800", "#808000",
                 "#ffe119", "#d2f53c", "#31b743", "#008080", "#aaffc3",
                 "#46f0f0", "#0082c8", "#000080", "#911eb4", "#f032e6",
                 "#e6beff", "#fabebe", "#fffac8", "#808080", "#000000")
    
    names(cellCols) <- c("B cells", "B cells memory", "B cells naive",
                         "Monocytes", "NK cells", "Plasma cells", "T cells",
                         "T cells CD4<sup>+</sup>", "T cells CD4<sup>+</sup> centralmemory",
                         "T cells CD4<sup>+</sup> effector memory", "T cells CD4<sup>+</sup> naive",
                         "T cells CD8<sup>+</sup>", "T cells CD8<sup>+</sup> central memory",
                         "T cells CD8<sup>+</sup> effector",
                         "T cells CD8<sup>+</sup> effector memory","T cells CD8<sup>+</sup> naive",
                         "T cells Gd", "T cells central memory",
                         "T cells effector memory", "T regs")
    
    imCellColScale <- scale_fill_manual(name = "Cell Type ", values = cellCols)
    
    imMethodCols<- imMethodCols[levels(reorder(xCellCorrs311$Method,
                                               xCellCorrs311$Tau,
                                               FUN = function(x) -(median(x, na.rm = TRUE))))]
    
    xCellBox311 <- ggplot(data = xCellCorrs311,
                          aes(x = reorder(Method, Tau,
                                          FUN = function(x) -(median(x, na.rm = TRUE))),
                              y = Tau,
                              text = sprintf(paste0("Method: %s <br>Cell Type: %s <br>Correlation:",
                                                    "%s <br>Adjusted P-Value: %s"),
                                             Method,
                                             `Cell Type`,
                                             format(round(Tau, 2), nsmall = 2),
                                             signif(q_val, digits = 2)
                              )
                          )
    ) +
      geom_beeswarm(aes(fill = `Cell Type`),
                    cex = 2,
                    alpha = 0.85,
                    size = 2.25
      ) + 
      geom_boxplot(aes(x = reorder(Method, Tau,
                                   FUN = function(x) -(median(x, na.rm = TRUE))),
                       y = Tau),
                   fill = imMethodCols[levels(reorder(xCellCorrs311$Method,
                                                      xCellCorrs311$Tau,
                                                      FUN = function(x) -(median(x, na.rm = TRUE))))],
                   outlier.shape = NA,
                   alpha = 0.55,
                   coef = 3,
                   inherit.aes = FALSE) +
      theme_minimal() +
      theme(axis.title.y = element_text(margin = margin(10, 20, 0, 0))) +
      ggtitle("<b>SDY311</b>") +
      xlab("Immune Estimation Method (<i>Sorted</i>)") +
      ylab("Correlation Between RNA & CyTOF Estimates <br> <sup>(Kendall's &#964; )</sup>") +
      imCellColScale
    
    ggplotly(xCellBox311,
             tooltip = "text") %>%
      config(displayModeBar = F) %>%
      layout(margin = list(l = 100))
  })
  
  output$xCellBoxplots420 <- renderPlotly({
    imMethodCols <- c("#F6EB16", "#6FCCDD", "#3953A4" , "#087022",
                      "#231F20", "#B9529F", "#ED2224", "#F7931D")
    
    names(imMethodCols) <- unique(xCellCorrs420$Method)
    
    imMethodColScale <- scale_fill_manual(name = "Method ", values = imMethodCols)
    
    cellCols = c("#800000", "#e6194b", "#f58231", "#914800", "#808000",
                 "#ffe119", "#d2f53c", "#31b743", "#008080", "#aaffc3",
                 "#46f0f0", "#0082c8", "#000080", "#911eb4", "#f032e6",
                 "#e6beff", "#fabebe", "#fffac8", "#808080", "#000000")
    
    names(cellCols) <- c("B cells", "B cells memory", "B cells naive",
                         "Monocytes", "NK cells", "Plasma cells", "T cells",
                         "T cells CD4<sup>+</sup>", "T cells CD4<sup>+</sup> centralmemory",
                         "T cells CD4<sup>+</sup> effector memory", "T cells CD4<sup>+</sup> naive",
                         "T cells CD8<sup>+</sup>", "T cells CD8<sup>+</sup> central memory",
                         "T cells CD8<sup>+</sup> effector",
                         "T cells CD8<sup>+</sup> effector memory","T cells CD8<sup>+</sup> naive",
                         "T cells Gd", "T cells central memory",
                         "T cells effector memory", "T regs")
    
    imCellColScale <- scale_fill_manual(name = "Cell Type ", values = cellCols)
    
    imMethodCols<- imMethodCols[levels(reorder(xCellCorrs420$Method,
                                               xCellCorrs420$Tau,
                                               FUN = function(x) -(median(x, na.rm = TRUE))))]
    
    xCellBox420 <- ggplot(data = xCellCorrs420,
                          aes(x = reorder(Method, Tau,
                                          FUN = function(x) -(median(x, na.rm = TRUE))),
                              y = Tau,
                              text = sprintf(paste0("Method: %s <br>Cell Type: %s <br>Correlation:",
                                                    "%s <br>Adjusted P-Value: %s"),
                                             Method,
                                             `Cell Type`,
                                             format(round(Tau, 2), nsmall = 2),
                                             signif(q_val, digits = 2)
                              )
                          )
    ) +
      geom_beeswarm(aes(fill = `Cell Type`),
                    cex = 2,
                    alpha = 0.85,
                    size = 2.25
      ) + 
      geom_boxplot(aes(x = reorder(Method, Tau,
                                   FUN = function(x) -(median(x, na.rm = TRUE))),
                       y = Tau),
                   fill = imMethodCols[levels(reorder(xCellCorrs420$Method,
                                                      xCellCorrs420$Tau,
                                                      FUN = function(x) -(median(x, na.rm = TRUE))))],
                   outlier.shape = NA,
                   alpha = 0.55,
                   coef = 3,
                   inherit.aes = FALSE) +
      theme_minimal() +
      ggtitle("<b>SDY420</b>") +
      xlab("Immune Estimation Method (<i>Sorted</i>)") +
      ylab("Correlation Between RNA & CyTOF Estimates <br> <sup>(Kendall's &#964; )</sup>") +
      imCellColScale
    
    ggplotly(xCellBox420,
             tooltip = "text") %>%
      config(displayModeBar = F) %>%
      layout(margin = list(l = 100))
  })
  
  #### HGSOC Benchmark ####
  
  load("./data/HGSOC_Benchmarking.RData")
  
  output$hgsocBar <- renderPlotly({
    imMethodCols <- c("#F6EB16", "#6FCCDD", "#3953A4" , "#087022",
                      "#231F20", "#B9529F", "#ED2224", "#F7931D")
    
    names(imMethodCols) <- unique(hgsocCorrs$Method)
    
    imMethodColScale <- scale_fill_manual(name = "Method ", values = imMethodCols)
    
    hgsocCorrs <- as.data.table(hgsocCorrs)
    
    hgsocCorrs[, ord := sprintf("%02i", frank(hgsocCorrs, -Tau, ties.method = "first"))]
    
    hgsocBar <- ggplot(data = hgsocCorrs,
                       aes(x = ord,
                           y = Tau,
                           fill = Method,
                           text = paste0("Cell Type: ",
                                         gsub("_"," ",Cell_Type),
                                         "<br>Method: ", Method,
                                         "<br> Correlation: ", format(round(Tau, 2), nsmall = 2),
                                         "<br> Adjusted P-Value: ", signif(q_val, digits = 2)
                           )
                       )
    ) +
      geom_col(colour = "#333333",
               size = 0.25,
               width = 0.75) + 
      facet_grid(~ Cell_Type,
                 scales = "free_x",
                 drop = TRUE) +
      scale_x_discrete(labels = hgsocCorrs[, setNames(as.character(Method), ord)]) + 
      theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=.5),
            panel.spacing.x=unit(2, "lines"),
            panel.background = element_rect(fill = "white", colour = "grey50"),
            plot.margin=unit(c(1,1,1.5,1.2),"cm")
      ) +
      xlab(" <br> <br> <br> Immune Estimation Methods (<i>Sorted</i>)") +
      ylab("Correlation Between RNA & IHC Estimates <br> <sup>(Kendall's &#964; )</sup>") +
      imMethodColScale
    
    ggplotly(hgsocBar,
             tooltip = "text") %>%
      config(displayModeBar = F)
  })
  
  #### Results Overview ####
  
  
  ## TODO: Calculate ranks to allow new method to be ranked.
  
  ## Plot Line Graph Of Results
  
  benchRanks <- readRDS("./data/Method_Ranks.rds")
  
  benchRanksShared <- SharedData$new(benchRanks, ~Method)
  
  output$overviewLinePlot <- renderPlotly({
    
    imMethodCols <- c("#6FCCDD", "#3953A4", "#087022" , "#231F20",
                      "#ED2224", "#F6EB16", "#F7931D", "#B9529F")
    
    names(imMethodCols) <- unique(benchRanks$Method)
    
    imMethodColScale <- scale_colour_manual(name = "Method ", values = imMethodCols)
    
    rankLine <- ggplot(data = benchRanksShared,
                       aes(x = Benchmark,
                           y = Rank,
                           colour = Method,
                           text = paste0("Method: ", Method,
                                         # Return name of benchmarking experiment except for mean rank 
                                         sapply(Benchmark, function(bench) {
                                           if (bench == "<b>Mean Rank</b>") {
                                             return("<br><b>Mean Rank: </b>")
                                             } else {
                                               paste0("<br>Benchmark: ", gsub("<br>", " ", bench), "<br>Rank: ")
                                               }
                                           }),
                                         Rank)
                           )
                       ) +
      geom_line(aes(group = Method),
                size = 1,
                alpha = 1) +
      geom_vline(xintercept = 7.5, linetype="dashed", size = 0.1, alpha = 0.5) +
      geom_vline(xintercept = 10.5, linetype="dashed", size = 0.1, alpha = 0.5) +
      annotate(geom="text", x=4, y=0.5, label="TCGA Inferences",
               color="black") +
      annotate(geom="text", x=9, y=0.5, label="PBMCs",
               color="black") +
      annotate(geom="text", x=12, y=0.5, label="Bulk Tumour",
               color="black") +
      
      scale_y_continuous(trans = "reverse", breaks = 1:8) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45)) +
      xlab("Benchmarking Experiment") +
      imMethodColScale
    
    rankLineLy <- ggplotly(rankLine, tooltip = "text") %>%
      layout(showlegend = FALSE) %>%
      config(displayModeBar = F)
    highlight(rankLineLy, "plotly_hover",
              opacityDim = getOption("opacityDim", 0.1))
    
  })
  
  
  
  #### In Progress ####

}

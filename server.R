source('./useful_functions.R')

function(input, output, session) {
  session$onSessionEnded(function() {
    stopApp()
  })
  #### Run ConsensusTME In App ####
  
  output$rnaFile <- renderUI({fileInput(inputId = "rnaFile",
                                        "Choose bulk RNA-seq file (.txt or .csv):",
                                        accept = c(
                                          "text/plain",
                                          "text/csv"
                                        )
  )
  })
  
  # Set limit of upload file size to 400MB
  options(shiny.maxRequestSize=400*1024^2) 
  
  # Define reactive expressions
  rnaData <- reactiveVal(NULL)
  rnaCheckPass <- reactiveVal(FALSE)
  consMat <- reactiveVal(NULL)
  
  observeEvent(input$rnaFile, {
    
    # Stop if no file is selected
    inFile <- input$rnaFile
    if (is.null(inFile)) {
      return(NULL)
    }
    
    # Load data
    tryCatch({
      geneNames <- fread(inFile$datapath, header = T, select = 1)
      
      rnaMat <- as.matrix(fread(inFile$datapath, header = T, drop = 1))
      row.names(rnaMat) <- as.character(unlist(geneNames[,1]))
      rm(geneNames)
      gc(full =T)
      
      # Additional check to ensure the data is numerical
      if (!is.numeric(rnaMat)) {
        stop("Data should be a numeric matrix")
      }

      
      # Check for HUGO symbols by initially comparing against ConsensusTME gene sets
      allConsGenes <- unique(unlist(consensusGeneSets$Unfiltered))
      if (!any(row.names(rnaMat) %in% allConsGenes)) {
        stop("No genes from RNA rownames present in ConsensusTME gene sets\n Please ensure HUGO gene symbols are being used")
      }
      # Assign data to the reactive variable
      rnaData(rnaMat)
      rnaCheckPass(TRUE)
      rm(rnaMat)
      gc(full = T)
    }, error = function(e) {
      # Send sweet alert with error message
      sendSweetAlert(
        session = session,
        title = "Error in gene expression upload",
        text = e$message,
        type = "error"
      )
      # Clear the reactive variable
      rnaData(NULL)
      rnaCheckPass(FALSE)
    })
    
    # Check if less than 50% of necessary genes are present
    if (rnaCheckPass()) {
      date <- Sys.Date()
      formattedDate <- format(date, "%d%b%y")
      
      fileNameExt <- inFile$name
      fileNameNoExt <- sub("\\..*$", "", fileNameExt)
      # Initiate download handler
      output$downConsBttn <- downloadHandler(
        filename = function() {
          paste0("ConsensusTME_", formattedDate, "_", fileNameNoExt, ".csv")
        },
        content = function(con) {
          write.csv(consMat(), con)
        }
      )
      
      allConsGenes <- unique(unlist(consensusGeneSets$Unfiltered))
      genesPresent <- sum(allConsGenes %in% row.names(rnaData()))
      if (genesPresent < length(allConsGenes) * 0.5) {
        percPresent <- (genesPresent/length(allConsGenes)) * 100
        # Send sweet alert with warning
        confirmSweetAlert(
          session = session,
          inputId = "confirmLowGenes",
          title = "Warning",
          text = HTML(sprintf("Only %.1f%% of Consensus%s genes are present in the uploaded dataset. This may effect performance", percPresent, tags$sup("TME"))),
          html = T,
          type = "warning",
          btn_labels = c(HTML("Cancel <span style='font-size: 15px;'>(Upload new file)</span>"), "Proceed"),
          btn_colors = c("#ad423d", "#138f36"),
          cancelOnDismiss = TRUE,
          closeOnClickOutside = FALSE
        )
      } 
    }
  })
  
  # Once progression changes update
  observeEvent(input$confirmLowGenes, {
    if (is.null(input$confirmLowGenes)) {
      return(NULL)
    } else if (isFALSE(input$confirmLowGenes)) {
      rnaData(NULL)
      rnaCheckPass(FALSE)
      output$rnaFile <- renderUI({fileInput(inputId = "rnaFile",
                                            "Choose bulk RNA-seq file (.txt or .csv):",
                                            accept = c(
                                              "text/plain",
                                              "text/csv"
                                            )
      )
      })
    } else {
      return(NULL)
    }
  })
  
  # Once checks are complete render panel for running ConsensusTME
  
  output$consButtonPan <- renderUI({
    if(isTruthy(rnaCheckPass())){
      fluidRow(
        column(2,
               offset = 5,
               actionBttn(
                 inputId = "runConsBttn",
                 label = HTML(sprintf("Run Consensus%s", tags$sup("TME"))),
                 style = "unite",
                 color = "primary"
               )
        )
      )
      
    }
  })
  
  observeEvent(input$runConsBttn, {
    ## Check for cancer selection
    if(is.null(input$CancSelect)) {
      sendSweetAlert(
        session = session,
        title = "Cancer Type",
        text = "Cancer type not selected, for non-specific gene sets select Unfiltered",
        type = "error"
      )
      return(NULL)
    }
    shinybusy::show_modal_spinner(
      spin = "fingerprint",
      color = "#29a0f0",
      text = HTML(sprintf("Running Consensus%s", tags$sup("TME")))
    )
    
    if (FALSE) {
      consMat(consensusTMEAnalysis(rnaData(), cancerType = input$CancSelect))
      rnaData(NULL)
      gc(full = T)
    } else {
      # Number of chunks
      numChunks <- ceiling(ncol(rnaData()) / 5)
      
      runConsSplit <- function(mat, numChunks) {
        splitCols <- split(1:ncol(mat), cut(1:ncol(mat), numChunks))
        
        # Apply function to each chunk
        results <- lapply(splitCols, function(cols) {
          # Extract matrix chunk
          matChunk <- mat[, cols]

          #gsvaResult <- GSVA::gsva(expr = matChunk, gset.idx.list = consensusGeneSets[[input$CancSelect]], method = "ssgsea", ssgsea.norm = FALSE, parallel.sz=FALSE)
          gsvaResult <- cppssGSEAskinny(X = matChunk, geneSetList = ConsensusTME::consensusGeneSets[[input$CancSelect]])
          
          # Remove chunk from memory
          rm(matChunk)
          
          # Explicitly call garbage collection
          gc(full = T)
          
          return(gsvaResult)
        })
        
        # Combine results
        finalResult <- do.call(cbind, results)
        
        # Apply Barbie et al. 2009 normalisation 
        numCol <- ncol(finalResult)
        finalResult <- finalResult[, 1:numCol, drop=FALSE] / (range(finalResult)[2] - range(finalResult)[1])
        return(finalResult)
        
      }
      
      # Use the function
      consMat(runConsSplit(rnaData(), numChunks))
      unloadNamespace("GSVA")
      rnaData(NULL)
      gc(full = T)
    }
    
    shinybusy::remove_modal_spinner()
    sendSweetAlert(
      session = session,
      title ="Cell Type Estimates Generated!",
      type = "success"
    )
    
    if(isTruthy(rnaCheckPass())){
      output$consButtonPan <- renderUI({
        fluidRow(
          column(2,
                 offset = 3,
                 actionBttn(
                   inputId = "runConsBttn",
                   label = HTML(sprintf("Run Consensus%s", tags$sup("TME"))),
                   style = "unite",
                   color = "primary"
                 )
          ),
          column(3,
                 offset = 1,
                 shinyWidgets::downloadBttn(
                   outputId = "downConsBttn",
                   label = HTML(sprintf("Download Consensus%s Results", tags$sup("TME"))),
                   style = "unite",
                   color = "success"
                 )
          )
        )
      })
    }
  })
  
  observeEvent(consMat(), {
    if(!is.null(consMat())) {
      
      output$consHeat <- plotly::renderPlotly({
        consPlotMat <- consMat()
        
        ## Subset to 3dp for plotting only
        consPlotMat <- round(consPlotMat, 3)
        
        ## Move immune score to separate row colour
        immRow <- as.data.frame(consPlotMat["Immune_Score", ])
        colnames(immRow) <- "Immune Score"
        
        consPlotMat <- consPlotMat[row.names(consPlotMat) != "Immune_Score", ]
        row.names(consPlotMat) <- gsub("_", " ", row.names(consPlotMat))
        
        immScorePal <- wes_palette("Zissou1", nrow(immRow), type = "continuous")
        names(immScorePal) <- factor(sort(immRow$`Immune Score`))
        
        heatmaply(consPlotMat,
                  scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
                    low = "#762a83",
                    mid = "white",
                    high = "#1b7837",
                    midpoint = 0, 
                    limits = c(min(consPlotMat), max(consPlotMat))
                  ),
                  method = "plotly",
                  dendrogram = "both",
                  grid_gap = 0.5,
                  branches_lwd = 0.3,
                  subplot_heights = c(0.15, 0.05, 0.8),
                  midpoint = 0,
                  col_side_colors = immRow,
                  col_side_palette = immScorePal,
                  key.title = "Consensus<sup>TME</sup> NES",
                  dend_hoverinfo = FALSE,
                  ylab = "Cell Types",
                  xlab = "Bulk Samples \n <sup> (hover for details) </sup>",
                  showticklabels = c(F,T),
                  label_names = c("Immune Cell", "Sample ID", "NES"),
                  heatmap_layers = theme(axis.line = element_blank())) %>%
          layout(height=700) %>%
          config(displayModeBar = T) %>%
          layout(paper_bgcolor='transparent') %>%
          layout(plot_bgcolor='transparent')
      })
      
      
      output$consTab <- renderFormattable({
        consTabMat <- as.data.frame(t(consMat()))
        consTabMat <- round(consTabMat, 2)
        colnames(consTabMat) <- gsub("_", " ", colnames(consTabMat))
        formattable(consTabMat, list(area(col = colnames(consTabMat)) ~ normalize_bar("#fcac6f")))
        
      })
      output$exploreRes <- renderUI({
        wellPanel(
          tags$style(type="text", "#string {text-align:center}"),
          HTML("<h4> 2. Explore cell type estimates  </h4>"),
          tags$hr(),
          tags$br(),
          tabsetPanel(
            tabPanel("Heatmap",
                     tags$br(),
                     withSpinner(
                       plotly::plotlyOutput("consHeat", height = 700),
                       type = 6),
                     tags$hr()
                     ),
            tabPanel("Table",
                     withSpinner(
                     formattableOutput("consTab", height = 700),
                     type = 6)
            )
          )
        )
      })
    }
  })
  
  
  
  
  #### Tumour Purity Benchmark ####
  
  load("./data/Tumour_Purity_Corrs.RData")
  load("./data/TCGA_Annotation.RData")
  leukSideCols <- side_cols[!row.names(side_cols) %in% c("THYM", "DLBC"), ]
  
  ## Plot Heatmap 
  output$purCorrPlotly <- plotly::renderPlotly({
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
  output$purCorrBoxPlotly <- plotly::renderPlotly({
    purLongCorrs <- longCorrs
    purLongCorrs <- purLongCorrs %>% dplyr::distinct()
    purLongCorrs$Method <- gsub("_Immune_Score", "", purLongCorrs$Method)
    purLongCorrs$Method <- gsub("MCP", "MCP-Counter", purLongCorrs$Method)
    purLongCorrs$Method <- gsub("ConsensusTME", "Consensus<sup>TME</sup>", purLongCorrs$Method)
    
    methodCols <- c("#F6EB16", "#6FCCDD", "#3953A4" , "#087022",
                     "#231F20", "#A65628", "#B9529F", "#ED2224", "#F7931D")
    
    names(methodCols) <- unique(purLongCorrs$Method)
    
    colVector = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666", "#1B9E77",
                  "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#1F78B4",
                  "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99",
                  "#B15928", "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC",
                  "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC",
                  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999",
                  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7",
                  "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD","#CCEBC5",
                  "#FFED6F")
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
    
    p <- plotly::ggplotly(plot, tooltip = c("Tau", "Method", "Cancer"))
    
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
  
  output$leukR2FitHeatmap <- plotly::renderPlotly({
    
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
  
  output$leukAicFitHeatmap <- plotly::renderPlotly({
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
  
  output$leukBicFitHeatmap <- plotly::renderPlotly({
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
  output$leukR2FitBoxplot <- plotly::renderPlotly({
    
    leukLong <- longFitMetric %>% dplyr::distinct()
    
    leukLong$Method <- gsub("MCPcounter", "MCP-Counter", leukLong$Method)
    
    leukLong$Method <- gsub("consensusTME_v2", "Consensus<sup>TME</sup>", leukLong$Method)
    
    leukMethodCols <- c("#ED2224", "#6FCCDD", "#F6EB16" , "#231F20",
                        "#087022", "#F7931D", "#B9529F", "#3953A4")
    
    names(leukMethodCols) <- unique(leukLong$Method)
    
    colVector = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666", "#1B9E77",
                           "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#1F78B4",
                           "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99",
                           "#B15928", "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC",
                           "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC",
                           "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999",
                           "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7",
                           "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD","#CCEBC5",
                           "#FFED6F")
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
    
    pR2 <- plotly::ggplotly(leukR2plot, tooltip = c("Adjusted R<sup>2</sup>", "Method", "Cancer"))
    
    pR2$x$data[1:8] <- lapply(pR2$x$data[1:8], FUN = function(x){
      x$marker = list(opacity = 0)
      return(x)
    })
    
    pR2 %>%
      config(displayModeBar = F)
  })
  
  ## AIC
  output$leukAicFitBoxplot <- plotly::renderPlotly({
    
    leukLong <- longFitMetric %>% dplyr::distinct()
    
    leukLong$Method <- gsub("MCPcounter", "MCP-Counter", leukLong$Method)
    
    leukLong$Method <- gsub("consensusTME_v2", "Consensus<sup>TME</sup>", leukLong$Method)
    
    leukMethodCols <- c("#ED2224", "#6FCCDD", "#F6EB16" , "#231F20",
                        "#087022", "#F7931D", "#B9529F", "#3953A4")
    
    names(leukMethodCols) <- unique(leukLong$Method)
    
    colVector = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666", "#1B9E77",
                           "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#1F78B4",
                           "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99",
                           "#B15928", "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC",
                           "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC",
                           "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999",
                           "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7",
                           "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD","#CCEBC5",
                           "#FFED6F")
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
    
    pAic <- plotly::ggplotly(leukAicplot, tooltip = c("AIC Z Score", "Method", "Cancer"))
    
    pAic$x$data[1:8] <- lapply(pAic$x$data[1:8], FUN = function(x){
      x$marker = list(opacity = 0)
      return(x)
    })
    
    pAic %>%
      config(displayModeBar = F)
  })
  
  ## BIC
  
  output$leukBicFitBoxplot <- plotly::renderPlotly({
    
    leukLong <- longFitMetric %>% dplyr::distinct()
    
    leukLong$Method <- gsub("MCPcounter", "MCP-Counter", leukLong$Method)
    
    leukLong$Method <- gsub("consensusTME_v2", "Consensus<sup>TME</sup>", leukLong$Method)
    
    
    leukMethodCols <- c("#ED2224", "#6FCCDD", "#F6EB16" , "#231F20",
                        "#087022", "#F7931D", "#B9529F", "#3953A4")
    
    names(leukMethodCols) <- unique(leukLong$Method)
    
    colVector = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666", "#1B9E77",
                           "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#1F78B4",
                           "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99",
                           "#B15928", "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC",
                           "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC",
                           "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999",
                           "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7",
                           "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD","#CCEBC5",
                           "#FFED6F")
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
    
    pBic <- plotly::ggplotly(leukBicplot, tooltip = c("BIC Z Score", "Method", "Cancer"))
    
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
  
  output$imageR2FitHeatmap <- plotly::renderPlotly({
    
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
  
  output$imageAicFitHeatmap <- plotly::renderPlotly({
    
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
  
  output$imageBicFitHeatmap <- plotly::renderPlotly({
    
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
  
  output$imageR2FitBoxplot <- plotly::renderPlotly({
    
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
    
    pR2 <- plotly::ggplotly(imR2plot, tooltip = c("Adjusted R<sup>2</sup>", "Method", "Cancer"))
    
    pR2$x$data[1:8] <- lapply(pR2$x$data[1:8], FUN = function(x){
      x$marker = list(opacity = 0)
      return(x)
    })
    
    pR2 %>%
      config(displayModeBar = F)
    
  })
  
  ## AIC
  
  output$imageAicFitBoxplot <- plotly::renderPlotly({
    
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
    
    pAic <- plotly::ggplotly(imAicplot, tooltip = c("AIC Z Score", "Method", "Cancer"))
    
    pAic$x$data[1:8] <- lapply(pAic$x$data[1:8], FUN = function(x){
      x$marker = list(opacity = 0)
      return(x)
    })
    
    pAic %>%
      config(displayModeBar = F)
    
  })
  
  ## BIC
  
  output$imageBicFitBoxplot <- plotly::renderPlotly({
    
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
    
    pBic <- plotly::ggplotly(imBicplot, tooltip = c("BIC Z Score", "Method", "Cancer"))
    
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
  output$mcpScatter <- plotly::renderPlotly({
    
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
    
    plotly::ggplotly(mcpScatter,
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
  
  output$mcpBar <- plotly::renderPlotly({
    
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
    
    plotly::ggplotly(mcpBar,
             tooltip = "text") %>%
      config(displayModeBar = F)
    
  })
  
  #### TIMER Boxplots ####
  
  timerLong <- readRDS("./data/TIMER_Benchmark_Data.rds")
  
  output$timerBoxplots <- plotly::renderPlotly({
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
    
    timPlotly <- plotly::ggplotly(timerPlot,
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
  
  output$cibersortBoxplots <- plotly::renderPlotly({
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
    
    plotly::ggplotly(cibBox,
             tooltip = "text") %>%
      config(displayModeBar = F)
  })
  
  #### xCell Benchmarks ####
  
  load("./data/xCell_Benchmarking_Correlations.RData")
  
  output$xCellBoxplots311 <- plotly::renderPlotly({
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
    
    plotly::ggplotly(xCellBox311,
             tooltip = "text") %>%
      config(displayModeBar = F) %>%
      layout(margin = list(l = 100))
  })
  
  output$xCellBoxplots420 <- plotly::renderPlotly({
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
    
    plotly::ggplotly(xCellBox420,
             tooltip = "text") %>%
      config(displayModeBar = F) %>%
      layout(margin = list(l = 100))
  })
  
  #### HGSOC Benchmark ####
  
  load("./data/HGSOC_Benchmarking.RData")
  
  output$hgsocBar <- plotly::renderPlotly({
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
                 scales = "free",
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
    
    plotly::ggplotly(hgsocBar,
             tooltip = "text") %>%
      config(displayModeBar = F)
  })
  
  #### Results Overview ####
  
  
  ## TODO: Calculate ranks to allow new method to be ranked.
  
  ## Plot Line Graph Of Results
  
  benchRanks <- readRDS("./data/Method_Ranks.rds")
  
  benchRanksShared <- SharedData$new(benchRanks, ~Method)
  
  output$overviewLinePlot <- plotly::renderPlotly({
    
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
    
    rankLineLy <- plotly::ggplotly(rankLine, tooltip = "text") %>%
      layout(showlegend = FALSE) %>%
      config(displayModeBar = F)
    highlight(rankLineLy, "plotly_hover",
              opacityDim = getOption("opacityDim", 0.1))
    
  })
  
  #### Current Approaches Table ####
  currDf <- read.csv("./data/Deconvolution_Methods.csv", fileEncoding = "UTF-8", check.names = F, stringsAsFactors = T)
  
  currDf$`Publication Link` <- as.character(currDf$`Publication Link`)
  
  currDf$`Link to Tool` <- mapply(function(datLink, extraInf){
    datLink <- as.character(datLink)
    extraInf <- paste0("- ", as.character(extraInf))
    extraInf <- gsub("Dead Link -", "Dead Link <br><sub>", extraInf)
    
    if (datLink == "N/A") {
      return(sprintf('<i class="fa fa-exclamation-triangle" aria-hidden="true" title="No link available"></i> %s</sub>', extraInf))
    } else {
      linkSingle = sprintf('<a href="%s" target="_blank" title="%s"><i class="fa fa-external-link" aria-hidden="true"></i></a> %s', datLink, datLink, extraInf)
      return(linkSingle)
    }
  }, currDf$`Link to Tool`, currDf$`Tool Info`, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  currDf$`Publication Link` <- sprintf('<a href="%s" target="_blank"><i class="fa fa-external-link" aria-hidden="true"></i></a>', currDf$`Publication Link`)
  currDf$`Tool Info` <- NULL
  
  output$currAppTab <- renderDT({
    datatable(
      currDf,
      escape = FALSE,  # interpret the content in the column as HTML
      selection = "none",
      filter = 'top',
      rownames = FALSE,
      options = list(
        autoWidth = TRUE,
        pageLength = 10,
        fixedHeader = FALSE, 
        scrollX = FALSE,
        scrollY = '1000px',
        scrollCollapse = TRUE,
        paging = FALSE,
        columnDefs = list(list(searchable = FALSE, targets = c(6, 7)))
      )
    )
  })
  
  #### Benchmarking Datasets Tables ####
  
  benchResDf <- read.csv("./data/Benchmarking_Resources_Data.csv", fileEncoding = "UTF-8", check.names = F, stringsAsFactors = T)
  
  benchResDf$`Data Link` <- mapply(function(datLink, extraInf){
    datLink <- as.character(datLink)
    extraInf <- as.character(extraInf)
    
    if (extraInf == "Yes") {
      extraInf <- ""
    } else {
      extraInf <- paste0(" - ", extraInf)
    }
    
    if (datLink == "N/A") {
      return(sprintf('<i class="fa fa-exclamation-triangle" aria-hidden="true" title="No link available"></i> %s', extraInf))
    } else {
      if (grepl(";", datLink)) {
        # Split multiple links
        links <- strsplit(datLink, ";")[[1]]
        # Generate HTML for each link
        linkStrings <- sapply(links, function(link) {
          sprintf('<a href="%s" target="_blank" title="%s"><i class="fa fa-external-link" aria-hidden="true"></i></a>', link, link)
        })
        # Combine back into a single string
        linkCombined = paste(linkStrings, collapse = " ")
        return(paste0(linkCombined, extraInf))
      } else {
        linkSingle = sprintf('<a href="%s" target="_blank" title="%s"><i class="fa fa-external-link" aria-hidden="true"></i></a>', datLink, datLink)
        return(paste0(linkSingle, extraInf))
      }
    }
  }, benchResDf$`Data Link`, benchResDf$`Available Publicly`, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  benchResDf$`Available Publicly` <- NULL
  
  benchResDf$`Publication Link` <- mapply(function(pubLink, pubTitle){
    pubLink <- as.character(pubLink)
    pubTitle <- as.character(pubTitle)
    
    if (pubLink == "N/A") {
      pubTitle <- paste0(" - ", pubTitle)
      return(sprintf('<i class="fa fa-exclamation-triangle" aria-hidden="true" title="No link available"></i> %s', pubTitle))
    } else {
      if (grepl(";", pubLink)) {
        # Split multiple links
        links <- strsplit(pubLink, ";")[[1]]
        # Generate HTML for each link
        linkStrings <- sapply(links, function(link) {
          sprintf('<a href="%s" target="_blank" title="%s"><i class="fa fa-external-link" aria-hidden="true"></i></a>', link, pubTitle)
        })
        # Combine back into a single string
        linkCombined = paste(linkStrings, collapse = " ")
        return(linkCombined)
      } else {
        linkSingle = sprintf('<a href="%s" target="_blank" title="%s"><i class="fa fa-external-link" aria-hidden="true"></i></a>', pubLink, pubTitle)
        return(linkSingle)
      }
    }
  }, benchResDf$`Publication Link`, benchResDf$`Publication`, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  benchResDf$`Publication` <- NULL
  
  output$benchDsTab <- renderDT({
    datatable(
      benchResDf,
      escape = FALSE, 
      selection = "none",
      filter = 'top',
      rownames = FALSE,
      options = list(
        autoWidth = TRUE,
        pageLength = 10,
        fixedHeader = TRUE, 
        scrollY = '1000px',
        scrollX = FALSE,
        scrollCollapse = TRUE,
        paging = FALSE,
        columnDefs = list(list(searchable = FALSE, targets = c(8, 9)))
      )
    )
  })
  
  #### Review Article Tables ####
  
  revArtDf <- read.csv("./data/Review_Articles.csv", fileEncoding = "UTF-8", check.names = F, stringsAsFactors = T)
  
  revArtDf$`Publication Link` <- mapply(function(pubLink, pubTitle){
    pubTitle <- as.character(pubTitle)
    pubLink <- as.character(pubLink)
    
    return(sprintf('<a href="%s" target="_blank" title="%s"><i class="fa fa-external-link" aria-hidden="true"></i></a>', pubLink, pubTitle))
    
  }, revArtDf$`Publication Link`, revArtDf$Title, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  output$reviewTab <- renderDT({
    datatable(
      revArtDf,
      escape = FALSE,
      selection = "none",
      filter = 'top',
      rownames = FALSE,
      options = list(
        autoWidth = FALSE,
        pageLength = 10,
        fixedHeader = TRUE, 
        scrollY = '1000px',
        scrollX = FALSE,
        scrollCollapse = TRUE,
        paging = FALSE,
        columnDefs = list(list(searchable = FALSE, targets = 3))
      )
    )
  })
  
}

#### Useful Functions ####

removeBlanks <- function(list.in) {
# Removes blanks from vector of character or a list containing gene sets
#  
# Args:
#   list.in: Either character vector or list containing gene sets
#
# Returns:
#   The vector or list with any blanks removed
  
  if (class(list.in) == 'list') {
    list.out <- lapply(list.in[list.in != ''], function(x) x[x != ''])
    return(list.out)
  } else if (class(list.in) == 'character') {
    list.out <- list.in[list.in != '']
    return(list.out)
  } else {
    sprintf('%s is not a class this function can use', class(list.in))
  }
}

mergeDups <- function(list.in, removenames = TRUE) {
# When supplied with a list with duplicate names, merges elements of list
# with an identical name.
# 
# Args:
#  list.in: List containing elements with the same name
#  removenames: if true removes the name of each string within list elements
# Returns:
#  A list with the elements with the same name combined

    if (class(list.in) == 'list'){
    lnames <- names(list.in)
    ulnames <- unique(lnames)
    list.grouped <- lapply(ulnames, function(x) list.in[lnames == x])
    list.merged <- lapply(list.grouped, unlist, recursive = F)
    if (removenames) {
      list.merged <- lapply(list.merged, unname)
    }
    names(list.merged) <- ulnames
    return (list.merged)
  } else {
    sprintf('%s is not a class this function can currently use',class(list.in))
  }
}
# 
# radioTooltip <- function(id, choice, title, placement = "bottom", trigger = "hover", options = NULL){
# 
#   options = shinyBS:::buildTooltipOrPopoverOptionsList(title, placement, trigger, options)
#   options = paste0("{'", paste(names(options), options, sep = "': '", collapse = "', '"), "'}")
#   bsTag <- shiny::tags$script(shiny::HTML(paste0("
#                                                  $(document).ready(function() {
#                                                  setTimeout(function() {
#                                                  $('input', $('#", id, "')).each(function(){
#                                                  if(this.getAttribute('value') == '", choice, "') {
#                                                  opts = $.extend(", options, ", {html: true});
#                                                  $(this.parentElement).tooltip('destroy');
#                                                  $(this.parentElement).tooltip(opts);
#                                                  }
#                                                  })
#                                                  }, 500)
#                                                  });
#                                                  ")))
#   htmltools::attachDependencies(bsTag, shinyBS:::shinyBSDep)
# }



Rcpp_fasterRndWalk <- function(n, R, Ra, X, geneSets) {
  es <- lapply(as.list(1:n), function(j) {
    geneRanking <- order(R[, j], decreasing=TRUE)
    es_sample <- lapply(geneSets, fasterRndWalk, geneRanking, j, Ra)
    unlist(es_sample)
  })
  es <- do.call("cbind", es)
  return(es)
}

# cppssGSEA <- function(geneMat, geneSetList) {
#   X <- GSVA:::.filterFeatures(as.matrix(geneMat), "ssgsea")
#   n <- ncol(X)
#   R <- t(sparseMatrixStats::colRanks(X, ties.method = "average"))
#   mode(R) <- "integer"
#   Ra <- abs(R)^0.25
# 
#   geneSets2 <- GSVA:::.mapGeneSetsToFeatures(geneSetList, rownames(X))
#   nesMat <- Rcpp_fasterRndWalk(n, R, Ra, X, geneSets)
#   colnames(nesMat) <- colnames(geneMat)
#   return(nesMat)
# }

cppssGSEAskinny <- function(X, geneSetList) {
  n <- ncol(X)
  R <- t(sparseMatrixStats::colRanks(X, ties.method = "average"))
  mode(R) <- "integer"
  Ra <- abs(R)^0.25
  geneRowNames <- seq_along(row.names(X))
  names(geneRowNames) <- row.names(X)
  geneSets <- lapply(geneSetList, function(geneGroup) {
    rowNums <- geneRowNames[geneGroup]
    rowNums <- rowNums[!is.na(rowNums)]
    as.integer(rowNums)
  })
  nesMat <- Rcpp_fasterRndWalk(n, R, Ra, X, geneSets)
  colnames(nesMat) <- colnames(X)
  return(nesMat)
}














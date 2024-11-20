#' TESet class definition
#' @export
setClass(
  "TESet",
  slots = c(
    counts = "matrix",
    rowData = "DataFrame",
    colData = "DataFrame",
    metadata = "list"
  ),
  prototype = list(
    counts = matrix(0, 0, 0),
    rowData = S4Vectors::DataFrame(),
    colData = S4Vectors::DataFrame(),
    metadata = list()
  )
)

#' Create TESet object
#'
#' @param counts matrix Count matrix
#' @param rowData DataFrame Row annotations
#' @param colData DataFrame Column annotations
#' @param metadata list Additional metadata
#' @return TESet object
#' @export
createTESet <- function(
    counts,
    rowData = NULL,
    colData = NULL,
    metadata = list()
) {
  # Create default rowData if not provided
  if (is.null(rowData)) {
    rowData <- S4Vectors::DataFrame(row.names = rownames(counts))
  }
  
  # Create default colData if not provided
  if (is.null(colData)) {
    colData <- S4Vectors::DataFrame(row.names = colnames(counts))
  }
  
  # Create TESet object
  te_set <- new(
    "TESet",
    counts = counts,
    rowData = rowData,
    colData = colData,
    metadata = metadata
  )
  
  return(te_set)
}

#' Show method for TESet objects
#' @param object TESet object
setMethod(
  "show",
  "TESet",
  function(object) {
    cat("TESet object with", nrow(object@counts), "features and",
        ncol(object@counts), "cells\n")
    cat("rowData:", ncol(object@rowData), "columns\n")
    cat("colData:", ncol(object@colData), "columns\n")
    cat("metadata:", length(object@metadata), "elements\n")
  }
)
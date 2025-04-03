#' Generates sparse form matrix for ASReml-R
#'
#' \code{full2sparse} Generates the three-column form matrix order as a lower triangule row-wise,
#' where the first two columns are the row and column indices and the third holds the values.
#' 
#' @param X square matrix 
#' @param rowNames optional input with names of individuals
#' 
#' @references
#' Adapted from function published by Borgognone et al. (2016). Crop Science 56:2616-2628

full2sparse <- function (X, rowNames=NA) {
  
  if (length(rowNames)==1) { rowNames = dimnames(X)[[1]] }
  
  which <- (X != 0 & lower.tri(X, diag = TRUE))
  df <- data.frame(Row = t(row(X))[t(which)], Col = t(col(X))[t(which)],
                   Value = t(X)[t(which)])
  if (is.null(rowNames)) {
    rowNames <- as.character(1:nrow(X))
  } else {
    rowNames <- as.character(rowNames)
  }
  
  attr(df, "rowNames") <- rowNames
  df
  
}
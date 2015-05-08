#' Potentials
#' @param dfg     discrete factor graph object
#' @return A list of current factors
#'
#'
#'
potentials <- function(dfg){
    stopifnot(is.dfg(dfg))

    dfg$dfgmodule$getFactorPotentials()
}

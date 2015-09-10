#################################################
# Building custom optimization functions
#################################################

update.matrix <- function(pot, expCounts){
  sweep( expCounts, 1, STATS = rowSums(expCounts), FUN = '/')
}

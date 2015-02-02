library(PGMscore)
context("Tail approximations")

test_that("NA when out of range saddlepoint",{
  varDim <- 2L
  facPot <- c(list(matrix(c(0.5,0.5),1,2)) )
  facNbs <- c(list(c(1L)) )
  
  mydfg <- dfg(varDim, facPot, facNbs)
  
  #Foreground model
  facPotFg <- c(list(matrix(c(0.9,0.1),1,2)))
  
  #Score range [log(0.1/0.5);log(0.9/0.5)] = [-1.609, .588]
  dfsaddle <- tailSaddle( seq(-2,1,0.1), mydfg, facPotFg)
  expect_equal( which(is.na(dfsaddle$p)), c(1:4,27:31) )
})
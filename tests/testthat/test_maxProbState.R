library(PGMscore)
context("Max probable state")

test_that("Independent variables cpp level",{
  varDim <- rep(4,4)
  facPot <- c(list(matrix(c(0.05,0.05,0.20,0.70),1,4)),
              list(matrix(c(0.05,0.05,0.70,0.20),1,4)),
              list(matrix(c(0.70,0.20,0.05,0.05),1,4)),
              list(matrix(c(0.05,0.70,0.20,0.05),1,4)))
  facNbs <- c(list(c(1L)),
              list(c(2L)),
              list(c(3L)),
              list(c(4L)))
  
  mydfg <- dfg(varDim, facPot, facNbs)
  
  mps <- mydfg$dfgmodule$maxProbState(list(),integer(0), logical(0))
  expect_equal( mps, c(3,2,0,1)) 
})

test_that("Independent variables cpp level with observed variables",{
  varDim <- rep(4,4)
  facPot <- c(list(matrix(c(0.05,0.05,0.20,0.70),1,4)),
              list(matrix(c(0.05,0.05,0.70,0.20),1,4)),
              list(matrix(c(0.70,0.20,0.05,0.05),1,4)),
              list(matrix(c(0.05,0.70,0.20,0.05),1,4)))
  facNbs <- c(list(c(1L)),
              list(c(2L)),
              list(c(3L)),
              list(c(4L)))
  
  mydfg2 <- dfg(varDim, facPot, facNbs)
  mps <- mydfg2$dfgmodule$maxProbState(list(),c(0,0,1,0), c(F,F,T,F))
  expect_equal( mps, c(3,2,1,1) )
})
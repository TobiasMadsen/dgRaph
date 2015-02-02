library(PGMscore)
context("Factor Marginals")

test_that("Factor marginals in disconnected graph",{
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

  #Single observation
  facExp <- mydfg$dfgmodule$facExpCounts(list(), matrix(c(NA,1,2,3),1,4) )
  
  expect_equal( as.vector(facExp[[1]]), c(0.05,0.05,0.20,0.70) )
  expect_equal( as.vector(facExp[[2]]), c(0,1,0,0))
  
  #Multiple observations
  facExp <- mydfg$dfgmodule$facExpCounts(list(), matrix(c(NA,1,2,3,
                                                          2,NA,NA,NA), 2, 4, byrow=T) )
  
  expect_equal( as.vector(facExp[[1]]), c(0.05,0.05,0.20,0.70) + c(0,0,1,0) )
  expect_equal( as.vector(facExp[[4]]), c(0,0,0,1) + c(0.05,0.70,0.20,0.05) )
})

test_that("Factor marginals in connected graph",{
  varDim <- rep(2L, 2)
  facPot <- c(list(matrix(c(0.5, 0.5),1,2)),
              list(matrix(c(0.75,0.25,0.25,0.75),2,2)))
  facNbs <- c(list(c(1L)),
              list(c(1L,2L)) )
  
  mydfg <- dfg(varDim, facPot, facNbs)
  
  facExp <- mydfg$dfgmodule$facExpCounts(list(), matrix(c(NA,1), 1, 2) )
  expect_equal( as.vector(facExp[[1]]), c(0.25,0.75))
  expect_equal( facExp[[2]], matrix(c(0,0,0.25,0.75),2,2))
  
  facExp <- mydfg$dfgmodule$facExpCounts(list(), matrix(c(NA,1,
                                                          0,NA), 2, 2, byrow=T) )
  expect_equal( as.vector(facExp[[1]]), c(0.25,0.75)+c(1,0))
  expect_equal( facExp[[2]], matrix(c(0,0,0.25,0.75),2,2) + matrix(c(0.75,0,0.25,0),2,2) )
})
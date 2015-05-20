library(dgRaph)
context("Factor Marginals")

test_that("Factor Marginals Interface 1",{
  source("cases/fourIndependentVariables.R")
  mydfg <- fourIndependentVariables()
  
  # Pass a matrix
  data <- matrix(c(NA,2,3,4,
                   3,NA,NA,NA), 2, 4, byrow=T) 
  expCounts <- facExpectations(data, mydfg)
  expect_equal( expCounts[[1]], matrix(c(0.05,0.05,1.2,0.7),1,4) )
  expect_equal( expCounts[[2]], matrix(c(0.05,1.05,0.7,0.2),1,4))
  
  # Pass a dataframe
  expCounts <- facExpectations(as.data.frame(data), mydfg)
  expect_equal( expCounts[[1]], matrix(c(0.05,0.05,1.2,0.7),1,4) )
})

test_that("Factor Marginals Interface 2",{
  source("cases/twoDepedentVariables.R")
  mydfg <- twoDependentVariables()
  
  # Pass a dataframe
  data <- data.frame(O1 = c(NA, 1), O2 = c(2, NA))
  expCounts <- facExpectations(data, mydfg)
  expect_equal( expCounts[[1]], matrix(c(1.4375,0.5625),1,2) )
  expect_equal( expCounts[[2]], matrix(c(0.75, 0, 0.6875, 0.5625),2,2))
})

test_that("Factor Marginals Interface 3 Shared Potentials",{
  # IID variables
  source("cases/IIDVariables.R")
  mydfg <- IIDVariables(10)
  
  # Pass a dataframe
  expCounts <- facExpectations(data.frame(matrix(c(1,1,1,1,2,2,2,NA,1,1),1,10)), mydfg)
  expect_equal( expCounts[[1]], matrix(c(6.7, 3.3), 1, 2))  
})

test_that("Factor Marginals C++ level 1",{
  source("cases/fourIndependentVariables.R")
  mydfg <- fourIndependentVariables()

  #Single observation
  facExp <- mydfg$dfgmodule$facExpCounts(matrix(c(NA,2,3,4),1,4) )
  
  expect_equal( as.vector(facExp[[1]]), c(0.05,0.05,0.20,0.70) )
  expect_equal( as.vector(facExp[[2]]), c(0,1,0,0))
  
  #Multiple observations
  facExp <- mydfg$dfgmodule$facExpCounts(matrix(c(NA,2,3,4,
                                                  3,NA,NA,NA), 2, 4, byrow=T) )
  
  expect_equal( as.vector(facExp[[1]]), c(0.05,0.05,0.20,0.70) + c(0,0,1,0) )
  expect_equal( as.vector(facExp[[4]]), c(0,0,0,1) + c(0.05,0.70,0.20,0.05) )
})

test_that("Factor Marginals C++ level 2",{
  source("cases/twoDepedentVariables.R")
  mydfg <- twoDependentVariables()
  
  facExp <- mydfg$dfgmodule$facExpCounts(matrix(c(NA,2), 1, 2) )
  expect_equal( as.vector(facExp[[1]]), c(0.4375,0.5625))
  expect_equal( facExp[[2]], matrix(c(0,0,0.4375,0.5625),2,2))
  
  facExp <- mydfg$dfgmodule$facExpCounts(matrix(c(NA,2,
                                                  1,NA), 2, 2, byrow=T) )
  expect_equal( as.vector(facExp[[1]]), c(0.4375,0.5625)+c(1,0))
  expect_equal( facExp[[2]], matrix(c(0,0,0.4375,0.5625),2,2) + matrix(c(0.75,0,0.25,0),2,2) )
})

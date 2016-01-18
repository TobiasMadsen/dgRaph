library(dgRaph)
context("Data input check")

test_that("Wrong number of columns",{
  source("cases/twoDepedentVariables.R")
  mydfg <- twoDependentVariables()
  data <- matrix(c(NA), 1, 1)
  
  expect_error( .checkInputData(mydfg, data))
})

test_that("Data out of range",{
  source("cases/twoDepedentVariables.R")
  mydfg <- twoDependentVariables()
  data <- matrix(c(3,2), 1, 2)
  
  expect_error( .checkInputData(mydfg, data))
})

test_that("Data out of range 2",{
  varDim <- c(3,2)
  facPot <- list(multinomialPotential(c(1,3)), multinomialPotential(c(1,2)))
  facNbs <- list(1,2)
  mydfg <- dfg(varDim, facPot, facNbs)
  data <- matrix(c(3,3,NA,NA),2,2)
  
  expect_true({.checkInputData(mydfg, data); TRUE})
})

test_that("DataList",{
  source("cases/twoDepedentVariables.R")
  mydfg <- twoDependentVariables()
  data <- matrix(c(3,2), 1, 2)
  
  expect_error( .checkInputData(mydfg, data))
})
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

test_that("DataList",{
  source("cases/twoDepedentVariables.R")
  mydfg <- twoDependentVariables()
  data <- matrix(c(3,2), 1, 2)
  
  expect_error( .checkInputData(mydfg, data))
})
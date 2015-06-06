library(dgRaph)
context("Likelihood Calculation")

test_that("Likelihood Partial Observation", {
  source("cases/twoDepedentVariables.R")
  mydfg <- twoDependentVariables()
  
  data <- matrix(c(NA, 1), 1, 2)
  dataList <- list()
  dataList[[1]] <- list(c(0.8,0.2))
  
  expect_equal( likelihood(data = data, dfg = mydfg, dataList = dataList),
                0.7*0.8*0.75+0.3*0.2*0.25)
})

test_that("Likelihood Calculation 1",{
  source("cases/twoDepedentVariables.R")
  mydfg <- twoDependentVariables()
  
  data <- data.frame(O1 = c(1, 2, NA, NA), O2 = c(NA, 1, NA, 2))
  
  expect_equal(likelihood(data = data, dfg = mydfg), c(0.7, 0.3*0.25, 1, 0.7*0.25+0.3*0.75))
  expect_equal(likelihood(data = data, dfg = mydfg, log = T), log(c(0.7, 0.3*0.25, 1, 0.7*0.25+0.3*0.75)))
})

test_that("Likelihood Calculation 2",{
  source("cases/fourIndependentVariables.R")
  mydfg <- fourIndependentVariables()
  
  data <- data.frame(O1 = c(1,NA), O2 = c(NA,3), O3 = c(3,4), O4 = c(NA,NA))  
  expect_equal(likelihood(data = data, dfg = mydfg), c(0.05*1*0.05*1, 1*0.7*0.05*1))
  expect_equal(likelihood(data = data, dfg = mydfg, log = T), log(c(0.05*1*0.05*1, 1*0.7*0.05*1)))
})
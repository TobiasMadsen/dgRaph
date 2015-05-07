library(dgRaph)
context("Max probable state")

test_that("Dependent variables R level",{
  source("cases/twoDepedentVariables.R")
  mydfg <- twoDependentVariables()
  
  df <- data.frame(x = c(2,1,NA,NA,1), y = c(NA,NA,NA,2,2))
  df_mps <- data.frame(x = c(2,1,1,2,1), y = c(2,1,1,2,2))
  
  df_calc_mps <- mps(df, mydfg)
  expect_equal(df_calc_mps, df_mps)
})

test_that("Independent variables R level",{
  source("cases/fourIndependentVariables.R")
  mydfg <- fourIndependentVariables()
  
  df <- data.frame(x = 1:3, y = c(NA, 2, NA), z = 3:1, w = rep(2,3))
  df_mps <- data.frame(x = 1:3, y = c(3, 2, 3), z = 3:1, w = rep(2,3))
  df_calc_mps <- mps(df, mydfg)
  expect_equal(df_calc_mps, df_mps)
})

test_that("Independent variables C++ level",{
  source("cases/fourIndependentVariables.R")
  mydfg <- fourIndependentVariables()
  
  mpstates <- mydfg$dfgmodule$maxProbState(integer(0), logical(0))
  expect_equal( mpstates, c(4,3,1,2)) 
})

test_that("Independent variables cpp level with observed variables",{
  source("cases/fourIndependentVariables.R")
  mydfg <- fourIndependentVariables()
  
  mpstates <- mydfg$dfgmodule$maxProbState(c(1,1,2,1), c(F,F,T,F))
  expect_equal( mpstates, c(4,3,2,2) )
})

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

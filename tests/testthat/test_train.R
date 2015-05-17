library(dgRaph)
context("Training")

test_that("Normal Mixture",{
  expect_true(T)
})

test_that("Common Factors",{
  expect_true(T)
  
})

test_that("Markov Chain",{
  N <- 1000
  varDim <- rep(2, N)
  facPot <- c(list(matrix(0.5, 1, 2)),
              list(matrix(c(0.6,0.4,0.7,0.3), 2, 2)))
  potMap <- c(1, rep(2, N-1))
  facNbs <- c(list(c(1L)), lapply(2:N, FUN=function(i){c(i-1, i)}))
  mydfg <- dfg(varDim, facPot, facNbs, potMap = potMap)
  
  # Test
  df <- data.frame(matrix(c(1,NA,1,2),1,1000))
  sink("/dev/null") # Portability?
  train(df, mydfg)
  sink()
  facPotTrained <- potentials(mydfg)
  expect_equal(facPotTrained[[2]][1,1], 0, tolerance = 1e-3)
  expect_equal(facPotTrained[[2]][1,2], 1, tolerance = 1e-3)
  expect_equal(facPotTrained[[2]][2,1], 1, tolerance = 1e-3)
  expect_equal(facPotTrained[[2]][2,2], 0, tolerance = 1e-3)
})

test_that("IID Variables",{
  # IID variables
  source("cases/IIDVariables.R")
  mydfg <- IIDVariables(10)
  
  # Test
  df <- data.frame(matrix(c(2,1,1,2,2,2,2,NA,2,1),1,10))
  train(df, mydfg)
  facPotTrained <- potentials(mydfg)
  expect_equal(facPotTrained[[1]][1,1], 0.333333, tolerance = 1e-3)
  expect_equal(facPotTrained[[1]][1,2], 0.666667, tolerance = 1e-3)
})
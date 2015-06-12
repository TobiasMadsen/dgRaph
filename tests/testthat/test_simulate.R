library(dgRaph)
context("Sampling")

test_that("Set seed",{
  source("cases/twoDepedentVariables.R")
  mydfg <- twoDependentVariables()
  
  set.seed(1)
  sim1 <- simulate(mydfg, 50)
  
  set.seed(1)
  sim2 <- simulate(mydfg, 50)
  
  expect_equal(sim1, sim2)
})

test_that("Sampling from 2 state graph",{
  source("cases/twoDepedentVariables.R")
  mydfg <- twoDependentVariables()
  
  sim <- simulate(mydfg, 1000)
  expect_equal( nrow(sim), 1000)
  expect_equal( ncol(sim), 2)
  expect_equal( colnames(sim), c('x', 'y'))
  
  # Stochastic tests will fail with probability ~1e-4
  expect_less_than( sum(sim$x == 1) , 753 )
  expect_more_than( sum(sim$x == 1) , 644)
})

test_that("Sampling unnormalized",{
  varDim <- c(2)
  facPot <- list(matrix(c(4,1),1,2))
  facNbs <- list(1)
  mydfg <- dfg(varDim, facPot, facNbs)
  sim <- simulate(mydfg, 1000)
  
  # Stochastic test
  expect_less_than( sum(sim[,1] == 1), 900)
  expect_more_than( sum(sim[,1] == '1'), 700)
})
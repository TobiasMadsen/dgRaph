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
  expect_lt( sum(sim$x == 1) , 753 )
  expect_gt( sum(sim$x == 1) , 644)
})

test_that("Sampling unnormalized 1",{
  varDim <- c(2)
  facPot <- list(matrix(c(4,1),1,2))
  facNbs <- list(1)
  mydfg <- dfg(varDim, facPot, facNbs)
  sim <- simulate(mydfg, 1000)
  
  # Stochastic test
  expect_lt( sum(sim[,1] == 1), 900)
  expect_gt( sum(sim[,1] == '1'), 700)
})

test_that("Sampling unnormalized 2",{
  varDim <- rep(2,2)
  facPot <- c(list(matrix(c(0.75,0.25),1,2)),
              list(matrix(c(0.75,0.25,0.25,0.75),2,2)))
  facNbs <- c(list(c(1L)),
              list(c(1L,2L)))
  mydfg <- dfg(varDim, facPot, facNbs)  
  
  sim <- simulate(mydfg, 100000)
  
  expect_true( findInterval(sum(sim[,1] == 1 & sim[,2] == 1), 9/16*100000*c(0.975,1.025)) == 1 )
  expect_true( findInterval(sum(sim[,1] == 2 & sim[,2] == 1), 1/16*100000*c(0.975,1.025)) == 1 )
  expect_true( findInterval(sum(sim[,1] == 1 & sim[,2] == 2), 3/16*100000*c(0.975,1.025)) == 1 )
  expect_true( findInterval(sum(sim[,1] == 2 & sim[,2] == 2), 3/16*100000*c(0.975,1.025)) == 1 )
})
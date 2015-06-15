library(dgRaph)
context("Tail approximations")

test_that("IS sampling two dependent variables",{
  varDim <- rep(2,2)
  facPot <- c(list(matrix(0.5,1,2)),
              list(matrix(0.5,2,2)))
  facPotFg <- c(list(matrix(c(0.75,0.25),1,2)),
                list(matrix(c(0.75,0.25,0.25,0.75),2,2)))
  facNbs <- c(list(c(1L)),
              list(c(1L,2L)))
  
  mydfg <- dfg(varDim, facPot, facNbs)
  is <- tailIS(x = 0.81, n = 10000, alpha = 1.5, dfg = mydfg, facPotFg = facPotFg)
  naive <- tailIS(x = 0.81, n = 10000, alpha = 0.0, dfg = mydfg, facPotFg = facPotFg)
  
  # Test estimates are close
  expect_less_than( abs(is$p - naive$p), 0.02)
  # Test confidence intervals overlap
  expect_less_than( is$low, naive$high)
  expect_less_than( naive$low, is$high)
})

test_that("NA when out of range saddlepoint",{
  varDim <- 2L
  facPot <- c(list(matrix(c(0.5,0.5),1,2)) )
  facNbs <- c(list(c(1L)) )
  
  mydfg <- dfg(varDim, facPot, facNbs)
  
  #Foreground model
  facPotFg <- c(list(matrix(c(0.9,0.1),1,2)))
  
  #Score range [log(0.1/0.5);log(0.9/0.5)] = [-1.609, .588]
  dfsaddle <- tailSaddle( seq(-2,1,0.1), mydfg, facPotFg)
  expect_equal( which(is.na(dfsaddle$p)), c(1:4,27:31) )
})

test_that("IS sampling binomial",{
  varDim <- 2
  facPot <- list(matrix(c(0.01, 0.99),1 ,2))
  facNbs <- list(1)
  mydfg <- dfg(varDim, facPot, facNbs)
  
  facPotFg <- list(matrix(c(0.2, 0.8), 1, 2))
  
  tail_df <- tailIS(x = 2.98, n = 1000, alpha = 3, dfg = mydfg, facPotFg = facPotFg)
  
  # Stochastic test
  expect_less_than(tail_df$p[1], 0.011)
  expect_more_than(tail_df$p[1], 0.009)
})

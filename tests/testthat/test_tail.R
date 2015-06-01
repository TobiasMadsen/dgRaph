library(dgRaph)
context("Tail approximations")

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

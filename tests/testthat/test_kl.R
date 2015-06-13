library(dgRaph)
context("Kullback-Leibler")

test_that("Small probabilities",{
  varDim <- c(2)
  facPotBg <- list(matrix(c(1, 1e-14),1,2))
  facPotFg <- list(matrix(c(1,0),1,2))
  facNbs <- list(1)

  dfg1 <- dfg(varDim, facPotBg, facNbs)
  dfg2 <- dfg(varDim, facPotFg, facNbs)
  
  expect_false( is.infinite( kl(dfg1, dfg2)) )
})

test_that("Two independent variables",{
  varDim <- c(2,2)
  facPotBg <- list(matrix(c(0.5),1,2))
  facPotFg <- list(matrix(c(0.7,0.3),1,2))
  facNbs <- list(1,2)
  potMap <- c(1,1)
  dfg1 <- dfg(varDim, facPotBg, facNbs, potMap)
  dfg2 <- dfg(varDim, facPotFg, facNbs, potMap)
  
  expect_equal( unname(kl(dfg1, dfg2)) , ((log(5/7)+log(5/3))*0.5)*2)
})

test_that("Two dependent variables",{
  varDim <- rep(2,2)
  facPotBg <- list(matrix(c(0.5),1,2),
                   matrix(c(0.5),2,2))
  facPotFg <- list(matrix(c(0.75,0.25),1,2),
                   matrix(c(0.75,0.25,0.25,0.75),2,2))
  facNbs <- list(c(1L),
                 c(1L,2L))
  
  dfg1 <- dfg(varDim, facPotBg, facNbs, varNames = c('x', 'y'))
  dfg2 <- dfg(varDim, facPotFg, facNbs, varNames = c('x', 'y'))  
  
  expect_equal( unname(kl(dfg1, dfg2)), (log(4/9)+log(4/3)+log(4/3)+log(4))/4 )
})

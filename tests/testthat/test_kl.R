library(dgRaph)
context("Kullback-Leibler")

test_that("Remapping of structures 1",{
  varDim <- c(2,2)
  facPot1 <- list(matrix(0.5,1,2))
  facPot2 <- list(matrix(c(0.2,0.8),1,2),
                  matrix(c(0.8,0.2),1,2))
  facNbs <- list(1,2)
  potMap1 <- c(1,1)
  potMap2 <- c(1,2)
  
  dfg1 <- dfg(varDim, facPot1, facNbs, potMap1)
  dfg2 <- dfg(varDim, facPot2, facNbs, potMap2)
  
  expect_equal( unname(kl(dfg1, dfg2)), log(0.5/0.8)+log(0.5/0.2))
})

test_that("Remapping of structures 2",{
  varDim <- c(2,2,2)
  facPot <- list(matrix(c(0.2,0.8),1,2),
                 matrix(c(0.3,0.7),1,2),
                 matrix(c(0.8,0.2),1,2))
  facNbs1 <- list(1,2,3)
  facNbs2 <- list(2,3,1)
  
  dfg1 <- dfg(varDim, facPot, facNbs1)
  dfg2 <- dfg(varDim, facPot, facNbs2)
  
  expect_equal( unname(kl(dfg1, dfg2)), 
                (0.2*log(0.2/0.8)+0.8*log(0.8/0.2))+
                  (0.2*log(0.2/0.7)+0.8*log(0.8/0.3))+
                  (0.3*log(0.3/0.2)+0.7*log(0.7/0.8)))
})

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

library(PGMscore)
context("Expect algorithm")

test_that("Simple phylogenetic model",{
  varDim <- rep(2L,4)
  facPot <- c(list(matrix(c(0.8,0.2),1,2)),
              rep(list(matrix(c(0.8,0.2,0.2,0.8),2,2)),3))
  facExp <- c(list(matrix(0,1,2)),
              rep(list(matrix(c(0,1,1,0),2,2)),3))
  facNbs <- c(list(c(0L)),
              list(c(0L,1L)),
              list(c(0L,2L)),
              list(c(2L,3L)))
  res <- mgfDFG2(varDim,facPot, facExp, facNbs)
  expect_equal(length(res),2)
  expect_equal(res[1],1)
  expect_equal(res[2],0.6)
})

test_that("Components",{
  varDim <- rep(1L,2)
  facPot <- c(list(matrix(3,1,1)),
              list(matrix(2,1,1)))
  facExp <- c(list(matrix(1,1,1)),
              list(matrix(2,1,1)))
  facNbs <- c(list(c(0L)),
              list(c(1L)))
  res <- mgfDFG2(varDim,facPot, facExp, facNbs)
  expect_equal(res[1],6)
  expect_equal(res[2],18)
})

test_that("Bernoulli connected to ancestor",{
  N <- 8
  t <- 0.25
  varDim <- rep(2,N+1)
  facPot1 <- c( list(matrix(c(0.5,0.5),1,2)), rep(list(matrix(c(0.50,0.50,0.50,0.50),2,2)),N) )
  facPot2 <- c( list(matrix(c(0.5,0.5),1,2)), rep(list(matrix(c(0.75,0.25,0.25,0.75),2,2)),N) )
  facNbs <- c(list(c(0L)), lapply(1:N, function(x){c(0,x)}))
  res <- mgfDFG2(varDim, facPotToFunA(facPot1,facPot2,t), facPotToFunB(facPot1,facPot2), facNbs)
  mgf <- function(t){(((1/2)**t+(3/2)**t)/2)**N}
  dtlogmgf <- function(t){N/((1/2)**t+(3/2)**t)*(log(1/2)*(1/2)**t+log(3/2)*(3/2)**t )}
  
  expect_equal(mgf(t),res[1])
  expect_equal(dtlogmgf(t),res[2]/res[1])
})

test_that("Bernoulli disconnected",{
  N <- 4
  varDim <- rep(2L,N)
  facPot1 <- rep(list(matrix(c(0.5,0.5),1,2)),N)
  facPot2 <- rep(list(matrix(c(0.75,0.25),1,2)),N)
  facNbs <- lapply(0:(N-1),function(i){c(i)})
  t <- 0.25
  res <- mgfDFG2(varDim, facPotToFunA(facPot1, facPot2,t), facPotToFunB(facPot1, facPot2), facNbs)
  mgf <- function(t){(((1/2)**t+(3/2)**t)/2)**N}
  dtlogmgf <- function(t){N/((1/2)**t+(3/2)**t)*(log(1/2)*(1/2)**t+log(3/2)*(3/2)**t )}
  expect_equal(mgf(t),res[1])
  expect_equal(dtlogmgf(t), res[2]/res[1])
})
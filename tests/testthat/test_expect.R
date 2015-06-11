library(dgRaph)
context("Expect")

test_that("Components", {
  varDim <- rep(1L,2)
  facPot <- c(list(matrix(3,1,1)),
              list(matrix(2,1,1)))
  facExp <- c(list(matrix(1,1,1)),
              list(matrix(2,1,1)))
  facNbs <- c(list(c(1L)),
              list(c(2L)))
  
  mydfg <- dfg(varDim, facPot, facNbs)
  
  res <- unname( expect(mydfg, facExp) )
  expect_equal(unname(res[1]),6)
  expect_equal(unname(res[2]),18)
})

test_that("Simple phylogenetic model",{
  varDim <- rep(2L,4)
  facPot <- c(list(matrix(c(0.8,0.2),1,2)),
              rep(list(matrix(c(0.8,0.2,0.2,0.8),2,2)),3))
  facExp <- c(list(matrix(0,1,2)),
              rep(list(matrix(c(0,1,1,0),2,2)),3))
  facNbs <- c(list(c(1L)),
              list(c(1L,2L)),
              list(c(1L,3L)),
              list(c(3L,4L)))
  mydfg <- dfg(varDim, facPot, facNbs)
  
  res <- unname( expect(mydfg, facExp) )
  expect_equal(length(res),2)
  expect_equal(res[1],1)
  expect_equal(res[2],0.6)
})

test_that("Bernoulli connected to ancestor",{
  N <- 15
  t <- 0.25
  varDim <- rep(2,N+1)
  facPot1 <- c( list(matrix(c(0.5,0.5),1,2)), rep(list(matrix(c(0.50,0.50,0.50,0.50),2,2)),N) )
  facPot2 <- c( list(matrix(c(0.5,0.5),1,2)), rep(list(matrix(c(0.75,0.25,0.25,0.75),2,2)),N) )
  facNbs <- c(list(c(1L)), lapply(1:N, function(x){c(1,x+1)}))
  
  facPot <- .facPotToFunA(facPot1,facPot2,t)
  facExp <- .facPotToFunB(facPot1,facPot2)
  
  mydfg <- dfg(varDim, facPot, facNbs)
  res <- unname(expect(mydfg, facExp))
  
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
  facNbs <- lapply(1:N,function(i){c(i)})
  t <- 0.25
  
  facPot <- .facPotToFunA(facPot1,facPot2,t)
  facExp <- .facPotToFunB(facPot1,facPot2)
  
  mydfg <- dfg(varDim, facPot, facNbs)
  res <- unname(expect(mydfg, facExp))
  
  mgf <- function(t){(((1/2)**t+(3/2)**t)/2)**N}
  dtlogmgf <- function(t){N/((1/2)**t+(3/2)**t)*(log(1/2)*(1/2)**t+log(3/2)*(3/2)**t )}
  expect_equal(mgf(t),res[1])
  expect_equal(dtlogmgf(t), res[2]/res[1])
})

library(dgRaph)
context("Underflow")

test_that("Variable 10 neighbors",{
  N <- 10
  source("cases/nNeighbors.R")
  mydfg <- nNeighbors(N)
  
  # Calculate statistics
  data <- data.frame(matrix(rep(1,N+1),1,N+1))
  expect_equal(likelihood(data, dfg = mydfg,log=T), N*log(0.01))
})

test_that("Variable 100 neighbors",{
  N <- 100
  source("cases/nNeighbors.R")
  mydfg <- nNeighbors(N)
  
  # Calculate statistics
  data <- data.frame(matrix(rep(1,N+1),1,N+1))
  expect_equal(likelihood(data, dfg = mydfg,log=T), N*log(0.01))
})

test_that("Variable 1000 neighbors",{
  N <- 1000
  source("cases/nNeighbors.R")
  mydfg <- nNeighbors(N)
  
  # Calculate statistics
  data <- data.frame(matrix(rep(1,N+1),1,N+1))
  expect_equal(likelihood(data, dfg = mydfg,log=T), N*log(0.01))
})

test_that("Graph 10 Nodes",{
  N <- 10
  source("cases/nDependentVariables.R")
  mydfg <- nDependentVariables(N)

  # Calculate statistics
  data <- data.frame(matrix(rep(1,N),1,N))
  expect_equal( likelihood(data = data, dfg = mydfg), 0.5*0.75**(N-1) )
})

test_that("Graph 100 Nodes",{
  N <- 100
  source("cases/nDependentVariables.R")
  mydfg <- nDependentVariables(N)
  
  # Calculate statistics
  data <- data.frame(matrix(rep(1,N),1,N))
  expect_equal( likelihood(data = data, dfg = mydfg), 0.5*0.75**(N-1) )
})

test_that("Graph 1000 Nodes",{
  N <- 1000
  source("cases/nDependentVariables.R")
  mydfg <- nDependentVariables(N)
  
  # Calculate statistics
  data <- data.frame(matrix(rep(1,N),1,N))
  expect_equal( likelihood(data = data, dfg = mydfg), 0.5*0.75**(N-1) )
})

test_that("Graph 10000 Nodes",{
  N <- 10000
  source("cases/nDependentVariables.R")
  mydfg <- nDependentVariables(N)
  
  # Calculate statistics
  data <- data.frame(matrix(rep(1,N),1,N))
  expect_equal( likelihood(data = data, dfg = mydfg, log = T), log(0.5)+(N-1)*log(0.75) )
})
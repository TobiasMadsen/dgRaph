library(dgRaph)
context("Training")

test_that("Training returns dfg",{
  mydfg <- IIDVariables(1)
  data <- data.frame(matrix(c(2,1,1,2,2,2,2,NA,2,1),10,1))
  
  optDfg <- train(data, mydfg)

  expect_true(is.dfg(optDfg))
  expect_false( identical(optDfg, mydfg))
})

test_that("Custom Optimization", {
  mypotential <- list(mat = matrix(c(0.2,0.2,0.3,0.3), 1, 4),
                      param = NULL)
  class(mypotential) <- c("uniform","potential")
  update.uniform <- function(pot, expCounts){
    ret <- list(mat = matrix(0.25,1,4))
    class(ret) <- c("uniform", "potential")
    ret
  }
  
  varDim <- c(4)
  facPot <- list(mypotential)
  facNbs <- list(1)
  mydfg  <- dfg(varDim, facPot, facNbs)
  
  # Custom training function
  df <- data.frame(O1 = c(1,2,3))
  
  optDfg <- train(df, mydfg)
  
  # Test
  expect_equal( as.matrix(potentials(optDfg)[[1]]), matrix(0.25,1,4))
})

test_that("Beta Distribution", {
  varDim <- c(100)
  set.seed(1)
  facPot <- c(list(betaPotential(c(1,100))))
  facNbs <- c(list(1))
  mydfg <- dfg(varDim, facPot, facNbs)
  
  # Generate data
  set.seed(1)
  df <- data.frame(O1 = ceiling(rbeta(500, 4,6)*100))
  
  # Train 
  optDfg <- train(df, mydfg)
    
  # Compare Kolmogorov Smirnoff style
  expect_less_than(
    max(abs(pbeta(1:100/100, 4, 6) - cumsum(as.matrix(potentials(optDfg)[[1]])[1,]))),
    0.02)
  
})

test_that("Normal Mixture",{
  varDim <- c(2,100)
  set.seed(1)
  facPot <- c(list(multinomialPotential(c(1,2))),
              list(normalPotential(c(2, 100))) )
  facNbs <- c(list(c(1)),
              list(c(1,2)))
  mydfg <- dfg(varDim, facPot, facNbs)
  
  # Generate data
  set.seed(1)
  r1 <- round(rnorm(400, 40, 8))
  r2 <- round(rnorm(800, 70, 6))
  df <- data.frame(H1 = NA, O1 = sample(c(r1,r2)))
  
  optDfg <- train(df, mydfg, threshold = 1e-5)

  # Check that decision boundary is in the 50-60's
  newFacPot <- potentials(optDfg)
  classification <- apply( sweep(as.matrix(newFacPot[[2]]), 1, STATS = as.matrix(newFacPot[[1]]), FUN = "*"), 2, which.max)
  expect_equal( length(unique(classification[1:50])), 1)
  expect_equal( length(unique(classification[65:100])), 1)
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
  optDfg <- train(df, mydfg)
    
  facPotTrained <- potentials(optDfg)
  expect_equal(facPotTrained[[2]][1,1], 0, tolerance = 1e-5)
  expect_equal(facPotTrained[[2]][1,2], 1, tolerance = 1e-5)
  expect_equal(facPotTrained[[2]][2,1], 1, tolerance = 1e-5)
  expect_equal(facPotTrained[[2]][2,2], 0, tolerance = 1e-5)
})

test_that("IID Variables",{
  # IID variables
  source("cases/IIDVariables.R")
  mydfg <- IIDVariables(10)
  
  # Test
  df <- data.frame(matrix(c(2,1,1,2,2,2,2,NA,2,1),1,10))
  
  optDfg <- train(df, mydfg)
  
  facPotTrained <- potentials(optDfg)
  expect_equal(facPotTrained[[1]][1,1], 0.333333, tolerance = 1e-5)
  expect_equal(facPotTrained[[1]][1,2], 0.666667, tolerance = 1e-5)
})
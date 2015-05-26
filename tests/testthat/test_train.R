library(dgRaph)
context("Training")

test_that("Custom Optimization", {
  varDim <- c(4)
  facPot <- list(matrix(c(0.2,0.2,0.3,0.3), 1, 4))
  facNbs <- list(1)
  mydfg  <- dfg(varDim, facPot, facNbs)
  
  # Custom training function
  optimList <- list("uniform" = function(expCounts){return(list(pot = matrix(0.25,1,4), str = "Uniform\n"))})
  df <- data.frame(O1 = c(1,2,3))
  
  tryCatch({
    sink("/dev/null")
    train(df, mydfg, optim = 'uniform', optimFun = optimList)
  },
  finally = {sink()})
  
  # Test
  expect_equal( potentials(mydfg)[[1]], matrix(0.25,1,4))
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
  tryCatch({
    sink("/dev/null")
    train(df, mydfg, optim = 'beta')
  },
  finally = {sink()})
  
  # Compare Kolmogorov Smirnoff style
  expect_less_than(
    max(abs(pbeta(1:100/100, 4, 6) - cumsum(potentials(mydfg)[[1]][1,]))),
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
  
  tryCatch({
    sink("/dev/null")
    train(df, mydfg, optim = c('row', 'norm'), threshold = 1e-5)
  }, 
  finally={sink()})

  # Check that decision boundary is in the 50-60's
  newFacPot <- potentials(mydfg)
  classification <- apply( sweep(newFacPot[[2]], 1, STATS = newFacPot[[1]], FUN = "*"), 2, which.max)
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
  tryCatch({
    sink("/dev/null") # Portability?
    train(df, mydfg)
  },
  finally = {sink()})
  
  facPotTrained <- potentials(mydfg)
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
  tryCatch({
    sink("/dev/null")
    train(df, mydfg)
  },
  finally = {sink()})
  facPotTrained <- potentials(mydfg)
  expect_equal(facPotTrained[[1]][1,1], 0.333333, tolerance = 1e-5)
  expect_equal(facPotTrained[[1]][1,2], 0.666667, tolerance = 1e-5)
})
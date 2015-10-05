library(dgRaph)
context("Tail approximations")

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
  
  # Most extreme score log(0.8/0.5)*2 ~= 0.94
  # With probability 0.5*0.5 = 0.25
  naive <- tailIS(x = 0.93, n = 4000, alpha = 0, dfg1 = dfg1, dfg2 = dfg2)
  is    <- tailIS(x = 0.93, n = 4000, alpha = 1.5, dfg1 = dfg1, dfg2 = dfg2)
  expect_less_than( abs(naive$p-0.25), 0.02)
  expect_less_than( abs(is$p-0.25), 0.02)
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
  
  # Most extreme score log(0.8/0.2)+log(0.8/0.7)+log(0.7/0.2) ~= 2.773
  # With probability 0.2*0.7*0.2 = 0.028
  naive <- tailIS(x = 2.7, n = 4000, alpha = 0, dfg1 = dfg1, dfg2 = dfg2)
  is    <- tailIS(x = 2.7, n = 4000, alpha = 1.5, dfg1 = dfg1, dfg2 = dfg2)
  expect_less_than( abs(naive$p-0.028), 0.01)
  expect_less_than( abs(is$p-0.028), 0.001) 
  
})

test_that("IS sampling two dependent variables",{
  varDim <- rep(2,2)
  facPot <- c(list(matrix(0.5,1,2)),
              list(matrix(0.5,2,2)))
  facPotFg <- c(list(matrix(c(0.75,0.25),1,2)),
                list(matrix(c(0.75,0.25,0.25,0.75),2,2)))
  facNbs <- c(list(c(1L)),
              list(c(1L,2L)))
  
  mydfg <- dfg(varDim, facPot, facNbs)
  mydfgFg <- dfg(varDim, facPotFg, facNbs)
  is <- tailIS(x = 0.81, n = 10000, alpha = 1.5, dfg1 = mydfg, dfg2= mydfgFg)
  naive <- tailIS(x = 0.81, n = 10000, alpha = 0.0, dfg = mydfg, dfg2 = mydfgFg)
  
  # Test estimates are close
  expect_less_than( abs(is$p - naive$p), 0.02)
  # Test confidence intervals overlap
  expect_less_than( is$low, naive$high)
  expect_less_than( naive$low, is$high)
})

test_that("NA when out of range saddlepoint",{
  varDim <- 2L
  facPot <- c(list(matrix(c(0.5,0.5),1,2)) )
  facPotFg <- c(list(matrix(c(0.9,0.1),1,2)))
  facNbs <- c(list(c(1L)) )
  
  dfg1 <- dfg(varDim, facPot, facNbs)
  dfg2 <- dfg(varDim, facPotFg, facNbs)
  
  #Score range [log(0.1/0.5);log(0.9/0.5)] = [-1.609, .588]
  dfsaddle <- tailSaddle( seq(-2,1,0.1), dfg1 = dfg1, dfg2 = dfg2)
  expect_equal( which(is.na(dfsaddle$p)), c(1:4,27:31) )
})

test_that("Saddlepoint single variable",{
  # Binomial distributed score B(1/2,4)
  varDim <- 5
  facPot <- list(matrix(dbinom(0:4,4,1/2),1,5))
  facPotFg <- list(exp(matrix(0:4,1,5))*facPot[[1]])
  facNbs <- list(1)
  mydfg <- dfg(varDim, facPot, facNbs)
  mydfgFg <- dfg(varDim, facPotFg, facNbs)

  tail_df <- tailSaddle(x = c(2,3), dfg1 = mydfg, dfg2 = mydfgFg)

  # Approximation at mean
  expect_equal(tail_df$p[1], 0.5)

  # Saddlepoint in 3
  n <- 4
  theta <- log(3)
  v <- 3/4
  sp <- ((1+exp(theta))/2)**n*exp(-3*theta)*exp((theta**2*v)/2)*(1-pnorm(theta*sqrt(v)))

  expect_equal( tail_df$p[2], sp, tolerance = 1e-5)
})

test_that("Normal single variable",{
  # Binomial distributed score B(1/2,4)
  varDim <- 5
  facPot <- list(matrix(dbinom(0:4,4,1/2),1,5))
  facPotFg <- list(exp(matrix(0:4,1,5))*facPot[[1]])
  facNbs <- list(1)
  mydfg <- dfg(varDim, facPot, facNbs)
  mydfgFg <- dfg(varDim, facPotFg, facNbs)
  
  tail_df <- tailNormal(x = c(2,3), dfg1 = mydfg, dfg2 = mydfgFg)
  
  # Approximation at mean
  expect_equal(tail_df$p[1], 0.5)
  
  # Normal approx in 3
  np <- pnorm(3, mean = 2, sd = 1, lower.tail = F)
  
  expect_equal( tail_df$p[2], np)
})

test_that("Saddlepoint multiple variables",{
  # Binomial distributed score B(1/2,8)
  varDim <- c(5,5)
  facPot <- list(matrix(dbinom(0:4,4,1/2),1,5),
                 matrix(dbinom(0:4,4,1/2),5,5,byrow =T))
  facPotFg <- list(exp(matrix(0:4,1,5))*facPot[[1]],
                   exp(matrix(0:4,5,5,byrow=T))*facPot[[2]])
  facNbs <- list(1,c(1,2))
  mydfg <- dfg(varDim, facPot, facNbs)
  mydfgFg <- dfg(varDim, facPotFg, facNbs)

  tail_df <- tailSaddle(x = c(4,6), dfg1 = mydfg, dfg2 = mydfgFg)
  
  # Approximation at mean
  test_that(tail_df$p[1], 0.5)

  # Saddlepoint in 6
  n <- 8
  theta <- log(3)
  v <- 3/2
  sp <- ((1+exp(theta))/2)**n*exp(-6*theta)*exp((theta**2*v)/2)*(1-pnorm(theta*sqrt(v)))

  expect_equal( tail_df$p[2], sp, tolerance = 1e-5)
})

test_that("Normal multiple variables",{
  # Binomial distributed score B(1/2,8)
  varDim <- c(5,5)
  facPot <- list(matrix(dbinom(0:4,4,1/2),1,5),
                 matrix(dbinom(0:4,4,1/2),5,5,byrow =T))
  facPotFg <- list(exp(matrix(0:4,1,5))*facPot[[1]],
                   exp(matrix(0:4,5,5,byrow=T))*facPot[[2]])
  facNbs <- list(1,c(1,2))
  mydfg <- dfg(varDim, facPot, facNbs)
  mydfgFg <- dfg(varDim, facPotFg, facNbs)
  
  tail_df <- tailNormal(x = c(4,6), dfg1 = mydfg, dfg2 = mydfgFg)
  
  # Approximation at mean
  test_that(tail_df$p[1], 0.5)
  
  # Normal Approx in 6
  np <- pnorm(6, mean = 4, sd = sqrt(2), lower.tail = F)
  
  expect_equal( tail_df$p[2], np, tolerance = 1e-5)
})

test_that("Lattice correction saddlepoint",{
  # Binomial distributed score B(1/2,4)
  varDim <- 5
  facPot <- list(matrix(dbinom(0:4,4,1/2),1,5))
  facPotFg <- list(exp(matrix(0:4,1,5))*facPot[[1]])
  facNbs <- list(1)
  mydfg <- dfg(varDim, facPot, facNbs)
  mydfgFg <- dfg(varDim, facPotFg, facNbs)

  tail_df <- tailSaddle(x = c(2,3), dfg1 = mydfg, dfg2 = mydfgFg)
  tail_df_lattice <- tailSaddle(x = c(2,3), dfg1 = mydfg, dfg2 = mydfgFg, lattice = 1)

  # No lattice correction at mean as saddlepoint is zero
  expect_equal(tail_df_lattice$p[1], tail_df$p[1])

  # Saddlepoint at 3 is 
  expect_equal(tail_df_lattice$p[2], tail_df$p[2]*log(3)*3/2)
})

test_that("IS sampling binomial",{
  varDim <- 2
  facPot <- list(matrix(c(0.01, 0.99),1 ,2))
  facPotFg <- list(matrix(c(0.2, 0.8), 1, 2))
  facNbs <- list(1)
  mydfg <- dfg(varDim, facPot, facNbs)
  mydfgFg <- dfg(varDim, facPotFg, facNbs)
  
  tail_df <- tailIS(x = 2.98, n = 1000, alpha = 3, dfg1 = mydfg, dfg2 = mydfgFg)
  
  # Stochastic test
  expect_less_than(tail_df$p[1], 0.011)
  expect_more_than(tail_df$p[1], 0.009)
})

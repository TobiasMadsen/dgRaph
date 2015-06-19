library(dgRaph)
context("Optimizing")

test_that("Beta optimizer",{
  # Generate data
  set.seed(1)
  dat <- rbeta(1000, 4, 140)
  expCounts <- matrix( hist(dat, breaks = seq(0, 0.1, 0.001), plot = FALSE)$counts, 1, 100)

  # Optimize
  betaOpt <- betaOptimize(range = c(0,0.1))(expCounts)
  pot <- betaOpt$pot
  alpha <- betaOpt$alphas[1]
  beta  <- betaOpt$betas[1]
  
  # Tests
  expect_equal( alpha, 4, 0.005)
  expect_equal( beta, 140, 0.005)
  
  # Reflection
  alpha <- betaOpt$alphas[1]
  beta  <- betaOpt$betas[1]
  expect_equal( pot, 
                betaPotential(dim = c(1,100), 
                              range = c(0,0.1),
                              alphas = alpha, 
                              betas = beta))
})

test_that("Linreg optimizer 1",{
  # Generate data
  dat <- matrix(c(0,0,1,0,0,
                  0,1,2,0,0,
                  1,2,1,0,0,
                  2,1,0,0,0,
                  1,0,0,0,0), 5,5)
  
  # Optimize
  linreg <- linregOptimize(range1 = c(0.5,5.5), range2 = c(0.5,5.5))(dat)
  
  # Tests
  expect_equal(linreg$alpha, -1)
  expect_equal(linreg$beta, 5)
  expect_equal(linreg$var, 0.5)
  
  # Reflection
  expect_equal(linreg$pot, 
               linregPotential(dim = dim(dat), 
                               range1 = c(0.5,5.5), 
                               range2 = c(0.5,5.5), 
                               alpha = -1, 
                               beta = 5, 
                               var = 0.5))
})

test_that("Linreg optimizer 2",{
  # Test scaling
  # Generate data
  dat <- matrix(c(0,0,1,0,
                  0,1,2,0,
                  1,2,1,0,
                  2,1,0,0,
                  1,0,0,0), 4,5)
  
  # Optimize
  # Bins x corresponds to 1,2,3,4
  # Bins y corresponds to 1,2,3,4,5
  linreg <- linregOptimize(range1 = c(0.5,4.5), range2 = c(0.5,5.5))(dat) 
  
  # Tests
  expect_equal(linreg$alpha, -1)
  expect_equal(linreg$beta, 5)
  expect_equal(linreg$var, 0.5)
  
  # Reflection
  expect_equal(linreg$pot, 
               linregPotential(dim = dim(dat), 
                               range1 = c(0.5,4.5), 
                               range2 = c(0.5,5.5), 
                               alpha = -1, 
                               beta = 5, 
                               var = 0.5))
  
  # Optimize
  # Bins x corresponds to 2,4,6,8
  # Bins y corresponds to 3,7,11,15,19
  linreg <- linregOptimize(range1 = c(1,9), range2 = c(1,21))(dat)
  
  # Tests
  expect_equal(linreg$alpha, -2)
  expect_equal(linreg$beta, 19)
  expect_equal(linreg$var, 8)
  
  # Reflection
  expect_equal(linreg$pot, 
               linregPotential(dim = dim(dat), 
                               range1 = c(1,9), 
                               range2 = c(1,21), 
                               alpha = -2, 
                               beta = 19, 
                               var = 8))
})

test_that("Linreg/Normal optimizer comparison 1",{
  # Generate data
  dat <- matrix(c(0,0,1,2,1,0),6,6, byrow = T)
  
  # Optimize
  linreg <- linregOptimize()(dat)
  
  # Rows identical
  expect_equal(linreg$pot[1,], linreg$pot[2,])
  
  # Check equal to similar one dimensional update
  datRow <- dat[1,,drop = F]
  normopt <- normOptimize()(datRow)
  expect_equal(normopt$pot[1,], linreg$pot[1,])
})

test_that("Normal optimizer 1",{
  # Generate data
  dat <- matrix(c(0,0,1,2,1,0,
                  0,0.5,4,0.5,0,0), 2, 6, byrow = T)
  
  # Optimize
  normopt <- normOptimize(range = c(1,13))(dat)
  
  # Tests
  expect_equal(normopt$means, c(8,6) )
  expect_equal(normopt$vars, c(2,0.8) )
  expect_equal(which.max(normopt$pot[1,]), 4)
  expect_equal(which.max(normopt$pot[2,]), 3)
  
  # Reflection
  expect_equal(normopt$pot, 
               normalPotential(dim = c(2,6), 
                               range = c(1,13), 
                               means = c(8,6), 
                               vars = c(2,0.8)))
})

test_that("Normal optimizer 2",{
  # Generate data
  dat <- matrix(rep(0,200), 1, 200)
  dat[1, 99] <- 2
  dat[1,100] <- 3
  dat[1,101] <- 3
  dat[1,102] <- 2
  
  # Optimize
  normopt <- normOptimize(range = c(100,300))(dat)
  
  # Tests
  expect_equal(normopt$means[1], 200)
  expect_equal(normopt$vars[1], 1.05) # Uses MLE for variance

  # Reflection
  expect_equal(normopt$pot, 
               normalPotential(dim = c(1,200), 
                               range = c(100,300), 
                               means = 200, 
                               vars = 1.05))
})

test_that("Fixedlink optimizer 1",{
  # Generate data
  dat <- matrix(c(1,1,0,0,0,
                  0,1,1,0,0,
                  0,0,1,1,0,
                  0,0,0,1,1), 4, 5, byrow = T)
  
  # Optimize
  fixedopt <- fixedlinkOptimize(range1 = c(0.5,4.5), range2 = c(0.5,5.5), alpha = 1, beta = 0)(dat)
  
  # Tests
  expect_equal(fixedopt$var, 0.5)
  
  # Reflection
  expect_equal(fixedopt$pot, 
               linregPotential(dim = dim(dat), 
                               range1 = c(0.5,4.5), 
                               range2 = c(0.5,5.5), 
                               alpha = 1, 
                               beta = 0, 
                               var = 0.5))
})

test_that("Fixedlink optimizer 2",{
  # Generate data
  dat <- matrix(c(0,1,1,0,0,
                  0,0,0,0,0,
                  0,0,1,1,0,
                  0,0,0,0,0), 4, 5, byrow = T)
  
  # Optimize
  fixedopt <- fixedlinkOptimize(range1 = c(2.5,6.5), range2 = c(0,10), alpha = 2, beta = -3)(dat)
  
  # Tests
  expect_equal(fixedopt$var, 2)
  
  # Reflection
  expect_equal(fixedopt$pot, 
               linregPotential(dim = dim(dat), 
                               range1 = c(2.5,6.5), 
                               range2 = c(0,10), 
                               alpha = 2, 
                               beta = -3, 
                               var = 2))
})
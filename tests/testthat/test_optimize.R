library(dgRaph)
context("Optimizing")

test_that("Beta optimizer",{
  # Generate data
  set.seed(1)
  dat <- rbeta(1000, 4, 140)
  discrete_dat <- as.numeric(cut(dat, breaks = seq(0, 0.1, 0.001), labels = 1:100))
  expCounts <- matrix(table(c(discrete_dat, 1:100))-1, 1, 100)
  
  # Optimize
  betaOpt <- betaOptimize(range = c(0,0.1))(expCounts)
  pot <- betaOpt$pot
  
  # Tests
  expect_less_than(
    max(abs(pbeta(1:100/1000, 4, 140) - cumsum(pot[1,]))),
    0.005)
  
  expect_true( 4*0.95 < betaOpt$alphas[1] && betaOpt$alphas[1] < 4*1.05 )
  expect_true( 140*0.95 < betaOpt$betas[1] && betaOpt$betas[1] < 140*1.05 )
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
  
  # Optimize
  # Bins x corresponds to 2,4,6,8
  # Bins y corresponds to 3,7,11,15,19
  linreg <- linregOptimize(range1 = c(1,9), range2 = c(1,21))(dat)
  
  # Tests
  expect_equal(linreg$alpha, -2)
  expect_equal(linreg$beta, 19)
  expect_equal(linreg$var, 8)
})

test_that("Normal optimizer",{
  # Generate data
  dat <- matrix(c(0,0,1,2,1,0), 1, 6)
  
  # Optimize
  normopt <- normOptimize(range = c(1,13))(dat)
  
  # Tests
  expect_equal(normopt$means[1], 8)
  expect_equal(normopt$vars[1], 2)
})
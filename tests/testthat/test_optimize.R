library(dgRaph)
context("Optimizing")

test_that("Beta optimizer",{
  # Generate data
  set.seed(1)
  dat <- rbeta(1000, 4, 140)
  expCounts <- matrix( hist(dat, breaks = seq(0, 0.1, 0.001), plot = FALSE)$counts, 1, 100)

  # Optimize
  betaPot <- betaPotential(dim = c(1,100), range = c(0,0.1))
  betaPotUpdated <- update(betaPot, expCounts)
  alpha <- betaPotUpdated$param$alphas[1]
  beta  <- betaPotUpdated$param$betas[1]
  
  # Tests
  expect_equal( alpha, 4, tolerance = 0.005)
  expect_equal( beta, 140, tolerance = 0.005)
  
  # Reflection
  expect_equal( betaPotUpdated, 
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
  linregPot <- linregPotential(dim = c(5,5), range1 = c(0.5,5.5), range2 = c(0.5,5.5))
  linregPotUpdated <- update(linregPot, dat)
  
  # Tests
  expect_equal(linregPotUpdated$param$alpha, -1)
  expect_equal(linregPotUpdated$param$beta, 5)
  expect_equal(linregPotUpdated$param$var, 0.5)
  
  # Reflection
  expect_equal(linregPotUpdated, 
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

test_that("Log regression 1",{
  # Generate data
  set.seed(1)
  expCounts <- t(sapply( log(seq(0.5,99.5,1)), FUN=function(x){
    dat <- rnorm(400, 10*x+20, 5)
    as.vector(unname(table(cut(dat, 0:100))))
  }))
  
  # Optimize
  reg <- .logregOptimize(expCounts)  
  
  # Tests
  expect_equal( reg$alpha, 10, tolerance = 0.005)
  expect_equal( reg$beta, 20, tolerance = 0.005)
  expect_equal( reg$var, 25, tolerance = 0.05)
  
  # Reflection
  expect_equal( reg$pot, 
                logregPotential(dim = c(100,100),
                                alpha = reg$alpha, 
                                beta = reg$beta,
                                var = reg$var
                                ))
})

test_that("Linreg/Normal optimizer comparison 1",{
  # Generate data
  dat <- matrix(c(0,0,1,2,1,0),6,6, byrow = T)
  
  # Optimize
  linregPot <- linregPotential(dim = c(6,6))
  linregPotUpdated <- update(linregPot, dat)
  
  # Rows identical
  expect_equal(as.matrix(linregPotUpdated)[1,], as.matrix(linregPotUpdated)[2,])
  
  # Check equal to similar one dimensional update
  datRow <- dat[1,,drop = F]
  normPot <- normalPotential(dim = c(1,6))
  normPotUpdated <- update(normPot,datRow)
  expect_equal( as.matrix(normPotUpdated)[1,], as.matrix(linregPotUpdated)[1,])
})

test_that("Normal optimizer 1",{
  # Generate data
  dat <- matrix(c(0,0,1,2,1,0,
                  0,0.5,4,0.5,0,0), 2, 6, byrow = T)
  
  # Optimize
  normPot <- normalPotential(range = c(1,13), dim = c(2,6))
  normPotUpdated <- update(normPot, dat)
  
  # Tests
  expect_equal(normPotUpdated$param$means, c(8,6) )
  expect_equal(normPotUpdated$param$vars, c(2,0.8) )
  expect_equal(which.max(as.matrix(normPotUpdated)[1,]), 4)
  expect_equal(which.max(as.matrix(normPotUpdated)[2,]), 3)
  
  # Reflection
  expect_equal(normPotUpdated, 
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
  normPot <- normalPotential(range = c(100,300), dim = c(1,200))
  normPotUpdated <- update(normPot, dat)
  
  # Tests
  expect_equal(normPotUpdated$param$means[1], 200)
  expect_equal(normPotUpdated$param$vars[1], 1.05) # Uses MLE for variance

  # Reflection
  expect_equal(normPotUpdated, 
               normalPotential(dim = c(1,200), 
                               range = c(100,300), 
                               means = 200, 
                               vars = 1.05))
})

test_that("Meanlink optimizer 1",{
  # Generate data
  dat <- matrix(c(1,1,0,0,0,
                  0,1,1,0,0,
                  0,0,1,1,0,
                  0,0,0,1,1), 4, 5, byrow = T)
  
  # Optimize
  meanlinkPot <- meanlinkPotential(dim = c(4,5),
                                   range1 = c(0.5,4.5),
                                   range2 = c(0.5,5.5),
                                   alpha = 1, beta = 0)
  meanlinkPotUpdated <- update(meanlinkPot, dat)
  
  # Tests
  expect_equal(meanlinkPotUpdated$param$var, 0.5)
  
  # Reflection
  expect_equal(meanlinkPotUpdated, 
               meanlinkPotential(dim = dim(dat), 
                                 range1 = c(0.5,4.5), 
                                 range2 = c(0.5,5.5), 
                                 alpha = 1, 
                                 beta = 0, 
                                 var = 0.5))
})

test_that("Meanlink optimizer 2",{
  # Generate data
  dat <- matrix(c(0,1,1,0,0,
                  0,0,0,0,0,
                  0,0,1,1,0,
                  0,0,0,0,0), 4, 5, byrow = T)
  
  # Optimize
  meanlinkPot <- meanlinkPotential(dim = c(4,5),
                                   range1 = c(2.5,6.5),
                                   range2 = c(0,10),
                                   alpha = 2,
                                   beta = -3)
  meanlinkPotUpdated <- update(meanlinkPot, dat)
  
  # Tests
  expect_equal(meanlinkPotUpdated$param$var, 2)
  
  # Reflection
  expect_equal(meanlinkPotUpdated, 
               meanlinkPotential(dim = dim(dat), 
                                 range1 = c(2.5,6.5), 
                                 range2 = c(0,10), 
                                 alpha = 2, 
                                 beta = -3, 
                                 var = 2))
})

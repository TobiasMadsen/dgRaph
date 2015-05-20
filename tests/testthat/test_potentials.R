library(dgRaph)
context("Potentials")

test_that("Beta potential",{
  # Generate a beta potential
  for(i in 1:10){
    set.seed(i)
    pot <- betaPotential(dim = c(2,100))
    expect_true(is.matrix(pot))
    expect_equal(dim(pot), c(2,100))
    expect_equal(rowSums(pot), c(1,1))
    expect_true(min(pot) >= 0)
  }
})

test_that("Multinomial potential",{
  # Generate a multinomial potential
  for(i in 1:10){
    set.seed(i)
    pot <- multinomialPotential(dim = c(2,4))
    expect_true(is.matrix(pot))
    expect_equal(dim(pot), c(2,4))
    expect_equal(rowSums(pot), c(1,1))
    expect_true(min(pot) >= 0)
  }
})

test_that("Normal potential", {
  # Generate a normal potential
  for(i in 1:10){
    set.seed(i)
    pot <- normalPotential(dim = c(3,100))
    expect_true(is.matrix(pot))
    expect_equal(dim(pot), c(3,100))
    expect_equal(rowSums(pot), c(1,1,1))
    expect_true(min(pot) >= 0)
  }
})

test_that("Get potentials",{
  source("cases/twoDepedentVariables.R")
  mydfg <- twoDependentVariables()
  
  pot <- potentials(mydfg)
  facPot <- c(list(matrix(c(0.7,0.3),1,2)),
              list(matrix(c(0.75,0.25,0.25,0.75),2,2)))
  
  expect_equal(pot, facPot)
})

test_that("Set potentials",{
  source("cases/twoDepedentVariables.R")
  mydfg <- twoDependentVariables()
  
  facPot <- c(list(matrix(c(0.5,0.5),1,2)),
              list(matrix(c(0.4,0.6,0.7,0.3),2,2)))
  mydfg$dfgmodule$resetPotentials( facPot )
  expect_equal( facPot, potentials(mydfg))
})
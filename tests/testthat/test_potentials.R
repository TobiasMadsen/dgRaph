library(dgRaph)
context("Potentials")

test_that("Beta potential",{
  # Generate a beta potential
  for(i in 1:3){
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
  for(i in 1:3){
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
  for(i in 1:3){
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
  
  newFacPot <- c(list(matrix(c(0.5,0.5),1,2)),
              list(matrix(c(0.4,0.6,0.7,0.3),2,2)))
  
  potentials(mydfg) <- newFacPot
  expect_equal( newFacPot, potentials(mydfg))
  expect_equal( newFacPot, mydfg$facPot)
})

test_that("Set potentials no update",{
  source("cases/twoDepedentVariables.R")
  mydfg <- twoDependentVariables()
  pot2 <- potentials(mydfg)[[2]]
  
  newFacPot <- c(list(matrix(c(0.5,0.5),1,2)),
                 list(matrix(0,0,0)))
  mydfg$dfgmodule$resetPotentials( newFacPot )
  expect_equal( newFacPot[[1]], potentials(mydfg)[[1]])
  expect_equal( pot2, potentials(mydfg)[[2]] )
})

test_that("Set potential error",{
  source("cases/twoDepedentVariables.R")
  mydfg <- twoDependentVariables()
  
  newFacPot <- c(list(matrix(c(0.5,0.5),1,2)),
                 list(matrix(c(0.4,0.6,0.7,0.3),1,2)))
  
  expect_error({potentials(mydfg) <- newFacPot})
})
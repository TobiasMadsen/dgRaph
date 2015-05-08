library(dgRaph)
context("Potentials")

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
  mydfg$dfgmodule$resetFactorPotentials( facPot )
  expect_equal( facPot, potentials(mydfg))
})
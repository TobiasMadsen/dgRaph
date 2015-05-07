library(dgRaph)
context("Likelihood Calculation")

test_that("Likelihood Calculation 1",{
  varDim <- rep(2L, 2)
  facPot <- c(list(matrix(c(0.5, 0.5),1,2)),
              list(matrix(c(0.75,0.25,0.25,0.75),2,2)))
  facNbs <- c(list(c(1L)),
              list(c(1L,2L)) )
  
  mydfg <- dfg(varDim, facPot, facNbs)
  
  data <- data.frame(O1 = c(1, 2, NA, NA), O2 = c(NA, 1, NA, 2))
  
  expect_equal(likelihood(data = data, dfg = mydfg), c(0.5, 0.125, 1, 0.5))
})
library(PGMscore)
context("PGM class")

test_that("Check potential dimensions",{
  varDim <- rep(2:5)
  facPot <- c(list(matrix(1,1,2)),
              list(matrix(1,2,3)),
              list(matrix(1,2,4)),
              list(matrix(1,4,5)))
  facNbs <- c(list(c(1L)),
              list(c(1L,2L)),
              list(c(1L,3L)),
              list(c(3L,4L)))
  
  expect_true({pgm <- PGM(varDim, facPot, facNbs); TRUE}) #If first commmand fails second command won't be evaluated
})

test_that("Wrong potential dimensions",{
  varDim <- rep(2:5)
  facPot <- c(list(matrix(1,1,2)),
              list(matrix(1,2,3)),
              list(matrix(1,2,4)),
              list(matrix(1,4,4)))
  facNbs <- c(list(c(1L)),
              list(c(1L,2L)),
              list(c(1L,3L)),
              list(c(3L,4L)))
  
  expect_error({pgm <- PGM(varDim, facPot, facNbs)}) #If first commmand fails second command won't be evaluated 
})

test_that("Too many neighbors",{
  varDim <- rep(2:5)
  facPot <- c(list(matrix(1,1,2)),
              list(matrix(1,2,3)),
              list(matrix(1,2,4)),
              list(matrix(1,4,5)))
  facNbs <- c(list(c(1L)),
              list(c(1L,2L)),
              list(c(1L,3L)),
              list(c(3L,4L,5L)))
  
  expect_error({pgm <- PGM(varDim, facPot, facNbs)}) #If first commmand fails second command won't be evaluated
})

test_that("Check acyclic",{
  varDim <- c(2,2)
  facPot <- c(list(matrix(1,2,2)),
              list(matrix(1,2,2)))
  facNbs <- c(list(c(1,2)),
              list(c(2,1)))
  
  expect_error({pgm <- PGM(varDim, facPot, facNbs)})
})

test_that("is.pgm",{
  x <- list()
  
  class(x) <- c('pgm')
  expect_true(is.pgm(x))
  
  class(x) <- c('troll','pgm')
  expect_true(is.pgm(x))
})

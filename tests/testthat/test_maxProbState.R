library(PGMscore)
context("Max probable state")

test_that("Dependent variables R level",{
  varDim <- rep(2,2)
  facPot <- c(list(matrix(c(0.7,0.3),1,2)),
              list(matrix(c(0.75,0.25,0.25,0.75),2,2)))
  facNbs <- c(list(c(1L)),
              list(c(1L,2L)))
  mydfg <- dfg(varDim, facPot, facNbs, varNames = c('x', 'y'))
  
  df <- data.frame(x = c(1,0,NA,NA,0), y = c(NA,NA,NA,1,1))
  df_mps <- data.frame(x = c(1,0,0,1,0), y = c(1,0,0,1,1))
  df_calc_mps <- mps(df, mydfg)
  expect_equal(df_calc_mps, df_mps)
})

test_that("Independent variables R level",{
  varDim <- rep(4,4)
  facPot <- c(list(matrix(c(0.05,0.05,0.20,0.70),1,4)),
              list(matrix(c(0.05,0.05,0.70,0.20),1,4)),
              list(matrix(c(0.70,0.20,0.05,0.05),1,4)),
              list(matrix(c(0.05,0.70,0.20,0.05),1,4)))
  facNbs <- c(list(c(1L)),
              list(c(2L)),
              list(c(3L)),
              list(c(4L)))
  
  mydfg <- dfg(varDim, facPot, facNbs, varNames = c('x', 'y', 'z', 'w'))
  
  df <- data.frame(x = 0:2, y = c(NA, 1, NA), z = 2:0, w = rep(1,3))
  df_mps <- data.frame(x = 0:2, y = c(2, 1, 2), z = 2:0, w = rep(1,3))
  df_calc_mps <- mps(df, mydfg)
  expect_equal(df_calc_mps, df_mps)
})

test_that("Independent variables cpp level",{
  varDim <- rep(4,4)
  facPot <- c(list(matrix(c(0.05,0.05,0.20,0.70),1,4)),
              list(matrix(c(0.05,0.05,0.70,0.20),1,4)),
              list(matrix(c(0.70,0.20,0.05,0.05),1,4)),
              list(matrix(c(0.05,0.70,0.20,0.05),1,4)))
  facNbs <- c(list(c(1L)),
              list(c(2L)),
              list(c(3L)),
              list(c(4L)))
  
  mydfg <- dfg(varDim, facPot, facNbs)
  
  mps <- mydfg$dfgmodule$maxProbState(integer(0), logical(0))
  expect_equal( mps, c(3,2,0,1)) 
})

test_that("Independent variables cpp level with observed variables",{
  varDim <- rep(4,4)
  facPot <- c(list(matrix(c(0.05,0.05,0.20,0.70),1,4)),
              list(matrix(c(0.05,0.05,0.70,0.20),1,4)),
              list(matrix(c(0.70,0.20,0.05,0.05),1,4)),
              list(matrix(c(0.05,0.70,0.20,0.05),1,4)))
  facNbs <- c(list(c(1L)),
              list(c(2L)),
              list(c(3L)),
              list(c(4L)))
  
  mydfg2 <- dfg(varDim, facPot, facNbs)
  mps <- mydfg2$dfgmodule$maxProbState(c(0,0,1,0), c(F,F,T,F))
  expect_equal( mps, c(3,2,1,1) )
})
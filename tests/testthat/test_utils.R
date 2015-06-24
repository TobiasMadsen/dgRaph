library(dgRaph)
context("Utils")

test_that("commonMap 1",{
  expect_equal(.commonMap(c(1,2,3),c(1,2,3)), c(1,2,3))
  
  expect_equal(.commonMap(c(1,2,3),c(1,1,2)), c(1,2,3))
  
  expect_equal(.commonMap(c(2,2,1),c(1,1,2)), c(1,1,2))
  
  expect_equal(.commonMap(c(1,2,1),c(2,1,2)), c(1,2,1))
  
  expect_equal(.commonMap(c(2,1,2),c(1,2,1)), c(1,2,1))
  
  expect_equal(.commonMap(c(1,1,1,3),c(2,2,3,1)), c(1,1,2,3))
})
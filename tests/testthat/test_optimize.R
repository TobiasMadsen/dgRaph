library(dgRaph)
context("Optimizing")

test_that("Beta optimizer",{
  # Generate data
  set.seed(1)
  dat <- rbeta(1000, 4, 140)
  discrete_dat <- as.numeric(cut(dat, breaks = seq(0, 0.1, 0.001), labels = 1:100))
  expCounts <- matrix(table(c(discrete_dat, 1:100))-1, 1, 100)
  
  # Optimize
  pot <- betaOptimize(range = c(0,0.1))(expCounts)
  
  expect_less_than(
    max(abs(pbeta(1:100/1000, 4, 140) - cumsum(pot[1,]))),
    0.05)
})
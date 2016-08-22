library(dgRaph)
context("Kullback-Leibler")

test_that("Poibin Binomial I",{
  N <- 10000
  p <- 0.01
  
  res_SA <- poibinSA(200, rep(p, N))
  res_BI <- pbinom(199, size = N, prob = p, lower.tail = F, log.p = T)
  
  expect_lt( abs(res_SA - res_BI), log(1.01)) # Maximally 1% difference in prob.-space
  
  res_SA <- poibinSA(150, rep(p, N))
  res_BI <- pbinom(149, size = N, prob = p, lower.tail = F, log.p = T)
  
  expect_lt( abs(res_SA - res_BI), log(1.01))
})
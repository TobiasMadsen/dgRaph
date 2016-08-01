library(dgRaph)
context("Expect")

test_that("Components", {
  varDim <- rep(1L,2)
  facPot <- c(list(matrix(3,1,1)),
              list(matrix(2,1,1)))
  facExp <- c(list(matrix(1,1,1)),
              list(matrix(2,1,1)))
  facNbs <- c(list(c(1L)),
              list(c(2L)))
  
  mydfg <- dfg(varDim, facPot, facNbs)
  
  res <- unname( expect(mydfg, facExp) )
  expect_equal(unname(res[1]),6)
  expect_equal(unname(res[2]),18)
  
  res <- unname( expect2(mydfg, facExp) )
  expect_equal(res[1],6)
  expect_equal(res[2],18)
  expect_equal(res[3],54)
})

test_that("Simple phylogenetic model",{
  varDim <- rep(2L,4)
  facPot <- c(list(matrix(c(0.8,0.2),1,2)),
              rep(list(matrix(c(0.8,0.2,0.2,0.8),2,2)),3))
  facExp <- c(list(matrix(0,1,2)),
              rep(list(matrix(c(0,1,1,0),2,2)),3))
  facNbs <- c(list(c(1L)),
              list(c(1L,2L)),
              list(c(1L,3L)),
              list(c(3L,4L)))
  mydfg <- dfg(varDim, facPot, facNbs)
  
  res <- unname( expect(mydfg, facExp) )
  expect_equal(length(res),2)
  expect_equal(res[1],1)
  expect_equal(res[2],0.6)
})

test_that("Bernoulli connected to ancestor",{
  N <- 15
  t <- 0.25
  varDim <- rep(2,N+1)
  facPot1 <- c( list(matrix(c(0.5,0.5),1,2)), rep(list(matrix(c(0.50,0.50,0.50,0.50),2,2)),N) )
  facPot2 <- c( list(matrix(c(0.5,0.5),1,2)), rep(list(matrix(c(0.75,0.25,0.25,0.75),2,2)),N) )
  facNbs <- c(list(c(1L)), lapply(1:N, function(x){c(1,x+1)}))
  
  facPot <- .facPotToFunA(facPot1,facPot2,t)
  facExp <- .facPotToFunB(facPot1,facPot2)
  
  mydfg <- dfg(varDim, facPot, facNbs)
  res <- unname(expect2(mydfg, facExp))
  
  mgf <- function(t){(((1/2)**t+(3/2)**t)/2)**N}
  dtlogmgf <- function(t){N/((1/2)**t+(3/2)**t)*(log(1/2)*(1/2)**t+log(3/2)*(3/2)**t )}
  dt2logmgf <- function(t){ N*((log(2)**2/2**t+log(3/2)**2*(3/2)**t)/(1/2**t+(3/2)**t)
                               -(-log(2)/2**t+log(3/2)*(3/2)**t)**2/(3**t/2**t+1/2**t)**2) }

  expect_equal(mgf(t),res[1])
  expect_equal(dtlogmgf(t),res[2]/res[1])
  expect_equal(dt2logmgf(t), (res[3]*res[1]-res[2]**2)/(res[1]**2))
})

test_that("Bernoulli disconnected",{
  N <- 4
  varDim <- rep(2L,N)
  facPot1 <- rep(list(matrix(c(0.5,0.5),1,2)),N)
  facPot2 <- rep(list(matrix(c(0.75,0.25),1,2)),N)
  facNbs <- lapply(1:N,function(i){c(i)})
  t <- 0.25
  
  facPot <- .facPotToFunA(facPot1,facPot2,t)
  facExp <- .facPotToFunB(facPot1,facPot2)
  
  mydfg <- dfg(varDim, facPot, facNbs)
  res <- unname(expect2(mydfg, facExp))
  
  mgf <- function(t){(((1/2)**t+(3/2)**t)/2)**N}
  dtlogmgf <- function(t){N/((1/2)**t+(3/2)**t)*(log(1/2)*(1/2)**t+log(3/2)*(3/2)**t )}
  dt2logmgf <- function(t){ N*((log(2)**2/2**t+log(3/2)**2*(3/2)**t)/(1/2**t+(3/2)**t)
                               -(-log(2)/2**t+log(3/2)*(3/2)**t)**2/(3**t/2**t+1/2**t)**2) }
  
  expect_equal(mgf(t),res[1])
  expect_equal(dtlogmgf(t), res[2]/res[1])
  expect_equal(dt2logmgf(t), (res[3]*res[1]-res[2]**2)/(res[1]**2))
})

#################################################
# Second order moment
#################################################

test_that("Constant variable", {
  varDim <- c(1L)
  facPot <- list(matrix(3,1,1))
  facExp <- list(matrix(2,1,1))
  facNbs <- list(1L)
  
  mydfg <- dfg(varDim, facPot, facNbs)
  
  res <- unname( expect2(mydfg, facExp) )
  expect_equal(res[1], 3)
  expect_equal(res[2], 6)
  expect_equal(res[3], 12)
})

test_that("Pair of variables",{
  varDim <- rep(2L,2)
  facPot <- list(matrix(1:2,1,2),
                 matrix(1:4,2,2, byrow = T))
  facExp <- list(matrix(1,1,2),
                 matrix(1,2,2))
  facNbs <- list(1,1:2)
  mydfg <- dfg(varDim, facPot, facNbs)
  
  res <- unname(expect2(mydfg, facExp))
  
  expect_equal(res[1], 17)
  expect_equal(res[2], 34)
  expect_equal(res[3], 68)
})

test_that("Multiple neighbors",{
  varDim <- 1L
  facPot <- lapply(1:8,function(i){matrix(1)})
  facExp <- facPot
  facNbs <- list(1,1,1,1,1,1,1,1)
  mydfg <- dfg(varDim, facPot, facNbs)
  
  res <- unname(expect2(mydfg, facExp))
  expect_equal(res[1], 1)
  expect_equal(res[2], 8)
  expect_equal(res[3], 64)
})

test_that("Components 2",{
  varDim <- rep(2L, 2)
  facPot <- list(matrix(1:2,1,2),
                 matrix(3:4,1,2))
  facExp <- facPot
  facNbs <- list(1,2)
  
  mydfg <- dfg(varDim, facPot, facNbs)
  res <- unname(expect2(mydfg, facExp))
  
  expect_equal(res[1], 21)
  expect_equal(res[2], 110)
  expect_equal(res[3], 586)
})

test_that("Larger structure I",{
  varDim <- rep(2L, 7)
  facPot <- list(matrix(1:2/3,1,2), matrix( c(1,3,5,3)/6,2,2), matrix( c(1,1,2,2)/3,2,2) )
  facExp <- list(matrix(0:1,1,2), matrix( 0,2,2), matrix( c(0,0,1,1),2,2) )
  facNbs <- list(1, c(1,2), c(1,3), c(1,4), c(2,5), c(3,6), c(4,7))
  potMap <- c(1,2,2,2,3,3,3)
  
  mydfg <- dfg(varDim, facPot, facNbs, potMap)
  res <- unname(expect2(mydfg, facExp))
  
  expect_equal(res[1], 1)
  expect_equal(res[2], 8/3)
  expect_equal(res[3], 8)
  
  varDim <- rep(2L, 7)
  facPot <- list(matrix(1:2/3,1,2)/100, matrix( c(1,3,5,3)/6,2,2)/100, matrix( c(1,1,2,2)/3,2,2)/100 )
  facExp <- list(matrix(0:1,1,2), matrix( 0,2,2), matrix( c(0,0,1,1),2,2) )
  facNbs <- list(1, c(1,2), c(1,3), c(1,4), c(2,5), c(3,6), c(4,7))
  potMap <- c(1,2,2,2,3,3,3)
  
  mydfg <- dfg(varDim, facPot, facNbs, potMap)
  res <- unname(expect2(mydfg, facExp))
  
  expect_equal(res[1], 1e-14)
  expect_equal(res[2], 8e-14/3)
  expect_equal(res[3], 8e-14)
})

test_that("Larger structure II", {
  varDim <- rep(2L, 7)
  facPot1 <- list(matrix(1:2/3,1,2), matrix( c(1,1,2,2)/3,2,2) )
  facPot2 <- list(matrix(3:4/7,1,2), matrix( c(3,3,4,4)/7,2,2) )
  facNbs <- list(1, c(1,2), c(1,3), c(1,4), c(2,5), c(3,6), c(4,7))
  potMap <- c(1,2,2,2,2,2,2)
  
  x <- -3
  facPot <- lapply(seq_along(facPot1), function(i){
    facPot1[[i]]*exp(x*(log(facPot2[[i]])-log(facPot1[[i]])))
  })
  facExp <- lapply(seq_along(facPot1), function(i){
    log(facPot2[[i]])-log(facPot1[[i]])
  })
  
  mydfg <- dfg(varDim, facPot, facNbs, potMap)
  res <- unname(expect2(mydfg, facExp))
  
  mgf <- function(x){
    b <- log(6/7)*7
    a <- log(3/2)
    p <- 1/3
    n <- 7
    ret <- rep(0,3)
    ret[1] <- exp(x*b)*(1-1/3+exp(x*a)/3)**7
    ret[2] <- a*n*p*(p*exp(1)^(a*x)-p+1)^(n-1)*exp(1)^(b*x+a*x)+b*(p*exp(1)^(a*x)-p+1)^n*exp(1)^(b*x)
    ret[3] <- a^2*(n-1)*n*p^2*(p*exp(1)^(a*x)-p+1)^(n-2)*exp(1)^(b*x+2*a*x)+
      a*(b+a)*n*p*(p*exp(1)^(a*x)-p+1)^(n-1)*exp(1)^(b*x+a*x)+
      a*b*n*p*(p*exp(1)^(a*x)-p+1)^(n-1)*exp(1)^(b*x+a*x)+b^2*(p*exp(1)^(a*x)-p+1)^n*exp(1)^(b*x)
    ret
  }
  
  analytic_res <- mgf(x)
  expect_equal(analytic_res[1], res[1])
  expect_equal(analytic_res[2], res[2])
  expect_equal(analytic_res[3], res[3])
})

test_that("Larger structure III", {
  varDim <- rep(2L, 3)
  facPot1 <- list(matrix(1:2,1,2)/3, matrix( c(1,1,2,2)/3,2,2) )
  facNbs <- list(1, c(1,2), c(2,3))
  potMap <- c(1,2,2)
  
  mydfg <- dfg(varDim, facPot1, facNbs, potMap)
  
  facExp <- list(matrix(0, 1,2), matrix(1:2,2,2,byrow = T))
  
  res <- unname(expect2(mydfg, facExp))
  expect_equal(res[1], 1)
  expect_equal(res[2], 10/3)
  expect_equal(res[3], 104/9)
})

test_that("Larger structure IV",{
  varDim <- rep(2L, 4)
  facPot1 <- list(matrix(1:2,1,2), matrix( c(1,1,2,2),2,2) )
  facNbs <- list(1, c(1,2), c(1,3), c(2,4))
  potMap <- c(1,2,2,2)
  
  mydfg <- dfg(varDim, facPot1, facNbs, potMap)
  
  facExp <- list(matrix(1:2, 1,2), matrix(1:2,2,2,byrow = T))
  
  res <- unname(expect2(mydfg, facExp))
  
  expect_equal(res[1], 3**4)
  expect_equal(res[2], 20/3*3**4)
  expect_equal(res[3], 408/9*3**4)
})

test_that("Unobserved data", {
  varDim <- c(2,2)
  facPot <- list(matrix(0.5, 1,2), matrix(0.5, 2,2))
  facNbs <- list(1, c(1,2))

  mydfg <- dfg(varDim, facPot, facNbs)

  facExp <- list(matrix(0,1,2), matrix(c(1,1,1,0), 2,2))

  res <- unname(expect(mydfg, facExp))

  expect_equal(res[1], 1)
  expect_equal(res[2], 0.75)
  
  res <- unname(expect.dfg(mydfg, facExp, data = matrix(c(NA, 1, 2, 2), 2, 2, byrow = T)))
  
  expect_equal(res[1,2] / res[1,1], 1)
  expect_equal(res[2,2] / res[2,1], 0)
})

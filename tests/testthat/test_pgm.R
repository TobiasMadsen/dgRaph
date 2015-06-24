library(dgRaph)
context("PGM class")

test_that("Remapping potentials 1",{
  varDim <- c(2,2)
  facPot <- list(matrix(1,1,2),
                 matrix(1,2,2))
  facNbs <- list(c(1,2),1)
  potMap <- c(2,1)
  
  mydfg <- dfg(varDim, facPot, facNbs, potMap)
  mydfgRemap <- .remapPotMapDfg(mydfg, c(1,2))
  
  expect_equal( mydfgRemap$potMap, c(1,2))
  expect_equal( mydfg$facPot[[1]], mydfgRemap$facPot[[2]])
  expect_equal( mydfg$facPot[[2]], mydfgRemap$facPot[[1]])
})

test_that("Remapping facNbs 1",{
  varDim <- c(2,2,2)
  facPot <- list(matrix(1,1,2),
                 matrix(1,2,2),
                 matrix(1,2,2))
  facNbs <- list(c(1),c(1,2),c(2,3))
  
  mydfg <- dfg(varDim, facPot, facNbs)
  mydfgRemap <- .remapFacNbsDfg(mydfg, c(2,3,1))
  
  expect_equal( mydfgRemap$potMap[2], mydfg$potMap[1])
  expect_equal( mydfgRemap$potMap[3], mydfg$potMap[2])
  expect_equal( mydfgRemap$potMap[1], mydfg$potMap[3])
  
  expect_equal( mydfgRemap$facNbs[[2]], mydfg$facNbs[[1]])
  expect_equal( mydfgRemap$facNbs[[3]], mydfg$facNbs[[2]])
  expect_equal( mydfgRemap$facNbs[[1]], mydfg$facNbs[[3]])
  
  # Going back
  mydfgRemap2 <- .remapFacNbsDfg(mydfgRemap, match(mydfgRemap$facNbs, mydfg$facNbs))
  expect_equal( mydfgRemap2$facNbs, mydfg$facNbs)
  expect_equal( mydfgRemap2$potMap, mydfg$potMap)
})

test_that("Variable and factor names",{
  varDim <- rep(4,6)
  facPot <- list(matrix(0.25,1,4),
                 matrix(0.25,4,4),
                 matrix(0.25,4,4),
                 matrix(0.25,4,4),
                 matrix(0.25,4,4),
                 matrix(0.25,4,4))
  facNbs <- c(list(c(1L)),
              list(c(1L,2L)),
              list(c(1L,3L)),
              list(c(3L,4L)),
              list(c(3L,5L)),
              list(c(5L,6L)) )
  facNames = c("Prior",rep("Int",5))
  varNames = c("Foo","Bar","Baz","Do","Re","Mi")
  expect_true({mydfg <- dfg(varDim, facPot, facNbs, varNames = varNames, facNames = facNames); TRUE})
})

test_that("Shared potentials",{
  # Markov chain
  N <- 1000
  varDim <- rep(2, N)
  facPot <- c(list(matrix(0.5, 1, 2)),
              list(matrix(c(0.6,0.4,0.7,0.3), 2, 2)))
  potMap <- c(1, rep(2, N-1))
  facNbs <- c(list(c(1L)), lapply(2:N, FUN=function(i){c(i-1, i)}))
  expect_true({mydfg <- dfg(varDim, facPot, facNbs, potMap = potMap); TRUE})
})

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
  
  expect_true({mydfg <- dfg(varDim, facPot, facNbs); TRUE}) #If first commmand fails second command won't be evaluated
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
  
  expect_error({mydfg <- dfg(varDim, facPot, facNbs)})
})

test_that("Wrong potential dimensions",{
  varDim <- c(2,2)
  facPot <- c(list(matrix(0,2,2)),
              list(matrix(0,2,2)))
  facNbs <- c(list(1),
              list(c(1,2)))
  
  expect_error({mydfg <- dfg(varDim, facPot, facNbs)}) 
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
  
  expect_error({mydfg <- dfg(varDim, facPot, facNbs)}) 
})

test_that("Check acyclic",{
  varDim <- c(2,2)
  facPot <- c(list(matrix(1,2,2)),
              list(matrix(1,2,2)))
  facNbs <- c(list(c(1,2)),
              list(c(2,1)))
  
  expect_error({mydfg <- dfg(varDim, facPot, facNbs)})
})

test_that("Check acyclic two components", {
  varDim <- c(2,2,2)
  facPot <- c(list(matrix(1,1,2)),
              list(matrix(1,2,2)),
              list(matrix(1,2,2)))
  facNbs <- c(list(c(1)),
              list(c(1,2)),
              list(c(2,1)))
  
  expect_error({mydfg <- dfg(varDim, facPot, facNbs)})
})

test_that("is.dfg",{
  x <- list()
  
  class(x) <- c('dfg')
  expect_true(is.dfg(x))
  
  class(x) <- c('troll','dfg')
  expect_true(is.dfg(x))
})

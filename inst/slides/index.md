---
title       : Factor Graphs
subtitle    : More fun with factor graphs
author      : Tobias Madsen
job         : 
framework   : revealjs       # {io2012, html5slides, shower, dzslides, ...}
revealjs: {theme: Solarized, transition: none, center: "true"}
hitheme     : zenburn      # 
bootstrap:
  theme: amelia
widgets     : mathjax            # {mathjax, quiz, bootstrap}
mode        : selfcontained # {standalone, draft}
---

<script src="http://ajax.googleapis.com/ajax/libs/jquery/1.9.1/jquery.min.js"></script>



# Factorgraph fun
### Tobias Madsen
<img src="../fig/AU_logo.png"  alt="some_text" width='50%' style="background:none; border:none; box-shadow:none;">

--- &vertical

## Discrete Factor Graphs

- Class of graphical models
- Specify relationship among variables

***

<img src="../fig/mixture_fac.svg"  alt="some_text" width='50%' margin='100px 15% 0% 15%' style="background:none; border:none; box-shadow:none;">

- Consists of variables and factors
- Factors have potentials
- The likelihood is the product of potentials

***

<img src="../fig/mixture_fac.svg"  alt="some_text" width='50%' margin='100px 15% 0% 15%' style="background:none; border:none; box-shadow:none;">

$$
p(x) = f_a(x_1)f_b(x_2)
$$

--- &vertical

## Known models
* Cast known models into DFGs

***

Hidden Markov Model

<img src="../fig/hmm.svg"  alt="some_text" width='50%' margin='100px 15% 0% 15%' style="background:none; border:none; box-shadow:none;">

***

Phylogenetic Models

<img src="../fig/phylo.svg"  alt="some_text" width='50%' margin='100px 15% 0% 15%' style="background:none; border:none; box-shadow:none;">


--- &vertical

## Local Examples
- Michals normal/cancer model
- SCFG
- Augmented motifs
- Malene, Johanna et al

***

# Micha&#322;
<img src="../fig/michal.png"  alt="some_text" width='50%' margin='100px 15% 0% 15%'>

***

# Sudhakar
<img src="../fig/sudhakar.png"  alt="some_text" width='50%' margin='100px 15% 0% 15%'>

***

# Myself
<img src="../fig/cons.svg" alt="some_text" width='50%' margin='100px 15% 0% 15%'>

***

# Myself

<img src="../fig/pwm.png"  alt="some_text" width='50%' margin='100px 15% 0% 15%'>

--- &vertical

## R Interface
### NB! Very imature

***

## Specifying factor graph

<!---
NB! Caching code with Rcpp modules causes crashes
-->


```r
library(dgRaph)
varDim <- rep(4,6)
facPot <- c(list(matrix(0.25,1,4)),
            list(matrix(0.25,4,4)),
            list(matrix(0.25,4,4)),
            list(matrix(0.25,4,4)),
            list(matrix(0.25,4,4)),
            list(matrix(0.25,4,4)))
facNbs <- c(list(c(1L)),
            list(c(1L,2L)),
            list(c(1L,3L)),
            list(c(3L,4L)),
            list(c(3L,5L)),
            list(c(5L,6L)) )
facNames = c("Prior",rep("Int",5))
varNames = c("Foo","Bar","Baz","Do","Re","Mi")
mydfg <- dfg(varDim, facPot, facNbs, varNames, facNames)
```

***

## Plotting


```r
plot(mydfg)
```

![plot of chunk unnamed-chunk-3](assets/fig/unnamed-chunk-3-1.png) 

---

## Most Probable State


```r
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
mydfg2$dfgmodule$maxProbState(list(),c(0,0,1,0), c(F,F,T,F))+1
```

```
## Error in eval(expr, envir, enclos): could not find valid method
```

---

## Factor Expectation Counts

```r
varDim <- rep(2L, 2)
facPot <- c(list(matrix(c(0.5, 0.5),1,2)),
            list(matrix(c(0.75,0.25,0.25,0.75),2,2)))
facNbs <- c(list(c(1L)),
            list(c(1L,2L)) )

mydfg3 <- dfg(varDim, facPot, facNbs)

mydfg3$dfgmodule$facExpCounts(matrix(c(NA,1,
                                       0,1,
                                       1,NA), 3, 2) )
```

```
## [[1]]
##      [,1] [,2]
## [1,] 1.25 1.75
## 
## [[2]]
##      [,1] [,2]
## [1,] 0.75 0.50
## [2,] 0.00 1.75
```

--- &vertical

## Comparing two models
Use bayes factor, the ratio between likelihoods

$$
K = \frac{P(X|M_2)}{P(X|M_1)}
$$

***

## Tail Approximations

Same as comparing loglikelihood-ratio 

$$
S(x) = \log \frac{\prod_{\mathcal{A}}f_a^\prime(x_a)}{\prod_{\mathcal{A}}f_a(x_a)} = \sum_{\mathcal{A}} \left[\log f_a^\prime(x_a) - \log f_a(x_a)\right]
$$

***

## Finding p-value
Assess the FDR given a score treshold
 
> - naive sampling
> - importance sampling
> - saddlepoint approximation

--- &vertical

## Naive Sampling

***


```r
varDim <- rep(4,6)
facPot <- c(list(matrix(0.25,1,4)),
            list(matrix(0.25,4,4)),
            list(matrix(0.25,4,4)),
            list(matrix(0.25,4,4)),
            list(matrix(0.25,4,4)),
            list(matrix(0.25,4,4)))
facNbs <- c(list(c(1L)),
            list(c(1L,2L)),
            list(c(1L,3L)),
            list(c(3L,4L)),
            list(c(3L,5L)),
            list(c(5L,6L)) )
mydfg <- dfg(varDim, facPot, facNbs)
```

***


```r
plot(mydfg)
```

![plot of chunk unnamed-chunk-7](assets/fig/unnamed-chunk-7-1.png) 

***


```r
rndDist <- function(row,col){
  x <- rexp(row*col,1/(c(1:row)+2))
  mat <- matrix(x, row, col, byrow = T)
  mat/rowSums(mat)
}
set.seed(1)
facPotFg <- c(list(rndDist(1,4)),
              list(rndDist(4,4)),
              list(rndDist(4,4)),
              list(rndDist(4,4)),
              list(rndDist(4,4)),
              list(rndDist(4,4)))
```

***


```r
x <- seq(-3,4,0.01)
dfnaive <- tailIS(x, n=1000, alpha=0, dfg=mydfg, facPotFg=facPotFg)
```

***


```r
ggplot(dfnaive, aes(x=x,y=p)) + geom_line() + theme_bw() +
  scale_y_log10() + 
  geom_ribbon(aes(ymin=pmax(low,0.0001),ymax=high),alpha=0.3,fill="blue") +
  annotation_logticks(sides="l")
```

![plot of chunk unnamed-chunk-10](assets/fig/unnamed-chunk-10-1.png) 

***

## Naive sampling

Naive sampling s**cks at estimating tail probabilities because

.fragment <b>Rare events are rare</b>

--- &vertical

## Importance sampling

- Sample from another distribution
- Reweight samples

***




```r
#True value
pnorm(4, 0, 1, lower.tail = F)
```

```
## [1] 3.16712418331e-05
```

```r
#Naive Sampling
naive <- rnorm(100000, 0,1) > 4
mean(naive)
```

```
## [1] 1e-05
```

```r
#Importance Sampling
is    <- rnorm(100000, 5, 1)
is    <- (dnorm(is, 0, 1)/dnorm(is, 5, 1))*(is>4)
mean(is)
```

```
## [1] 3.15910435059e-05
```

***

## Importance Sampling Distributions

Use a class of distributions parameterized by $\alpha$

$$
P^{IS}(X) = P_{M_1}(X)^{(1-\alpha)}P_{M_2}(X)^{\alpha}
$$

***


```r
x <- seq(-3,4,0.01)
dfis <- tailIS(x, n=1000, alpha=1.5, dfg=mydfg, facPotFg=facPotFg)
```

***


```r
ggplot(dfis, aes(x=x,y=p)) + geom_line() + theme_bw() +
  scale_y_log10() + 
  geom_ribbon(aes(ymin=pmax(low,0.0001),ymax=high),alpha=0.3,fill="red") +
  annotation_logticks(sides="l")
```

![plot of chunk unnamed-chunk-14](assets/fig/unnamed-chunk-14-1.png) 

***

### Combining alphas


```r
x <- seq(-3,4,0.01)
dfis_combined <- tailIS(x, n=1000, alpha=c(0,0.5,1.5), dfg=mydfg, facPotFg=facPotFg)
```

***


```r
ggplot(dfis_combined, aes(x=x,y=p)) + geom_line() + theme_bw() +
  scale_y_log10() + 
  geom_ribbon(aes(ymin=pmax(low,0.0001),ymax=high),alpha=0.3,fill="red") +
  annotation_logticks(sides="l")
```

![plot of chunk unnamed-chunk-16](assets/fig/unnamed-chunk-16-1.png) 

--- &vertical

## Saddlepoint approximation

***

### Moment generating function

$$
P(S > s) = 
$$

***


```r
x <- seq(-3,4,0.01)
dfsaddle <- tailSaddle(x, mydfg, facPotFg)
```

***


```r
ggplot(dfsaddle, aes(x=x,y=p)) + geom_line() + theme_bw() + 
  scale_y_log10() +
  ylab("p-value") + xlab("score") +
  ggtitle("Saddlepoint approximation of the upper tail")
```

![plot of chunk unnamed-chunk-18](assets/fig/unnamed-chunk-18-1.png) 

---

### Comparison of all three methods


![plot of chunk unnamed-chunk-19](assets/fig/unnamed-chunk-19-1.png) 

---

## Acknowledgements
<img src="../fig/finish.jpg"  alt="some_text" width='50%' margin='100px 15% 0% 15%' style="background:none; border:none; box-shadow:none;">
### Jakob Skou
### Jens Ledet
### Asger Hobolth


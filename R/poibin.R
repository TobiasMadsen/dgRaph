#' Poison binomial tail probabilities using Saddlepoint Approximation
#'
#' Let \eqn{X_1,X_2,\ldots X_N} be independent Bernouilli trials
#' with \eqn{X_i \sim p_i}.
#' We evaluate the probability
#' \deqn{P(s_1X_1+s_2X_2\ldots s_NX_N > t)}.
#' 
#' @param t   Point where tail probability is evaluated (not this function is not vectorized over t)
#' @param p   Probabilities of each Bernouilli trial
#' @param s   Score associated with each trial
#' @param lattice Lattice size in minimal lattice
#' @param log Return log probability
#' 
#' @examples
#' poibinSA(300, rep(0.01, 10000))
#' pbinom(299, 10000, 0.01, lower.tail = F, log.p = T)
poibinSA <- function(t, p, s = rep(1, length(p)), lattice = 1L, log = T){
  stopifnot(length(t) == 1)
  stopifnot(t > 0)
  stopifnot(all(s > 0))
  # Giver saddelpunkts approksimationen paa log-skala
  # Use a combination of bisection and Newton raphson
  
  theta_old <- -Inf
  theta <- 0
  
  # 
  theta_largest_negative <- -Inf # Largest theta that gives a negative value so far
  theta_smallest_positive <- Inf # Smallest theta that gives a positive value so far
  
  iter <- 0
  while(abs(theta - theta_old) > 1e-8){
    # For convergence properties
    iter <- iter + 1
    theta_old <- theta
    
    # Newton raphson
    val <- cumulantD1(theta, p, s) - t
    if(val < 0 && theta > theta_largest_negative)
      theta_largest_negative <- theta
    if(val > 0 && theta < theta_smallest_positive)
      theta_smallest_positive <- theta
    deriv <- cumulantD2(theta, p, s)
    step <- - val / deriv 
    
    # Next theta if newton raphson
    theta <- theta + step
    if(theta < theta_largest_negative || theta > theta_smallest_positive){
      # Fall back to bisection
      theta <- (theta_largest_negative+theta_smallest_positive)/2
    }
  }
  
  # El approximazione
  v <- cumulantD2(theta, p, s)
  
  # Uden lattice-korrektion
  if(lattice == 0){
    ret <- cumulant(theta, p, s)-t*theta+theta^2*v/2+pnorm(theta*sqrt(v), lower.tail = F, log.p = T)
  } else{
    ret <- cumulant(theta, p, s)-t*theta+theta^2*v/2+pnorm(theta*sqrt(v), lower.tail = F, log.p = T) +
            log(abs(theta*lattice))- log((1-exp(-lattice*abs(theta))))
  }
  
  ifelse(log, ret, exp(ret))
}
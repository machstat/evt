# This file contains the following functions:
# gev.pplikMinima, gev.pprof

gev.pplikMinima <- function(x, theta, z0, p0) {
  # initial guesses
  sgm <- theta[1]
  xi <- theta[2]
  muTilde <- z0 + sgm/xi * (( - log(1 - p0))^( - xi) - 1)
  
  # normalizing
  y <- (x - muTilde)/sgm
  y <- 1 - xi * y
  
  # negative log-likelihood
  if(is.infinite(muTilde) || sgm <= 0 || any(y <= 0))
    l <- 10^6
  else l <- length(y) * log(sgm) + sum(y^(-1/xi)) + sum(log(
    y)) * (1/xi + 1)
  return(l)
}

gev.pprofMinima <- function(x, theta0, z0, p.low, p.up, nint = 1000, conf = 0.95, pplot = 0) {
  l <- numeric(nint)
  p <- seq(p.low, p.up, length = nint)
  
  # optimizing algorithm for the parameters
  for(i in 1:nint) {
    p0 <- p[i]
    opt <- optim(par = theta0, fn = gev.pplikMinima, x = x, p0 = p0, z0 = z0)
    theta0 <- opt$par 
    
    # negative log-likelihood for p
    l[i] <- opt$value 
  }
  
  # MLE of probability 
  mle.p <- p[which.max(-l)]
  
  # maximum value of log-likelihood function
  plik.max <- max(-l)
  
  # quantile of chi-square
  q <- plik.max - 0.5 * qchisq(conf, 1)
  
  # find  the  index  of -q in l for lower bound (LB) and upper bound (UB)
  plik.lb <- l[l <= -q ][1] ## bigger values
  ind.lb <- match(plik.lb, l) ## 21
  
  plik.ub <- l[l>=-q & p>mle.p][1]
  ind.ub <- match(plik.ub, l)
  
  # confidence interval
  ci.lb <- p[ind.lb] 
  ci.ub <- p[ind.ub] 
  
  # plotting of the values
  if (pplot){
    plot(p,  - l, type="l", xlab = "Probability", ylab = 
           " Profile Log-likelihood")
    abline(h = plik.max - 0.5 * qchisq(conf, 1), col = "red")
    abline(h = plik.max, col = 4)
    abline(v=ci.lb, col="red")
    abline(v=ci.ub, col="red")
  }
  
  return(data.frame("LB" = ci.lb,"UB" = ci.ub, "MLE" = mle.p))
}
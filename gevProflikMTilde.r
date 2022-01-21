gev.plikplot <- function(x, theta0, z0, p.low, p.up, nint, pplot=0) {
  l <- numeric(nint)
  p <- seq(p.low, p.up, length = nint)
  
  # Optimizing algorithm for the parameters
  for(i in 1:nint) {
    p0 <- p[i]
    opt <- optim(theta0, gev.plikMinima, x=x, p0=p0, z0=z0)
    theta0 <- opt$par 
    
    # negative log-likelihood for p
    l[i] <- opt$value 
  }
  
  # MLE of probability 
  mle.p <- p[which.max(-l)]
  
  # maximum value of log-likelihood function
  plik.max <- max(-l)
  
  # quantile of chi-square
  q <- plik.max- 0.5 * qchisq(0.95, 1)
  
  # find  the  index  of -q in l for lower bound (LB) and upper bound (UB)
  plik.LB <- l[l <= -q ][1] ## bigger values
  ind.LB <- match(plik.LB,l) ## 21
  
  plik.UB <- l[l>=-q & p>mle.p][1]
  ind.UB <- match(plik.UB,l)
  
  # Confidence interval
  CI.LB <- p[ind.LB] 
  CI.UB <- p[ind.UB] 
  
  # Plotting of the values
  if (pplot){
    plot(p,  - l, type="l", xlab = "Probability", ylab = 
           " Profile Log-likelihood")
    abline(h = plik.max - 0.5 * qchisq(0.95, 1), col = "red")
    abline(h = plik.max, col = 4)
    abline(v=CI.LB, col="red")
    abline(v=CI.UB, col="red")
  }
  
  return(data.frame("LB"=CI.LB,"UB"=CI.UB, "MLE"=mle.p))
}

gev.plikMinima <- function(x, theta, z0, p0) {
  sgm <- theta[1]
  xi <- theta[2]
  
  muTilde <- z0 + sgm/xi * (( - log(1 - p0))^( - xi) - 1)
  y <- (x - muTilde)/sgm
  y <- 1 - xi * y
  if(is.infinite(muTilde) || sgm <= 0 || any(y <= 0))
    l <- 10^6
  else l <- length(y) * log(sgm) + sum(y^(-1/xi)) + sum(log(
    y)) * (1/xi + 1)
  return(l)
}

dataBlock <- c(3.766, 0.000, 3.149, 3.041, 5.902, 3.906, 4.839, 4.032, 3.189,
               1.382, 3.120, 1.116, 2.423, 3.312, 1.828, 4.511, 4.223, 4.083, 3.292,
               4.263, 3.556, 5.415, 5.334, 1.407, 0.960, 3.652, 2.277, 2.743, 2.597,
               2.070, 3.449, 2.418, 3.902, 1.902, 4.149, 2.261,
               1.023, 1.511, 2.062, 2.458, 3.374, 2.498, 5.015, 3.635, 4.161)

p.low <- 0.000001
p.up <- 0.053

# Number of points between interval
nint <- 1000

# Profile likelihood
theta <- c(-3.562554, 1.278715, -0.264189) ## Estimated parameters
z0 <- 0

1-pgev(0,theta[1],theta[2],theta[3])

theta0 <- c(theta[2], theta[3])# Initial guesses for sigma and xi

gev.plikplot(x=dataBlock, theta0=theta0, z0, p.low=p.low, p.up=p.up, nint)


# This file contains the following functions:
# gpd.pplikMinima, gpd.pprof

gpd.plikMinima <- function(x, theta, u, eta, p0, z0) {
  # initial guesses
  xi <- theta
  sgm <- ((u - z0) * theta)*((p0 / eta)^(-theta) - 1)^(-1)
  
  # normalizing
  y <- (u - x) / sgm
  y <- 1 + theta * y
  
  # negative log-likelihood
  if(any(y <= 0))
    l <- 10^6
  else l <- -(length(x) * (log(eta) - log(sgm)) - (1/theta + 1) * sum(log(y)))
  l
}

gpd.pprofMinima <- function(x, theta0, u, z0, p.low, p.up, nint = 1000, conf = 0.95, pplot = 0){
  l <- numeric(nint)
  p <- seq(p.low, p.up, length = nint)
  
  theta <- theta0[2]
  eta <- sum(x < u)/length(x)
  x <- x[x < u]
  
  # optimizing algorithm for the parameters
  for (i in 1:nint) {
    p0 <- p[i]
    opt <- optim(theta, fn=gpd.plikMinima, method = "BFGS", x = x, u = u, eta = eta, p0 = p0, z0 = z0)
    theta <- opt$par
    
    # negative log -likelihood  for p
    l[i] <- opt$value
  }
  
  # MLE of  probability
  mle.p <- p[which.max(-l)]
  
  # maximum  value  of log -likelihood  function
  plik.max <- max(-l)
  
  # quantile  of chi -square
  q <- plik.max- 0.5 * qchisq(conf, 1)
  
  # find   the   index   of -q in l for  lower  bound (LB) and  upper  bound (UB)
  plik.lb <- l[l <= -q ][1] 
  ind.lb <- match(plik.lb,l) 
  
  plik.ub <- l[l >= -q & p > mle.p][1]
  ind.ub <- match(plik.ub ,l)
  
  # confidence  interval
  ci.lb <- p[ind.lb] 
  ci.ub <- p[ind.ub] 
  
  # plotting  of the  values
  if(pplot) {
    plot(p,  -l, type = "l", xlab = "p", ylab = 
           "Profile Log-likelihood", ylim = c(-704, plik.max))
    abline(h = q, col = "red")
    abline(h = plik.max, col = 4)
    abline(v = ci.lb, col = "red")
    abline(v = ci.ub, col = "red")
  }
  data.frame("LB" = ci.lb, "UB" = ci.ub,"MLE" = mle.p)
}
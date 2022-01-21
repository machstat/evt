gpd.plikplot <- function(x, theta0, u, eta, z0, p.low, p.up, nint, pplot=0){
  # Create a vector  consisting  of  likelihood  values
  l <- numeric(nint)
  
  # optimizing  algorithm  for  the  parameters
  p <- seq(p.low, p.up, length = nint)
  theta <- theta0[2]
  for (i in 1:nint) {
    p0 <- p[i]
    opt <- optim(theta, fn=gpd.plik, method = "BFGS", x = xdat, u = u, eta = eta, p0 = p0, z0 = z0)
    theta <- opt$par
    
    # negative log -likelihood  for p
    l[i] <- opt$value
  }
  
  # MLE of  probability
  mle.p <- p[which.max(-l)]
  
  # maximum  value  of log -likelihood  function
  plik.max <- max(-l)
  
  # quantile  of chi -square
  q <- plik.max- 0.5 * qchisq(0.95, 1)
  
  # find   the   index   of -q in l for  lower  bound (LB) and  upper  bound (UB)
  plik.LB <- l[l <= -q ][1] 
  ind.LB <- match(plik.LB,l) 
  
  plik.UB <- l[l >= -q & p > mle.p][1]
  ind.UB <- match(plik.UB ,l)
  
  # confidence  interval
  CI.LB <- p[ind.LB] 
  CI.UB <- p[ind.UB] 
  
  # plotting  of the  values
  if(pplot) {
    plot(p,  -l, type = "l", xlab = "p", ylab = 
           "Profile Log-likelihood", ylim = c(-704, plik.max))
    abline(h = plik.max - 0.5 * qchisq(0.95, 1), col = "red")
    abline(h = plik.max, col = 4)
    abline(v = CI.LB, col="red")
    abline(v = CI.UB, col="red")
  }
  data.frame("LB" = CI.LB, "UB" = CI.UB,"MLE"=mle.p)
}

gpd.plik <- function(x, theta, u, eta, p0, z0) {
  xi <- theta
  sgm <- ((u - z0) * theta)*((p0 / eta)^(-theta) - 1)^(-1)
  y <- (u - xdat) / sgm
  y <- 1 + theta * y
  if(any(y <= 0))
    l <- 10^6
  else l <- -(length(xdat) * (log(eta) - log(sgm)) - (1/theta + 1) * sum(log(y)))
  l
}



##### Fitting GPD model
threshold <- c(-2,-3,-4,-5,-6,-8)
u <- threshold[4]
data.gpd <- negDistance[negDistance > u]

model.gpd <- fevd(x = negDistance, threshold = u, type = "GP")
summary.gpd <- summary(model.gpd)
cov.gpd <- summary.gpd$cov.theta
theta.gpd <- summary.gpd$par
ci(model.gpd, type = "parameter")
eta <- model.gpd$rate


### TEST GPD PROFILE LIKELIHOOD
xdat <- c(3.766, 3.766, 2.928, 2.928, 0,0, 3.149, 3.149, 4.653, 4.653, 3.041, 3.041, 3.906, 3.906, 4.839, 4.839, 4.032, 4.032, 3.189, 3.189, 
          3.114, 3.114, 3.380, 3.380, 1.382, 1.382, 3.129, 3.129, 4.840, 4.840, 3.522, 3.522, 2.408, 2.408, 3.120, 3.120, 3.226, 3.226,
          3.918, 3.918, 4.753, 4.753, 1.116, 1.116, 3.561, 3.561, 2.603, 2.603, 3.810, 3.810, 4.478, 4.478, 2.423, 2.423, 4.965, 4.965, 4.274,
          4.274, 3.312, 3.312, 4.450, 4.450, 2.868, 2.868, 1.828, 1.828, 3.503, 3.503, 4.214, 4.214, 4.511, 4.511, 4.512, 4.512, 4.517, 4.517,
          4.322, 4.322, 4.223, 4.223, 4.626, 4.626, 4.386, 4.386, 4.083, 4.083, 3.292, 3.292, 4.263, 4.263, 4.319, 4.319, 3.556, 3.556, 3.859,
          3.859, 3.605, 3.605, 1.407, 1.407, 3.917, 3.917, 4.334, 4.334, 4.236, 4.236, 0.960, 0.960, 3.152, 3.152, 3.419, 3.419, 3.652, 3.652,
          4.751, 4.751, 3.559, 3.559, 2.277, 2.277, 4.641, 4.641, 2.757, 2.757, 2.743, 2.743, 2.597, 2.597, 3.712, 3.712, 2.070, 2.070, 2.469,
          2.469, 3.449, 3.449, 4.866, 4.866, 2.418, 2.418, 2.857, 2.857, 3.915, 2.873, 2.873, 3.915, 3.902, 3.902, 4.803, 4.803, 1.902, 1.902,
          4.158, 4.158, 4.605, 4.605, 2.792, 2.792, 4.149, 4.149, 3.280, 3.280, 2.261, 2.261, 4.237, 4.237, 4.859, 4.859, 1.023, 1.023, 2.772,
          2.772, 4.195, 4.195, 3.006, 3.006, 1.806, 1.806, 1.511, 1.511, 3.244, 3.244, 3.273, 3.273, 4.258, 4.258, 4.103, 4.103, 2.062, 2.062,
          4.789, 4.789, 4.588, 4.588, 2.458, 2.458, 3.374, 3.374, 4.50, 4.505, 2.498, 2.498, 4.302, 4.302, 4.570, 4.570, 3.635, 3.635, 4.511,
          4.511, 4.161, 4.161)

theta <- c(2.0691032, -0.3817382)
u <- 5
z0 <- 0
eta <- 0.1402116
p.low <- 0.000001
p.up <- 0.0013
nint <- 10000 

gpd.plikplot(xdat, theta, u, eta, z0, p.low, p.up, nint, pplot=1)
eta*(1 - pgpd(0+u, scale = theta[1], shape = theta[2])) #0.000173


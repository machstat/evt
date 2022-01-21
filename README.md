# Extreme Value theory
Extreme Value Theory

## gev.pprofMinima: Profile Log-likelihoods for GEV Models for minimas

### Description
Generate profile log-likelihoods wit the function `gev.pplikMinima` for the point estimate (probabilities) of GEV models. Note: does not work for Gumbel distribution.

### Usage
`gev.pprofMinima <- function(x, theta0, z0, p.low, p.up, nint = 1000, conf = 0.95, pplot = 0)`

### Arguments
- x: column vector consisting of data that are fitted to GEV
- theta0: column vector of the parameters of GEV
- z0: the value at which to estimate GEV
- p.low, p.up: the least and greates values at which to compute the profile likelihood
- nint: default value is 1000, it is the number of data points that should be generated between p.low and p.up
- conf: default value is 0.95, it is the confidence coefficient of the confidence interval
- pplot: default value is 0, but takes the values 0 or 1; 1 plots the profile log-likelihood for all the points between `p.low` and `p.up`, 0 does not plot anything

### Value
The confidence interval for the probabilities evaluated at `z0` together with a MLE of `p0` computed from the profile log-likelihood. If `pplot=1`the profile log-likelihood between the specified interval of `p` is plotted. The blue horizontal line corresponds to MLE and the blue corresponds to the value of chi-squared. The intersection between chi-square and the profile likelihood are computed as the lower and upper bound of the confidence interval.

### Examples
```
p.low <- 0.000001
p.up <- 0.053
nint <- 1000
theta0 <- c(-3.562554, 1.278715, -0.264189)
z0 <- 0

gev.pprofMinima(x=x, theta0=theta0, z0, p.low=p.low, p.up=p.up, nint, pplot = 1) # do not run, input data `x` is missing
```
![gevprof](https://github.com/machstat/evt/blob/main/gevprof.png)


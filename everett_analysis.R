
###############
# Perform analysis

library(ggplot2)
source("./everett_analysis_lib.R")
options(warn=1)

d.original <- read.delim("./everett_raw.txt", sep="")

# The lower limit of the pareto distribution can't be larger than
# the upper bound of the first bucket
x.min.upper <- d.original[1, "Max"]

breaks <- d.original[, "Max"][-nrow(d.original)]  # Exclude the last row from the breaks
counts <- d.original[, "Companies.n"]

# The number of different pareto distributions.  If this is 2, that means
# the distribution can have two different power laws.
regions <- 2


UnpackPar <- function(par) {
  #x.min <- LogitFunction(par[1:regions], logit.max=x.min.upper)
  x.min <- cumsum(ExpFunction(par[1:regions]))
  alpha <- ExpFunction(par[(regions + 1):(2 * regions)])
  return(list(x.min=x.min, alpha=alpha))
}


# Par should be the transformed parameters [ x.min, alpha ]
# (the transformed parameters should be able to take any value on the
# real line)
LikForOptim <- function(par) {
  par <- UnpackPar(par)
  x.min <- par$x.min
  alpha <- par$alpha
  #print(c(x.min, alpha))
  l <- CensoredPiecewiseParetoLogLik(alpha=alpha, x.min=x.min, breaks=breaks, counts=counts)
  #print(l)
  return(-l)
}


# Starting parameters
par.0 <- c(InvExpFunction(0.27), InvExpFunction(90),
           InvExpFunction(0.25), InvExpFunction(1.6))

# The number of bootstraps
B <- 1
bootstrap.results <- matrix(nrow=B, ncol=1 + length(par.0))
colnames(bootstrap.results) <- c("bootstrap", "alpha.0", "alpha.1", "x.min.0", "x.min.1")
boot.d <- list()
for (b in 1:B) {
  cat(sprintf("Bootstrap %d\n", b))
  
  # Draw a bootstrap
  if (b == 1) {
    d <- d.original
  } else {
    d <- BootstrapDataFrame(d.original)
  }
  boot.d[[b]] <- d
    
  # Fit the model
  res <- optim(par.0, LikForOptim, method="SANN", control=list(maxit=1e4))
  
  res.par <- UnpackPar(res$par)
  print(res$value)
  print(alpha <- res.par$alpha)
  print(x.min <- res.par$x.min)
  bootstrap.results[b,] <- c(b, res.par$alpha, res.par$x.min)
}
# Uncomment to save automatically when sourcing
#save(bootstrap.results, par.0, d.original, file="bootstrap_results.Rdata")


# Simulate data with the fitted parameters
d <- d.original
sim <- DrawPiecewisePareto(n=sum(d[,"Companies.n"]), x.min=x.min, alpha=alpha)
d$Simulated.Companies.n <- as.numeric(table(cut(sim, breaks=c(0, breaks, Inf))))

# Compare the fitted and raw distributions
ggplot(d, aes(x=log10(Avg))) + 
  geom_line(aes(y=log10(Companies.n), color="Actual"), lwd=2) +
  geom_line(aes(y=log10(Simulated.Companies.n), color="Sim"),lwd=2) +
  geom_vline(aes(xintercept=log10(x.min)), lwd=1)


# Draw a big, uncensored simulated data set for moment computation.
big.sim <- DrawPiecewisePareto(n=1e7, x.min=x.min, alpha=alpha)
sd(big.sim)
mean(big.sim)

risky.utilities <- CRRAUtilty(wealth=wealth, gain=scale * big.sim, gamma=gamma)
risky.utility <- mean(risky.utilities)
#sd(risky.utilities) / sqrt(length(big.sim))

EquivalentWealth <- function(deterministic.gain) {
  return(CRRAUtilty(wealth=wealth + deterministic.gain, gain=0, gamma=gamma) -
                    risky.utility)
}

##################
#
# Set these values and re-execute the uniroot command to 
# evaluate the cash value of different risky exposures.
#
##################

# Calculate the value of exposure to this risk.
wealth <- 3 # This many million in wealth to start
gamma <- 2  # CRRA Coefficient
scale <- 0.01   # The proportion of the company's exit value you actually get

uniroot(EquivalentWealth, interval=c(0, 50))$root

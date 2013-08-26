
# This script follows the notation for Pareto distributions from the
# Wikipedia page.

DrawPareto <- function(n, x.min, alpha) {
  x <- runif(n)
  return(x.min / (x ^ (1 / alpha)))
}

# Draws from a pareto distribution constrained to lie between
# x.min and x.max.
DrawTruncatedPareto <- function(n, x.min, x.max, alpha) {
  x <- runif(n)
  constant <- (x.min ^ (-alpha) - x.max ^ (-alpha))
  return(((x.min ^ (-alpha) - x * constant)) ^ (-1 / alpha))
}


# N draws, each a number from 1:length(probs)
DrawMultinoulli <- function(n, probs) {
  probs <- probs / sum(probs)
  q <- length(probs)  # Distinct number of categories
  
  x <- matrix(rep(runif(n), q), ncol=q)
  x.thresholds <- matrix(rep(cumsum(probs), each=n), ncol=q)

  x <- rowSums(x >= x.thresholds) + 1
  x[x > q] <- q # For floating point errors
  return(x)
}


# The log likelihood that an observation is between
# x.min and x.max (without normalizing constant)
ParetoLogLikInRange <- function(alpha, range.min, range.max) {
  return(log(range.min ^ (-alpha) - range.max ^ (-alpha)))
}


# The log likelihood (with normalizing constant) of
# <counts> observations within buckets defined by
#
# [ x.min, breaks[1] )
# [ breaks[1], breaks[2] )
# …
# [ breaks[len(breaks)], Inf )
#
# …where the lower limit of the Pareto distribution is x.min
# and the shape is alpha.
#
# <counts> should be in the same order as <breaks>, with
# <counts> having one more element than <breaks>.
#
CensoredParetoLogLik <- function(alpha, x.min, breaks, counts) {
  if (!x.min < breaks[1]) {
    stop("x.min must be less than the smallest break.")
  }
  breaks.left <- c(x.min, breaks)
  breaks.right <- c(breaks, Inf)
  log.lik <- breaks.left ^ (-alpha) - breaks.right ^ (-alpha)
  log.lik <- log(log.lik)
  log.lik <- sum(counts * log.lik) +           # Summed observations
             sum(counts) * alpha * log(x.min)  # Normalizing constant
  return(log.lik)
}

# <breaks> should be a sorted numeric vector.
# Returns the index of breaks that is the greatest
# element less than or equal to the corresponding
# element of x.
FirstIndexLessThan <- function(x, breaks) {
  nb <- length(breaks)
  nx <- length(x)
  res <- matrix(rep(x, nb), ncol=nb) >= matrix(rep(breaks, each=nx), ncol=nb)
  res <- rowSums(res)
  res[res == 0] <- NA
  return(res)
}


PiecewiseParetoAreaProbs <- function(alpha, x.min, continuous=TRUE) {
  # An array of the raw areas under the power law in 
  # each region
  
  if (length(alpha) != length(x.min)) {
    stop("alpha and x.min must be the same length")
  }
  
  n <- length(alpha)
  areas <- (x.min ^ (-alpha)) / alpha -
           c((x.min[-1] ^ (-alpha[-n])) / alpha[-n], 0)
  
  # Get an array of constants for each segment.
  if (continuous) {
    # Choose the constants to make the distribution
    # continuous at the breakpoints.  This won't make sense
    # without seeing a write-up.
    
    # The multiples of k0, the constant in the first area.
    k <- c(1, x.min[-1] ^ (diff(alpha)))
    k <- cumprod(k)
    k0 <- 1 / sum(k * areas)    
    constants <- k0 * k
  } else {
    constants <- rep(1/sum(areas), length(areas))
  }
  areas <- areas * constants
  return(list(areas=areas, constants=constants))
}


# The CDF of a pareto with multiple discontinuities of the slope.
# q is a vector of quantiles, alpha a vector of slopes,
# and x.min a vector of the same length indicating where the
# slopes start.  x.min should be sorted in ascending order and
# correspond componentwise to alpha.
pPiecewisePareto <- function(q, alpha, x.min, continuous=TRUE) {
  if (length(alpha) != length(x.min)) {
    stop("Alpha and x.min must have the same length.")
  }
  n <- length(x.min)
  
  pieces <- PiecewiseParetoAreaProbs(alpha, x.min, continuous=continuous)
  areas <- pieces$areas
  constants <- pieces$constants
  
  cumulative.areas <- cumsum(areas)
  cumulative.areas <- c(0, cumulative.areas[-length(cumulative.areas)])

  buckets <- FirstIndexLessThan(q, x.min)
  left.areas <- cumulative.areas[buckets]
  left.boundaries <- x.min[buckets]
  x.alphas <- alpha[buckets]
  x.constants <- constants[buckets]
  
  probs <- left.areas +
           x.constants * (left.boundaries ^ (-x.alphas) - q ^ (-x.alphas)) / x.alphas
  
  return(probs)
}

alpha <- c(1.1, 2.2, 3.3); x.min <- c(1, 100, 1000)
q <- exp(seq(log(0.1), log(10000), length.out=500))
p <- pPiecewisePareto(q, alpha, x.min, continuous=T)

plot(log10(q), c(0, log10(diff(p)))); abline(v=log10(x.min)); abline(h=1)
plot(log10(q), p); abline(v=log10(x.min)); abline(h=1)


CensoredPiecewiseParetoLogLik <- function(alpha, x.min, breaks, counts) {
  if (!x.min[1] < breaks[1]) {
    return(-Inf)
    #stop("x.min must be less than the smallest break.")
  }
  breaks.left <- c(x.min[1], breaks)
  breaks.right <- c(breaks, Inf)
  
  prob.in.range <- (pPiecewisePareto(breaks.right, alpha, x.min) -
                    pPiecewisePareto(breaks.left, alpha, x.min))
  #print(alpha)
  #print(x.min)
  #print(prob.in.range)
  #print("------")
  if (any(counts > 0 & prob.in.range <= 0)) {
    return(-Inf)
  }
  log.lik <- sum(counts * log(prob.in.range))
  return(log.lik)
}


DrawPiecewisePareto <- function(n, alpha, x.min) {
  area.probs <- PiecewiseParetoAreaProbs(alpha, x.min)$areas
  
  # First draw how many observations come from each region  
  x <- DrawMultinoulli(n, area.probs)
  x <- table(x)

  res <- list()
  for (area in names(x)) {
    index <- as.numeric(area)
    if (x[[area]] > 0) {
      x.lower <- x.min[index]
      x.upper <- ifelse(index == length(x.min), Inf, x.min[index + 1])
      res[[area]] <- DrawTruncatedPareto(n=x[[area]], x.min=x.lower, x.max=x.upper,
                                         alpha=alpha[index])
    }
  }
  res <- unlist(res)
  names(res) <- NULL
  return(res)
}
#x <- DrawPiecewisePareto(1000, alpha, x.min)


DrawPiecewiseParetoMoment <- function(n, alpha, x.min, MomentFun) {
  x <- DrawPiecewisePareto(n, alpha, x.min)
  return(mean(MomentFun(x)))
}

# Helper functions for constraining the MLE parameter estimates
LogitFunction <- function(x, logit.min=0, logit.max=1) {
  x <- exp(x) / (exp(x) + 1)
  return(x * (logit.max - logit.min) + logit.min)
}

InvLogitFunction <- function(x, logit.min=0, logit.max=1) {
  x <- (x - logit.min) / (logit.max - logit.min)
  return(log(x / (1 - x)))
}

ExpFunction <- function(x, exp.min=0) {
  return(exp(x) + exp.min)
}

InvExpFunction <- function(x, exp.min=0) {
  return(log(x - exp.min))
}

BootstrapDataFrame <- function(d) {
  count.col <- "Companies.n"
  n <- sum(d[[count.col]])
  p <- d[[count.col]] / n
  d[count.col] <- rmultinom(1, size=n, prob=p)
  return(d)
}


CRRAUtilty <- function(wealth, gain, gamma) {
  return(((wealth + gain) ^ (1 - gamma)) / (1 - gamma))
}

CRRARiskValue <- function(wealth, gain, gamma) {
  return(CRRAUtilty(wealth, gain, gamma) - CRRAUtilty(wealth, 0, gamma))
}










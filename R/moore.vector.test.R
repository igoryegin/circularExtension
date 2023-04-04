###########################################################
#                                                         #
#                                                         #
# MOORE'S MODIFICATION OF RAYLEIGH'S TEST                 #
# FOR VECTOR DATA OR WEIGHTED CIRCULAR OBSERVATIONS       #
# Author: Igor Yegin                                      #
#                                                         #
# Last update: 16/03/2023                                 #
#                                                         #
###########################################################

moore.vector.test <- function(x, w, p.value = c("asymptotic", "simulated")) {
  require(circular)
  if(inherits(x, "3x3"))
    stop("This function does not support 3x3 objects")
  if(length(x) != length(w))
    stop("Vector of angles (x) and vector of weights/lengths (w) must have equal number of elements")
  if(!is.numeric(x) | !is.numeric(w))
    stop("Vector of angles (x) and/or vector of weights/lengths (w) are not numeric")
  INPUT <- deparse(substitute(x))
  p.value <- match.arg(p.value)
  n <- length(x)
  w <- rank(w)
  statistic <- function(x, w) {
    ss <- sum(w * sin(x))
    cc <- sum(w * cos(x))
    (cc^2 + ss^2) / (n * (n + 1) * (2 * n + 1) / 12)
  }
  STATISTIC <- statistic(x = x, w = w)
  method.asymp <- function() {
    assign("METHOD", "Moore's test of circular uniformity for vector data", envir = parent.frame())
    assign("PVAL", 1 - pchisq(STATISTIC, 2), envir = parent.frame())
    assign("PARAMETER", 2, envir = parent.frame())
  }
  method.sim <- function() {
    x.sim <- matrix(rcircularuniform(n * 9999), ncol = 9999)
    w.sim <- replicate(9999, sample(1:n, n))
    sim.statistics <- sapply(seq_len(9999), function(i) statistic(x.sim[, i], w.sim[, i]))
    assign("METHOD", "Moore's test of circular uniformity for vector Data (simulated p-values)", envir = parent.frame())
    assign("PVAL", 1 / (length(sim.statistics) + 1) * (length(which(sim.statistics >= STATISTIC)) + 1), envir = parent.frame())
    assign("PARAMETER", NA, envir = parent.frame())
  }
  if(p.value == "asymptotic") {
    method.asymp()
  }
  else {
    method.sim()
  }
  names(STATISTIC) <- "X2"
  names(PARAMETER) <- "df"
  structure(list(method = METHOD, data.name = INPUT,
                 statistic = STATISTIC, parameter = PARAMETER,
                 p.value = PVAL), class = "htest")
}

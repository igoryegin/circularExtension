###########################################################
#                                                         #
#                                                         #
# MOORE'S MODIFICATION OF RAYLEIGH'S TEST                 #
# FOR VECTOR DATA OR WEIGHTED CIRCULAR OBSERVATIONS       #
# Author: Igor Yegin                                      #
#                                                         #
# Last update: 07/04/2023                                 #
#                                                         #
###########################################################

moore.vector.test <- function(x, w = NULL, p.value = c("asymptotic", "simulated"), rho.test = NULL) {
  require(circular)
  if(inherits(x, "3x3"))
    stop("This function does not support 3x3 objects")
  if(!is.null(w) & length(x) != length(w))
    stop("Vector of angles (x) and vector of weights/lengths (w) must have equal number of elements")
  INPUT <- deparse(substitute(x))
  p.value <- match.arg(p.value)
  if(inherits(x, "polar.coord")) {
    x.t <- x$theta
    n <- length(x$theta)
    w <- x$rho
    w.r <- rank(x$rho)
  }
  else {
    x.t <- x
    n <- length(x)
    w.r <- rank(w)
  }
  statistic <- function(x, w) {
    ss <- sum(w * sin(x))
    cc <- sum(w * cos(x))
    (cc^2 + ss^2) / (n * (n + 1) * (2 * n + 1) / 12)
  }
  STATISTIC <- statistic(x = x.t, w = w.r)
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
  rayleigh <- structure(list(method = METHOD, data.name = INPUT,
                             statistic = STATISTIC, parameter = PARAMETER,
                             p.value = PVAL), class = "htest")
  if(!is.null(rho.test)) {
    rho.success <- which(w > rho.test)
    binom <- binom.test(length(rho.success), length(w), p = 0.5, alternative = "greater")
    binom$data.name <- deparse(substitute(x))
    return(list(rayleigh = rayleigh, binom = binom))
  }
  else {
    return(rayleigh)
  }
}

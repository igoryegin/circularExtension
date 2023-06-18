###########################################################
#                                                         #
#                                                         #
# MOORE'S MODIFICATION OF RAYLEIGH'S TEST                 #
# FOR VECTOR DATA OR WEIGHTED CIRCULAR OBSERVATIONS       #
# Author: Igor Yegin                                      #
#                                                         #
# Last update: 18/06/2023                                 #
#                                                         #
###########################################################

moore.vector.test <- function(x, w = NULL, p.value = c("asymptotic", "simulated"), 
                              rho.CI = c("none", "asymptotic", "bootstrap"), rho.CI.conf = 0.95, 
                              rho.minlength = NULL, rho.minlength.conf = 0.95) {
  require(circular)
  require(boot)
  if(inherits(x, "3x3"))
    stop("This function does not support 3x3 objects")
  if(!is.null(w) & length(x) != length(w))
    stop("Vector of angles (x) and vector of weights/lengths (w) must have equal number of elements")
  INPUT <- deparse(substitute(x))
  p.value <- match.arg(p.value)
  rho.CI <- match.arg(rho.CI)
  if(is.list(x)) {
    X <- x$theta
    n <- length(x$theta)
    w <- x$rho
    W <- rank(x$rho)
  }
  else {
    X <- x
    n <- length(x)
    W <- rank(w)
  }
  statistic <- function(x, w) {
    ss <- sum(w * sin(x))
    cc <- sum(w * cos(x))
    (cc^2 + ss^2) / (n * (n + 1) * (2 * n + 1) / 12)
  }
  boot.fun <- function(data, ind) {
    Rbar <- sqrt(weighted.mean(cos(as.numeric(data$X[ind])), data$W[ind]) ^ 2 +
                   weighted.mean(sin(as.numeric(data$X[ind])), data$W[ind]) ^ 2)
    Rbar
  }
  STATISTIC <- statistic(x = X, w = W)
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
  if(rho.CI == "asymptotic") {
    n <- length(X)
    rhobar <- circular:::TrigonometricMomentRad(X, p = 1, center = FALSE)[2]
    a <- circular:::TrigonometricMomentRad(X, p = 1, center = FALSE)[3]
    a2 <- circular:::TrigonometricMomentRad(X, p = 2, center = FALSE)[3]
    b <- circular:::TrigonometricMomentRad(X, p = 1, center = FALSE)[4]
    b2 <- circular:::TrigonometricMomentRad(X, p = 2, center = FALSE)[4]
    rhocorr <- rhobar - (2 * n + 1) / (3 * n * (n + 1)) * (1 / rhobar ^ 3) * 
      (b^2 / 2 * (1 + a2 - 2 * a^2) - a * b * (b2 - 2*a*b) + a^2 / 2 * (1 - a2 - 2 * b^2))
    rhose <- qnorm(rho.CI.conf) * sqrt(
      (2 * n + 1) / (3 * n * (n + 1)) * (1 / rhobar ^ 2) * 
        (a^2 * (1 + a2 - 2 * a^2) + 2 * a * b * (b2 - 2 * a * b) + b^2 * (1 - a2 - 2 * b^2))
      )
    CINT <- c(pmax(0, rhocorr - rhose), 1)
    attr(CINT, "conf.level") <- rho.CI.conf
  }
  else if(rho.CI == "bootstrap") {
    CINT <- boot(data.frame(X = X, W = W), boot.fun, R = 10000)
    CINT <- c(pmax(0, boot:::basic.ci(CINT$t0, CINT$t, conf = 2 * rho.CI.conf - 1)[4]), 1)
    attr(CINT, "conf.level") <- rho.CI.conf
  }
  else {
    CINT <- NULL
  }
  names(STATISTIC) <- "X2"
  names(PARAMETER) <- "df"
  rayleigh <- structure(list(method = METHOD, data.name = INPUT,
                             statistic = STATISTIC, parameter = PARAMETER,
                             p.value = PVAL, conf.int = CINT), class = "htest")
  if(!is.null(rho.minlength)) {
    success <- which(w > rho.minlength)
    binom <- binom.test(length(success), length(w), p = 0.5, alternative = "greater", conf.level = rho.minlength.conf)
    binom$data.name <- deparse(substitute(x))
    return(list(rayleigh = rayleigh, binom = binom))
  }
  else {
    return(rayleigh)
  }
}

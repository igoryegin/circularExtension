###########################################################
#                                                         #
#                                                         #
# FISHER'S TEST FOR COMMON MEDIAN DIRECTION               #
#                                                         #
#                                                         #
# Author: Igor Yegin                                      #
#                                                         #
# Last update: 07/05/2023                                 #
#                                                         #
###########################################################

median.circular.test <- function(x, p.value = c("asymptotic", "simulated")) {
  p.value <- match.arg(p.value)
  x.unl <- unlist(x)
  nk <- lengths(x)
  N <- length(x.unl)
  theta.tilde <- circular:::MedianCircularRad(x.unl)
  statistic <- function(x) {
    mk <- sapply(x, function(X) length(which(X - theta.tilde < 0)))
    M <- sum(mk)
    N ^ 2 / (M * (N - M)) * sum(mk ^ 2 / nk) - (N * M) / (N - M)
  }
  STATISTIC <- statistic(x = x)
  method.asymp <- function() {
    assign("PVAL", 1 - pchisq(STATISTIC, length(x) - 1), envir = parent.frame())
    assign("METHOD", "Fisher's test for common circular median", envir = parent.frame())
    assign("PARAMETER", length(x) - 1, envir = parent.frame())
  }
  method.sim <- function() {
    vectrunif <- Vectorize(
      function(n)
        runif(n, -pi, pi),
      SIMPLIFY = FALSE
    )
    sim.vectors <- replicate(9999, vectrunif(nk), simplify = FALSE)
    sim.statistics <- sapply(sim.vectors, statistic)
    assign("PVAL", 1 / (length(sim.statistics) + 1) * (length(which(sim.statistics >= STATISTIC)) + 1), envir = parent.frame())
    assign("METHOD", "Fisher's test for common circular median (simulated p-values)", envir = parent.frame())
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
  INPUT <- deparse(substitute(x))
  structure(list(data.name = INPUT, method = METHOD, statistic = STATISTIC,
                 parameter = PARAMETER, p.value = PVAL), class = "htest")
}

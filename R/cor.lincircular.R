###########################################################
#                                                         #
#                                                         #
# TEST OF CIRCULAR-LINEAR CORRELATION (BOTH PARAMETRIC    #
#                                      AND RANK-BASED)    #
#                                                         #
# Author: Igor Yegin                                      #
#                                                         #
# Last update: 28/03/2023                                 #
#                                                         #
###########################################################


cor.lincircular <- function(lvar, cvar, method = c("jwm", "rank"), p.value = c("asymptotic", "simulated")) {
  if(length(lvar) != length(cvar))
    stop("x and y have different lengths")
  n <- length(lvar)
  method <- match.arg(method)
  p.value <- match.arg(p.value)
  if(method == "jwm") {
    statfun <- function(x, y) {
      rxc <- cor(x, cos(y))
      rxs <- cor(x, sin(y))
      rcs <- cor(cos(y), sin(y))
      (rxc ^ 2 + rxs ^ 2 - 2 * rxc * rxs * rcs) / (1 - rcs ^ 2)
    }
    STATISTIC <- statfun(x = lvar, y = cvar)
    method.asymp <- function() {
      assign("PARAMETER", sprintf("F(%d, %d)", 2, n + 3), envir = parent.frame())
      assign("PVAL", 1 - pf((n - 3) * STATISTIC / (2 * (1 - STATISTIC)), 2, n + 3), envir = parent.frame())
      assign("METHOD", "Johnson-Wehrly-Mardia test of linear-circular association", envir = parent.frame())
    }
    method.sim <- function() {
      mc.array <- array(replicate(10000, c(rnorm(n), circular:::RvonmisesRad(n, 0, 1))), dim = c(50, 2, 10000))
      mc.statistics <- apply(mc.array, 1:2, statistic)
      assign("PARAMETER", NA, envir = parent.frame())
      assign("PVAL", 1 / (length(mc.statistics) + 1) * (length(which(mc.statistics >= STATISTIC)) + 1),
             envir = parent.frame())
      assign("METHOD", "Johnson-Wehrly-Mardia test of linear-circular association (simulated p-values)",
             envir = parent.frame())
    }
    if(p.value == "asymptotic") {
      method.asymp()
    }
    else {
      method.sim()
    }
    names(STATISTIC) <- "R^2"
  }
  else {
    statfun <- function(x, y, raw.ranks = NULL) {
      if(!is.null(raw.ranks)) {
        rj <- raw.ranks
      }
      else {
        rj <- rank(y[order(x)])
      }
      Tc <- sum(seq_len(n) * cos(2 * pi * rj / n))
      Ts <- sum(seq_len(n) * sin(2 * pi * rj / n))
      (24 * (Tc ^ 2 + Ts ^ 2)) / (n ^ 3 + n)
    }
    STATISTIC <- statfun(x = lvar, y = cvar, raw.ranks = NULL)
    method.asymp <- function() {
      assign("PARAMETER", "X2(2)", envir = parent.frame())
      assign("PVAL", 1 - pchisq(STATISTIC, 2), envir = parent.frame())
      assign("METHOD", "Mardia's rank test of linear-circular association", envir = parent.frame())
    }
    method.sim <- function() {
      mc.ranks <- replicate(10000, sample(seq_len(n), n))
      mc.statistics <- apply(mc.ranks, 2, statistic, x = NULL, y = NULL)
      assign("PARAMETER", NA, envir = parent.frame())
      assign("PVAL", 1 / (length(mc.statistics) + 1) * (length(which(mc.statistics >= STATISTIC)) + 1),
             envir = parent.frame())
      assign("METHOD", "Mardia's rank test of linear-circular association (simulated p-values)",
             envir = parent.frame())
    }
    if(p.value == "asymptotic") {
      method.asymp()
    }
    else {
      method.sim()
    }
    names(STATISTIC) <- "U"
  }
  names(PARAMETER) <- "null.dist."
  INPUT <- paste0(deparse(substitute(lvar)), " and ", deparse(substitute(cvar)))
  structure(list(data.name = INPUT, method = METHOD, statistic = STATISTIC,
                 parameter = PARAMETER, p.value = PVAL), class = "htest")
}

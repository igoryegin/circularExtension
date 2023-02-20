###########################################################
#                                                         #
#                                                         #
# RAYLEIGH'S TEST OF CIRCULAR UNIFORMITY FOR GROUPED DATA #
#                                                         #
# Author: Igor Yegin                                      #
#                                                         #
# Last update: 20/02/2023                                 #
#                                                         #
###########################################################


grouped.rayleigh.test <- function(x.out, x.zero = NULL, p.value = c("auto", "asymptotic", "simulated"),
                                  template = c("none", "3x3")) {
  if(!is.double(x.out) & !is.integer(x.out))
    stop("non-numeric vector")
  INPUT <- deparse(substitute(x.out))
  p.value <- match.arg(p.value)
  template <- match.arg(template)
  template <- ifelse(template == "3x3" & length(x.out) != 8 & length(x.zero) != 1, "none", template)
  if(template == "3x3") {
    m <- sum(length(x.out), length(x.zero))
    n <- sum(x.out, x.zero)
    x <- c(x.out, x.zero)
    w <- rep(c(1, 0), c(length(x.out), length(x.zero)))
    cosj <- c(cos(2 * pi * seq_len(length(x.out)) / length(x.out)),
              cos(2 * pi * seq_len(length(x.zero)) / length(x.zero)))
    sinj <- c(sin(2 * pi * seq_len(length(x.out)) / length(x.out)),
              sin(2 * pi * seq_len(length(x.zero)) / length(x.zero)))
    coefmat <- w * cosj %*% t(w * cosj) + w * sinj %*% t(w * sinj)
  }
  else {
    m <- length(x.out)
    n <- sum(x.out)
    x <- x.out
    cd <- outer(1:length(x.out), 1:length(x.out), `-`)
    coefmat <- cos(2 * pi * cd / length(x.out))
  }
  statistic <- function(x) {
    as.numeric(
      round(2/length(x.out) * t((x - n/m) / sqrt(n/m)) %*% coefmat %*% (x - n/m) / sqrt(n/m), 5)
    )
  }
  method.asymp <- function() {
    assign("PARAMETER", 2, envir = parent.frame())
    assign("PVAL", 1 - pchisq(STATISTIC, 2), envir = parent.frame())
    assign("METHOD", "Rayleigh Test for Grouped Observations", envir = parent.frame())
  }
  method.mc <- function() {
    assign("PARAMETER", NA, envir = parent.frame())
    mc.vectors <- rmultinom(10000, n, rep(1/m, m))
    mc.statistics <- apply(mc.vectors, 2, statistic)
    assign("PVAL", 1 / length(mc.statistics) * (length(which(mc.statistics >= STATISTIC))), envir = parent.frame())
    assign("METHOD", "Rayleigh Test for Grouped Observations (Monte Carlo p-values)", envir = parent.frame())
  }
  STATISTIC <- statistic(x)
  if(p.value == "asymptotic") {
    method.asymp()
  }
  else if(p.value == "simulated") {
    method.mc()
  }
  else {
    if(any(x < 5)) {
      method.mc()
    }
    else {
      method.asymp()
    }
  }
  names(STATISTIC) <- "X2"
  names(PARAMETER) <- "df"
  structure(list(method = METHOD, data.name = INPUT, statistic = STATISTIC,
                 parameter = PARAMETER, p.value = PVAL), class = "htest")
}

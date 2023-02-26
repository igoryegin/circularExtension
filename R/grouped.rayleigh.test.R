###########################################################
#                                                         #
#                                                         #
# RAYLEIGH'S TEST OF CIRCULAR UNIFORMITY FOR GROUPED DATA #
#                                                         #
# Author: Igor Yegin                                      #
#                                                         #
# Last update: 27/02/2023                                 #
#                                                         #
###########################################################


grouped.rayleigh.test <- function(x.outer, x.zero = NULL, p.value = c("auto", "asymptotic", "simulated"),
                                  template = c("none", "3x3")) {
  if(!is.numeric(x.outer))
    stop("Please provide a vector of whole numbers")
  if(any(x.outer - trunc(x.outer)) > 0)
    warning("Decimal numbers are provided. Only the integer parts of these numbers will be considered")
  x.outer <- trunc(x.outer)
  INPUT <- deparse(substitute(x.outer))
  p.value <- match.arg(p.value)
  template <- match.arg(template)
  template <- ifelse(template == "3x3" & !(length(x.outer) == 8 & length(x.zero) == 1), "none", template)
  if(template == "3x3") {
    m <- sum(length(x.outer), length(x.zero))
    n <- sum(x.outer, x.zero)
    x <- c(x.outer, x.zero)
    w <- rep(c(1, 0), c(length(x.outer), length(x.zero)))
    cc <- c(cos(2 * pi * seq_len(length(x.outer)) / length(x.outer)),
              cos(2 * pi * seq_len(length(x.zero)) / length(x.zero)))
    ss <- c(sin(2 * pi * seq_len(length(x.outer)) / length(x.outer)),
              sin(2 * pi * seq_len(length(x.zero)) / length(x.zero)))
    coefmat <- w * cc %*% t(w * cc) + w * ss %*% t(w * ss)
  }
  else {
    m <- length(x.outer)
    n <- sum(x.outer)
    x <- x.outer
    cd <- outer(1:length(x.outer), 1:length(x.outer), `-`)
    coefmat <- cos(2 * pi * cd / length(x.outer))
  }
  statistic <- function(x) {
    as.numeric(
      round(2/length(x.outer) * t((x - n/m) / sqrt(n/m)) %*% coefmat %*% (x - n/m) / sqrt(n/m), 5)
    )
  }
  method.asymp <- function() {
    assign("PARAMETER", 2, envir = parent.frame())
    assign("PVAL", 1 - pchisq(STATISTIC, 2), envir = parent.frame())
    assign("METHOD", "Rayleigh Test for Grouped Observations", envir = parent.frame())
  }
  method.mc <- function() {
    assign("PARAMETER", NA, envir = parent.frame())
    mc.vectors <- rmultinom(9999, n, rep(1/m, m))
    mc.statistics <- apply(mc.vectors, 2, statistic)
    assign("PVAL", 1 / (length(mc.statistics) + 1) * (length(which(mc.statistics >= STATISTIC)) + 1), envir = parent.frame())
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
  structure(list(method = ifelse(template == "3x3", paste("3x3", METHOD), METHOD),
                 data.name = INPUT, statistic = STATISTIC,
                 parameter = PARAMETER, p.value = PVAL), class = "htest")
}

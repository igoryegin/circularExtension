###########################################################
#                                                         #
#                                                         #
# RAYLEIGH'S TEST OF CIRCULAR UNIFORMITY FOR GROUPED DATA #
#                                                         #
# Author: Igor Yegin                                      #
#                                                         #
# Last update: 21/03/2023                                 #
#                                                         #
###########################################################


grouped.rayleigh.test <- function(x, x.zero = NULL, sym.axes = 1, p.value = c("auto", "asymptotic", "simulated"),
                                  template = c("none", "3x3")) {
  if(!is.numeric(x))
    stop("Please provide a vector of whole numbers")
  if(any(x - trunc(x) > 0))
    warning("Decimal numbers are provided. Only the integer parts of these numbers will be considered")
  x <- trunc(x)
  if(sym.axes - trunc(sym.axes) > 0)
    warning(sprintf("The number of symmetry axis is not an integer. Only the integer part (%d) will be considered", trunc(sym.axes)))
  sym.axes <- trunc(sym.axes)
  INPUT <- deparse(substitute(x))
  p.value <- match.arg(p.value)
  template <- match.arg(template)
  template <- ifelse(template == "3x3" & !(length(x) == 8 & length(x.zero) == 1), "none", template)
  if(template == "3x3") {
    m <- sum(length(x), length(x.zero))
    n <- sum(x, x.zero)
    x <- c(x, x.zero)
    w <- rep(c(1, 0), c(length(x), length(x.zero)))
    cc <- c(cos(2 * sym.axes * pi * seq_len(length(x)) / length(x)),
              cos(2 * sym.axes * pi * seq_len(length(x.zero)) / length(x.zero)))
    ss <- c(sin(2 * sym.axes * pi * seq_len(length(x)) / length(x)),
              sin(2 * sym.axes * pi * seq_len(length(x.zero)) / length(x.zero)))
    coefmat <- w * cc %*% t(w * cc) + w * ss %*% t(w * ss)
  }
  else {
    m <- length(x)
    n <- sum(x)
    x <- x
    cd <- outer(1:length(x), 1:length(x), `-`)
    coefmat <- cos(2 * sym.axes * pi * cd / length(x))
  }
  statistic <- function(x) {
    as.numeric(
      round(2/length(x[-length(x)]) * t((x - n/m) / sqrt(n/m)) %*% coefmat %*% (x - n/m) / sqrt(n/m), 5)
    )
  }
  STATISTIC <- statistic(x = x)
  method.asymp <- function() {
    assign("PARAMETER", 2, envir = parent.frame())
    assign("PVAL", 1 - pchisq(STATISTIC, 2), envir = parent.frame())
    assign("METHOD", "Rayleigh test of circular uniformity for grouped observations", envir = parent.frame())
  }
  method.sim <- function() {
    assign("PARAMETER", NA, envir = parent.frame())
    mc.vectors <- rmultinom(9999, n, rep(1/m, m))
    mc.statistics <- apply(mc.vectors, 2, statistic)
    assign("PVAL", 1 / (length(mc.statistics) + 1) * (length(which(mc.statistics >= STATISTIC)) + 1), envir = parent.frame())
    assign("METHOD", "Rayleigh test of circular uniformity for grouped observations (simulated p-values)", envir = parent.frame())
  }
  if(p.value == "asymptotic") {
    method.asymp()
  }
  else if(p.value == "simulated") {
    method.sim()
  }
  else {
    if(any(x < 5)) {
      method.sim()
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

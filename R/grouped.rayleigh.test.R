###########################################################
#                                                         #
#                                                         #
# RAYLEIGH'S TEST OF CIRCULAR UNIFORMITY FOR GROUPED DATA #
#                                                         #
# Author: Igor Yegin                                      #
#                                                         #
# Last update: 22/03/2023                                 #
#                                                         #
###########################################################


grouped.rayleigh.test <- function(x, sym.axes = 1, p.value = c("auto", "asymptotic", "simulated")) {
  if(!is.numeric(x) & !inherits(x, "3x3"))
    stop("Please provide a vector of whole numbers or a 3x3 object")
  if(!(sym.axes %in% 1:floor(length(x) / 2)))
    stop(sprintf("The number of symmetry axes should be between 1 and %d", floor(length(x) / 2)))
  xunl <- unlist(x)
  sym.axes <- trunc(sym.axes)
  p.value <- match.arg(p.value)
  if(!inherits(x, "3x3")) {
    if(any(x - trunc(x) > 0))
      warning("Decimal numbers are provided. Only the integer parts of these numbers will be considered")
    x <- trunc(x)
    if(sym.axes - trunc(sym.axes) > 0)
      warning(sprintf("The number of symmetry axis is not an integer. Only the integer part (%d) will be considered", trunc(sym.axes)))
    m <- length(x)
    n <- sum(x)
    cd <- outer(1:length(x), 1:length(x), `-`)
    coefmat <- cos(2 * sym.axes * pi * cd / length(x))
  }
  else if(inherits(x, "3x3")) {
    m <- length(xunl)
    n <- sum(xunl)
    w <- rep(c(1, 0), c(length(x$outer), length(x$zero)))
    cc <- c(cos(2 * sym.axes * pi * seq_len(length(x$outer)) / length(x$outer)),
              cos(2 * sym.axes * pi * seq_len(length(x$zero)) / length(x$zero)))
    ss <- c(sin(2 * sym.axes * pi * seq_len(length(x$outer)) / length(x$outer)),
              sin(2 * sym.axes * pi * seq_len(length(x$zero)) / length(x$zero)))
    coefmat <- w * cc %*% t(w * cc) + w * ss %*% t(w * ss)
  }
  statistic <- function(x) {
    if(inherits(x, "3x3"))
      as.numeric(
        round(2/length(x$outer) * t((xunl - n/m) / sqrt(n/m)) %*% coefmat %*% (xunl - n/m) / sqrt(n/m), 5)
      )
    else
      as.numeric(
        round(2/length(x) * t((x - n/m) / sqrt(n/m)) %*% coefmat %*% (x - n/m) / sqrt(n/m), 5)
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
    if(any(xunl < 5)) {
      method.sim()
    }
    else {
      method.asymp()
    }
  }
  INPUT <- deparse(substitute(x))
  names(STATISTIC) <- "X2"
  names(PARAMETER) <- "df"
  structure(list(method = ifelse(inherits(x, "3x3"), paste("3x3", METHOD), METHOD),
                 data.name = INPUT, statistic = STATISTIC,
                 parameter = PARAMETER, p.value = PVAL, test.data = x), class = "htest")
}

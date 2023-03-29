###########################################################
#                                                         #
#                                                         #
# RAYLEIGH'S TEST OF CIRCULAR UNIFORMITY FOR GROUPED DATA #
#                                                         #
# Author: Igor Yegin                                      #
#                                                         #
# Last update: 29/03/2023                                 #
#                                                         #
###########################################################


grouped.rayleigh.test <- function(x, sym.axes = 1, p.value = c("auto", "asymptotic", "simulated"),
                                  boot.CI = NULL) {
  xunl <- unlist(x)
  if(!is.numeric(x) & !inherits(x, "3x3"))
    stop("Please provide a vector of whole numbers or a 3x3 object")
  if(!(sym.axes %in% 1:floor(length(xunl) / 2)))
    stop(sprintf("The number of symmetry axes should be between 1 and %d", floor(length(xunl) / 2)))
  sym.axes <- trunc(sym.axes)
  p.value <- match.arg(p.value)
  if(inherits(x, "3x3")) {
    x.r <- x$outer
    x.b <- x$zero
    m <- length(x.r)
    n <- sum(x.r)
    cd <- outer(1:length(x.r), 1:length(x.r), `-`)
    coefmat <- cos(2 * sym.axes * pi * cd / length(x.r))
  }
  else {
    if(any(x - trunc(x) > 0))
      warning("Decimal numbers are provided. Only the integer parts of these numbers will be considered")
    x <- trunc(x)
    if(sym.axes - trunc(sym.axes) > 0)
      warning(sprintf("The number of symmetry axis is not an integer. Only the integer part (%d) will be considered", trunc(sym.axes)))
    x.r <- x
    m <- length(x.r)
    n <- sum(x.r)
    cd <- outer(1:length(x.r), 1:length(x), `-`)
    coefmat <- cos(2 * sym.axes * pi * cd / length(x.r))
  }
  statistic <- function(x) {
    as.numeric(
        round(2/length(x) * t((x - n/m) / sqrt(n/m)) %*% coefmat %*% (x - n/m) / sqrt(n/m), 5)
      )
  }
  STATISTIC.R <- statistic(x = x.r)
  method.asymp <- function() {
    assign("PARAMETER.R", 2, envir = parent.frame())
    assign("PVAL.R", 1 - pchisq(STATISTIC.R, 2), envir = parent.frame())
    assign("METHOD.R", "Rayleigh test of circular uniformity for grouped observations", envir = parent.frame())
  }
  method.sim <- function() {
    assign("PARAMETER.R", NA, envir = parent.frame())
    mc.vectors <- rmultinom(9999, n, rep(1/m, m))
    mc.statistics <- apply(mc.vectors, 2, statistic)
    assign("PVAL.R", 1 / (length(mc.statistics) + 1) * (length(which(mc.statistics >= STATISTIC.R)) + 1),
           envir = parent.frame())
    assign("METHOD.R", "Rayleigh test of circular uniformity for grouped observations (simulated p-values)",
           envir = parent.frame())
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
  if(!is.null(boot.CI)) {
    CINT.R <- apply(rmultinom(9999, n, x.r),
                  2,
                  function(X) sqrt(statistic(X) / (2 * n)))
    CINT.R <- quantile(CINT.R, c((1 - boot.CI) / 2, 1 - (1 - boot.CI) / 2))
  }
  else {
    CINT.R <- NULL
  }
  INPUT.R <- deparse(substitute(x.r))
  names(STATISTIC.R) <- "X2"
  names(PARAMETER.R) <- "df"
  attr(CINT.R, "conf.level") <- boot.CI
  rayleigh <- structure(list(method = ifelse(inherits(x, "3x3"), paste("3x3", METHOD.R), METHOD.R),
                             data.name = INPUT.R, statistic = STATISTIC.R, conf.int = CINT.R,
                             parameter = PARAMETER.R, p.value = PVAL.R, test.data = x), class = "htest")
  if(inherits(x, "3x3")) {
    binom <- binom.test(x.b, sum(c(x.r, x.b)), p = 0.5, alternative = "less")
    binom$data.name <- deparse(substitute(x))
  }
  return(list(rayleigh = rayleigh, binom = binom))
}

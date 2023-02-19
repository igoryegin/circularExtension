# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hello <- function() {
  print("Hello, world!")
}

grouped.rayleigh.test <- function(x, p.value = c("auto", "asymptotic", "montecarlo"),
                                  template = c("none", "3x3"), neutral = 0) {
  if(!is.double(x) & !is.integer(x))
    stop("non-numeric vector")
  INPUT <- deparse(substitute(x))
  p.value <- match.arg(p.value)
  m <- length(x)
  template <- match.arg(template)
  template <- ifelse(template == "3x3" & m != 8, "none", template)
  m.pol <- ifelse(template == "3x3", m + 1, m)
  if(template == "3x3") {
    x <- c(x, neutral)
    n <- sum(x)
    cd <- outer(1:m, 1:m, `-`)
    cosmat <- cos(2 * pi * cd / m) |> rbind(rep(0, m)) |> cbind(rep(0, m.pol))
  }
  else {
    n <- sum(x)
    p <- 1 / m
    cd <- outer(1:m, 1:m, `-`)
    cosmat <- cos(2 * pi * cd / m)
  }
  statistic <- function(x) {
    round(2/m * t((x - n/m.pol) / sqrt(n/m.pol)) %*% cosmat %*% (x - n/m.pol) / sqrt(n/m.pol), 5) |>
      as.numeric()
  }
  method.asymp <- function() {
    assign("PARAMETER", 2, envir = parent.frame())
    assign("PVAL", 1 - pchisq(STATISTIC, 2), envir = parent.frame())
    assign("METHOD", "Rayleigh Test for Grouped Observations", envir = parent.frame())
  }
  method.mc <- function() {
    assign("PARAMETER", NA, envir = parent.frame())
    mc.vectors <- rmultinom(10000, n, rep(1/m.pol, m.pol))
    mc.statistics <- apply(mc.vectors, 2, statistic)
    assign("PVAL", 1 / length(mc.statistics) * (length(which(mc.statistics >= STATISTIC))), envir = parent.frame())
    assign("METHOD", "Rayleigh Test for Grouped Observations (Monte Carlo p-values)", envir = parent.frame())
  }
  STATISTIC <- statistic(x)
  if(p.value == "asymptotic") {
    method.asymp()
  }
  else if(p.value == "montecarlo") {
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

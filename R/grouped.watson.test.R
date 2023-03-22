###########################################################
#                                                         #
#                                                         #
# WATSON'S TEST OF UNIFORMITY FOR GROUPED CIRCULAR DATA   #
#                                                         #
# Author: Igor Yegin                                      #
#                                                         #
# Last update: 21/03/2023                                 #
#                                                         #
###########################################################

grouped.watson.test <- function(x, p.value = c("asymptotic", "simulated"),
                                template = c("none", "months")) {
  require(mgcv)
  if(inherits(x, "3x3"))
    stop("This function does not support 3x3 objects")
  if(!is.numeric(x))
    stop("Please provide a vector of whole numbers")
  if(any(x - trunc(x) > 0))
    warning("Decimal numbers are provided. Only the integer parts of these numbers will be considered")
  INPUT <- deparse(substitute(x))
  p.value <- match.arg(p.value)
  template <- match.arg(template)
  n <- sum(x)
  m <- length(x)
  if (template == "months" & m == 12) {
    p <- c(31, 28.2425, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31) / 365.2425
  }
  else if (template == "months" & m != 12) {
    warning("Template 'months' is selected, but the length of x is not 12. The template is switched back to 'none'.")
    p <- rep(1 / m, m)
  }
  else {
    p <- rep(1 / m, m)
  }
  statistic <- function(x) {
    Zj <- cumsum(x) - cumsum(n * p)
    Zbar <- sum(Zj * p)
    1 / n * sum(p * (Zj - Zbar) ^ 2)
  }
  X2.weights.fun <- function(i, k) {
    1 / (2 * k ^ 2 * (1 - cos(i * pi / k)))
  }
  if(m %% 2 == 1) {
    X2.weights <- X2.weights.fun(i = seq_len(m - 1) + seq_len(m - 1) %% 2,
                                 k = m)
  }
  else {
    X2.weights <- c(X2.weights.fun(i = seq_len(m - 2) + seq_len(m - 2) %% 2,
                                   k = m),
                    1 / (4 * m ^ 2))
  }
  STATISTIC <- statistic(x = x)
  method.asymp <- function() {
    assign("PVAL", psum.chisq(STATISTIC, X2.weights), envir = parent.frame())
    assign("METHOD", "Watson test of uniformity for grouped circular data", envir = parent.frame())
  }
  method.sim <- function() {
    mc.vectors <- rmultinom(9999, n, p)
    mc.statistics <- apply(mc.vectors, 2, statistic)
    assign("PVAL", 1 / (length(mc.statistics) + 1) * (length(which(mc.statistics >= STATISTIC)) + 1), envir = parent.frame())
    assign("METHOD", "Watson test of uniformity for grouped circular data (simulated p-values)", envir = parent.frame())
  }
  if(p.value == "asymptotic") {
    method.asymp()
  }
  else {
    method.sim()
  }
  names(STATISTIC) <- "U^2"
  structure(list(method = METHOD, data.name = INPUT,
                 statistic = STATISTIC, p.value = PVAL,
                 n.group = m, n.obs = n), class = "htest")
}

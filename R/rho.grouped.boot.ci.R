rho.grouped.boot.ci <- function(x, nrepl = 10000, ci = 0.95) {
  xunl <- unlist(x)
  n <- sum(xunl)
  if(inherits(x, "3x3")) {
    X2 <- sapply(apply(rmultinom(nrepl, n, xunl), 2, as.3x3),
                 function(X) grouped.rayleigh.test(X, p.value = "asymptotic")$statistic)
  }
  else {
    X2 <- apply(rmultinom(nrepl, n, xunl), 2,
                function(X) grouped.rayleigh.test(X, p.value = "asymptotic")$statistic)
  }
  X2 <- sqrt(X2 / (2 * n))
  INTERVAL <- quantile(X2, c((1 - ci) / 2, 1 - (1 - ci) / 2))
  list(x = X2,
       interval = INTERVAL,
       conf.int = ci)
  cat(paste0(100*ci, "% confidence interval: \n(", round(INTERVAL[1], 5), ", ", round(INTERVAL[2], 5), ")"))
}

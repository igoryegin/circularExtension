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


cor.lincircular <- function(lvar, cvar, method = c("jwm", "rank"), p.value = "asymptotic") {
  if(length(lvar) != length(cvar))
    stop("x and y have different lengths")
  n <- length(lvar)
  method <- match.arg(method)
  p.value <- match.arg(p.value)
  if(method == "jwm") {
    rxc <- cor(lvar, cos(cvar))
    rxs <- cor(lvar, sin(cvar))
    rcs <- cor(cos(cvar), sin(cvar))
    STATISTIC <- (rxc ^ 2 + rxs ^ 2 - 2 * rxc * rxs * rcs) / (1 - rcs ^ 2)
    PVAL <- 1 - pf((n - 3) * STATISTIC / (1 - STATISTIC), 2, n - 3)
    names(STATISTIC) <- "R^2"
    METHOD <- "Johnson-Wehrky-Mardia test of linear-circular association"
  }
  else {
    rj <- rank(cvar[order(lvar)])
    Tc <- sum(seq_len(n) * cos(2 * pi * rj / n))
    Ts <- sum(seq_len(n) * sin(2 * pi * rj / n))
    STATISTIC <- (24 * (Tc ^ 2 + Ts ^ 2)) / (n ^ 3 + n)
    PVAL <- 1 - pchisq(STATISTIC, 2)
    names(STATISTIC) <- "U"
    METHOD <- "Mardia's rank test of linear-circular association"
  }
  INPUT <- paste0(deparse(substitute(lvar)), " and ", deparse(substitute(cvar)))
  structure(list(data.name = INPUT, method = METHOD, statistic = STATISTIC, p.value = PVAL), class = "htest")
}

###########################################################
#                                                         #
#                                                         #
# MOORE'S MODIFICATION OF RAYLEIGH'S TEST                 #
# FOR VECTOR DATA OR WEIGHTED CIRCULAR OBSERVATIONS       #
# Author: Igor Yegin                                      #
#                                                         #
# Last update: 27/02/2023                                 #
#                                                         #
###########################################################

moore.test <- function(x, w, p.value = c("asymptotic", "simulated")) {
  if(length(x) != length(w))
    stop("Vector of angles (x) and vector of weights/lengths (w) must have equal number of elements")
  if(!is.numeric(x) | !is.numeric(w))
    stop("Vector of angles (x) and/or vector of weights/lengths (w) are not numeric")
  INPUT <- deparse(substitute(x))
  n <- length(x)
  w <- rank(w)
  ss <- sum(w * sin(x))
  cc <- sum(w * cos(x))
  STATISTIC <- (sqrt(ss^2 + cc^2)) / n ^ (3/2)
  PVAL <- exp(-3 * STATISTIC ^ 2)
  names(STATISTIC) <- "R*"
  structure(list(method = "Moore's Test of Uniformity for Vector Data", data.name = INPUT,
                 statistic = STATISTIC, p.value = PVAL), class = "htest")
}

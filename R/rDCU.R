###########################################################
#                                                         #
#                                                         #
# SAMPLE FROM DISCRETE CIRCULAR UNIFORM DISTRIBUTION      #
#                                                         #
# Author: Igor Yegin                                      #
#                                                         #
# Last update: 14/07/2023                                 #
#                                                         #
###########################################################

rDCU <- function(n, xi = 0, m) {
  dcuniloc <- seq_len(m) * 2 * pi / m + xi
  dcusamp <- circular(sample(dcuniloc, n, replace = TRUE))
  dcusamp
}
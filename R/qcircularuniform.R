###########################################################
#                                                         #
#                                                         #
# CIRCULAR UNIFORM DISTRIBUTION: QUANTILE FUNCTION        #
#                                                         #
# Author: Igor Yegin                                      #
#                                                         #
# Last update: 14/07/2023                                 #
#                                                         #
###########################################################

qcircularuniform <- function(u) {
  ifelse(u < 0, 0, 
         ifelse(u > 1, 2 * pi,
                2 * pi * u)
         )
}

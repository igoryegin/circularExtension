###########################################################
#                                                         #
#                                                         #
# CIRCULAR UNIFORM DISTRIBUTION: DISTRIBUTION FUNCTION    #
#                                                         #
# Author: Igor Yegin                                      #
#                                                         #
# Last update: 21/07/2023                                 #
#                                                         #
###########################################################

pcircularuniform <- function(theta) {
  ifelse(theta < 0, 0, 
         ifelse(theta > 2 * pi, 1,
                theta / (2 * pi))
  )
}

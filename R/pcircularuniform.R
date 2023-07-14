###########################################################
#                                                         #
#                                                         #
# CIRCULAR UNIFORM DISTRIBUTION: DISTRIBUTION FUNCTION    #
#                                                         #
# Author: Igor Yegin                                      #
#                                                         #
# Last update: 14/07/2023                                 #
#                                                         #
###########################################################

pcircularuniform <- function(q) {
  ifelse(q < 0, 0, 
         ifelse(q > 2 * pi, 1,
                q / (2 * pi))
  )
}

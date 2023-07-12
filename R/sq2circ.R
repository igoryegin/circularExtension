###########################################################
#                                                         #
#                                                         #
# SQUARE UNIFORM TO DISC UNIFORM DISTRIBUTION             #
#                                                         #
#                                                         #
# Author: Igor Yegin                                      #
#                                                         #
# Last update: 10/07/2023                                 #
#                                                         #
###########################################################

sq2circ <- function(x, y = NULL) {
  if(is.list(x) & length(x) == 2 & length(x[[1]]) == length(x[[2]])) {
    rho <- pmax(abs(x[[1]]), abs(x[[2]]))
    theta <- 2 * pi * patan2bunif(atan2(y = x[[2]], x = x[[1]])) - pi
  }
  else {
    rho <- pmax(abs(x), abs(y))
    theta <- 2 * pi * patan2bunif(atan2(y = y, x = x)) - pi
  }
  x <- rho * cos(theta)
  y <- rho * sin(theta)
  return(list(x = x, y = y))
}

### Auxiliary functions

datan2bunif <- Vectorize(
  function(x) {
    if((x >= -pi & x < -3 / 4 * pi) | (x >= -pi / 4 & x <= pi / 4) | (x >= 3 / 4 * pi & x <= pi))
      1 / (8 * cos(x) ^ 2)
    else if((x >= -3 / 4 * pi & x <= -pi / 4) | (x >= pi / 4 & x <= 3 / 4 * pi))
      1 / (8 * sin(x) ^ 2)
    else
      0
  }
)

patan2bunif <- Vectorize(
  function(q) {
    if(q <= -pi)
      0
    else if(q > -pi & q <= -3/4 * pi)
      tan(q) / 8
    else if(q > -3/4 * pi & q <= -1/4 * pi)
      1/4 - (1 / tan(q)) / 8
    else if(q > -1/4 * pi & q <= 1/4 * pi)
      1/2 + tan(q) / 8
    else if(q > 1/4 * pi & q <= 3/4 * pi)
      3/4 - (1 / tan(q)) / 8
    else if(q > 3/4 * pi & q <= pi)
      1 + tan(q) / 8
    else
      1
  }
)
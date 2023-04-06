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

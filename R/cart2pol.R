cart2pol <- function(x, y = NULL) {
  if(is.list(x) & length(x) == 2 & length(x[[1]]) == length(x[[2]])) {
    theta <- atan2(y = x[[1]], x = x[[2]])
    rho <- sqrt(x[[1]] ^ 2 + x[[2]] ^ 2)
  }
  else {
    theta <- atan2(y = y, x = x)
    rho <- sqrt(x ^ 2 + y ^ 2)
  }
  structure(list(rho = rho, theta = theta), class = "polar.coord")
}

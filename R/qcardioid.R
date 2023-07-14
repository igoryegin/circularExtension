qcardioid <- Vectorize(
  function(u, mu = circular(0), rho = 0) {
    if (abs(rho) > 0.5) 
      stop("rho must be between -0.5 and 0.5")
    eps <- 10*.Machine$double.eps
    if (u <= eps) {
      theta <- 0
    }
    else if (u >= 1-eps) {
      theta <- 2*pi-eps
    }
    else {
      dcardioidfun <- function(x) (x + 2 * rho * sin(x - mu))/(2 * pi) - u
      theta <- uniroot(dcardioidfun, lower = 0, upper = 2*pi - eps)$root
    }
    theta
  },
  vectorize.args = "u")

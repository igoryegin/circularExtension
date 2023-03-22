as.3x3 <- function(x) {
  structure(list(outer = x[1:8], zero = x[9]), class = "3x3")
}

as.3x3 <- function(x) {
  if(is.numeric(x))
    structure(list(outer = x[1:8], zero = x[9]), class = "3x3")
  else if(any(x - trunc(x) > 0)) {
    warning("Decimal numbers are provided. Only the integer parts of these numbers will be considered")
    x <- trunc(x)
    structure(list(outer = x[1:8], zero = x[9]), class = "3x3")
  }
}

grouped.watson.test <- function(x, p.value = c("auto", "asymptotic", "simulated"),
                                template = c("none", "3x3")) {
  template <- match.arg(template)
  n <- sum(x)
  m <- length(x)
  p <- ifelse(template == "months",
              c(31, 28.2425, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31) / 365.2425,
              n / m)
}

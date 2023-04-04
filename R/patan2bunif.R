qatan2bunif <- Vectorize(
  function(q) {
    ifelse(q >= pi, 1, integrate(datan2bunif, -pi, q)$value)
  }
)

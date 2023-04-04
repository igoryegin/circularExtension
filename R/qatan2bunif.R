qatan2bunif <- Vectorize(
  function(q) {
    integrate(datan2bunif, -pi, q)$value
  }
)

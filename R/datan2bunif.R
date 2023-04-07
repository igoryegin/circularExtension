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

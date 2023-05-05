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

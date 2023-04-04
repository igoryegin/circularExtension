datan2bunif <- function(x) {
  ifelse((x >= -pi & x < -3 / 4 * pi) | (x >= -pi / 4 & x <= pi / 4) | (x >= 3 / 4 * pi & x <= pi),
         1 / (8 * cos(x) ^ 2),
         ifelse((x >= -3 / 4 * pi & x <= -pi / 4) | (x >= pi / 4 & x <= 3 / 4 * pi),
                1 / (8 * sin(x) ^ 2),
                0)
  )
}

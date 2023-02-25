\name{grouped.rayleigh.test}
\alias{grouped.rayleigh.test}
\title{Grouped Rayleigh Test}
\usage{grouped.rayleigh.test()}
\description{
Rayleigh Test of Circular Uniformity for Grouped Data
}
\examples{
library(circular)
x <- rvonmises(100, mu = 0, kappa = 1)
x <- table(findInterval(x, seq(0, 2*pi, pi/6)))
grouped.rayleigh.test(x, p.value = "asymptotic")
}

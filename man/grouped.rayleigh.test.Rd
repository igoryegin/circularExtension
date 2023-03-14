\name{grouped.rayleigh.test}
\alias{grouped.rayleigh.test}
\title{Grouped Rayleigh Test}
\description{
Peforms the Rayleigh's test of circular uniformity for grouped data.
}
\arguments{
  \item{x.outer}{number of observations per group}
  \item{x.zero}{number of observations in the zero group. Only applicable to 3x3 design}
  \item{sym.axes}{number of axes of symmetry}
  \item{p.value}{method of calculating the p-value. If \code{automatic}, it is asymptotic if each group has a count of at least 5. Otherwise, the p-value is simulated with 10,000 repetitions.}
  \item{template}{data design. Currently, only "none" and "3x3" are supported.}
}
\value{A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the chi-square statistic.}
  \item{p.value}{the p-value of the test.}
  \item{parameter}{degrees of freedom. 2 if the p-value is asymptotic, NA if the p-value is simulated.}
  \item{method}{a character string indicating what type of test was
    performed.}
  \item{data.name}{a character string giving the name(s) of the data.}}
\examples{
library(circular)
x <- rvonmises(100, mu = 0, kappa = 0.5)
x <- table(findInterval(x, seq(0, 2*pi, pi/6)))
grouped.rayleigh.test(x, p.value = "asymptotic")
}

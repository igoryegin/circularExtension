\name{grouped.watson.test}
\alias{grouped.watson.test}
\title{Grouped Watson Test}
\description{
  Performs the Watson's test of uniformity for for grouped circular data.
}
\arguments{
  \item{x}{number of observations per group}
  \item{p.value}{method of calculating the p-value.
  Under \eqn{H_0}, the Watson test statistic is distributed as a weighted sum of chi-square distributions, which typically cannot be expressed in a closed form. Thererefore, if \code{asymptotic} method is chosen, the p-value is computed numerically using \code{mgcv::psum.chisq()} function.

  If \code{simulated}, the p-value is simulated with 10,000 repetitions.}
  \item{template}{data design. Currently, only "none" and "months" are supported.}
}
\value{A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the U-squared statistic.}
  \item{p.value}{the p-value of the test.}
  \item{method}{a character string indicating what type of test was
    performed.}
  \item{data.name}{a character string giving the name(s) of the data.}
  }
\examples{
set.seed(1)
library(circular)
x <- rvonmises(100, mu = 0, kappa = 0.5)
x <- table(findInterval(x, seq(0, 2*pi, pi/6)))
grouped.rayleigh.test(x, p.value = "asymptotic")
}
\references{
Moore, B. R. (1980). A Modification of the Rayleigh Test for Vector Data. Biometrika, 67(1), 175–180.
Pewsey, A., Neuhaeuser M., & Ruxton, G. D. (2013). Circular Statistics in R. Oxford: OUP Oxford.
}

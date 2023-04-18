\name{moore.vector.test}
\alias{moore.vector.test}
\title{Moore's Test}
\description{
Performs the Moore's test of circular uniformity for vector data or weighted angular data.
}
\arguments{
  \item{x}{vector of angles}
  \item{w}{vector of lengths or weights}
  \item{p.value}{method of calculating the p-value. \code{asymptotic} calculates the p-value using the formula in Moore (1980). With \code{simulated}, the p-value is simulated with 9,999 repetitions.}
  \item{rho.CI}{confidence level for mean resultant vector length. Confidence interval is computed using normal bootstrap with 10 thousand replications.}
  \item{rho.minlength}{mean resultant length threshold. Sign test (i.e., binomial test with \eqn{H_0:p=0.5}) is performed}
  \item{rho.minlength.conf}{confidence level for \code{rho.minlength}}
}
\value{A list with either one \code{"htest"} object or (if rho.minlength is not \code{NULL}) two such objects, one named \code{rayleigh} and the othe named \code{binom}. Here is a list of some components the ouptup may contain:
  \item{statistic}{the value of the chi-square statistic.}
  \item{conf.int}{confidence interval for the mean resultant vector length.}
  \item{p.value}{the p-value of the test.}
  \item{parameter}{degrees of freedom. 2 if the p-value is asymptotic, NA if the p-value is simulated.}
  \item{estimate}{the estimated probability of a zero observation.}
  \item{null.value}{the probability of success under the null, \code{p.binom}}
  \item{method}{a character string indicating what type of test was
    performed.}
  \item{data.name}{a character string giving the name(s) of the data.}
  }
\examples{
library(circular)
x <- rcircularuniform(100)
w <- runif(100, 0.3, 1)
moore.vector.test(x = x, w = w)

x <- sort(rcircularuniform(100))
w <- sort(runif(100, 0.3, 1))
moore.vector.test(x = x, w = w)
}
\references{
Moore, B. R. (1980). A Modification of the Rayleigh Test for Vector Data. Biometrika, 67(1), 175–180.
}

\name{grouped.rayleigh.test}
\alias{grouped.rayleigh.test}
\title{Grouped Rayleigh Test}
\description{
Performs the Rayleigh's test of circular uniformity for grouped data.
}
\arguments{
  \item{x}{a vector of the number of observations per group, or a '3x3' object}
  \item{sym.axes}{number of axes of symmetry}
  \item{p.value}{method of calculating the p-value. If \code{automatic}, it is asymptotic if each group has a count of at least 5. Otherwise, the p-value is simulated with 9,999 random samples from discrete uniform model.}
  \item{p.binom}{binomal distribution parameter, indicating the null hypothesis probability of a zero observation. Only applicable if a '3x3' object is passed.}
  \item{boot.CI}{bootstrap confidence interval size for mean resultant length (10,000 replications). Parametric bootstrap is used: assuming the counts follow multinomial distribution, its parameters are estimated as proportions of actual observations in each group. This is equivalent to taking the Maximum Likelihood Estimate of multinomial distribution parameters.}
  }
\value{A list with class \code{"htest"} containing some of the following components:
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
set.seed(1)
library(circular)
x <- rvonmises(100, mu = 0, kappa = 0.5)
x <- table(findInterval(x, seq(0, 2*pi, pi/6)))
grouped.rayleigh.test(x, p.value = "asymptotic", boot.CI = 0.9)
}
\references{
Pewsey, A., Neuhaeuser M., & Ruxton, G. D. (2013). Circular Statistics in R. Oxford: OUP Oxford.
}


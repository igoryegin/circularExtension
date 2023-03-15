\name{grouped.watson.test}
\alias{grouped.watson.test}
\title{Grouped Watson Test}
\description{
  Performs the Watson's test of uniformity for for grouped circular data.
}
\arguments{
  \item{x}{number of observations per group}
  \item{p.value}{method of calculating the p-value. If \code{automatic}, it is asymptotic if each group has a count of at least 5. Otherwise, the p-value is simulated with 10,000 repetitions.
  ATTENTION: \code{asymptotic} method requires computing a weighted sum of chi-square distributions, which usually cannot be represented in a closed form. As such, this method has not yet been developed, but we are working on it.}
  \item{template}{data design. Currently, only "none" and "3x3" are supported.}
}
\value{A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the U-squared statistic.}
  \item{p.value}{the p-value of the test.}
  \item{method}{a character string indicating what type of test was
    performed.}
  \item{data.name}{a character string giving the name(s) of the data.}}
\examples{
}

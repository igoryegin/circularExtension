\name{moore.test}
\alias{moore.test}
\title{Moore's Test}
\description{
Performs the Moore's test of circular uniformity for vector data or weighted angular data.
}
\arguments{
  \item{x}{vector of angles}
  \item{w}{vector of lengths or weights}
  \item{p.value}{method of calculating the p-value. \code{asymptotic} calculates the p-value using the formula in Moore (1980). With \code{asymptotic}, the p-value is simulated with 10,000 repetitions.}
}
\value{A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the chi-square statistic.}
  \item{p.value}{the p-value of the test.}
  \item{data.name}{a character string giving the name(s) of the data.}}
\examples{
}

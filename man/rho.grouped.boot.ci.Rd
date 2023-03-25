\name{rho.grouped.boot.ci}
\alias{rho.grouped.boot.ci}
\title{
Bootstrap estimate for mean resultant vector length
}
\description{
Generates a confidence interval for the mean resultant vector length (rho).
Parametric bootstrap is used, assuming the multinomial distribution with the
frequencies estimated by Maximum Likelihood method.
}
\usage{
rho.grouped.boot.ci(x, nrepl = 10000, ci = 0.95)
}
\arguments{
  \item{x}{a vector of whole numbers, or a 3x3 object}
  \item{nrepl}{the number of replications for a bootstrap estimate}
  \item{ci}{confidence interval size}
}
\value{
  \item{x}{A vector of mean resultant vector lengths from nrepl samples}
  \item{interval}{The lower and the upper bound for rho's confidence interval}
  \item{ci}{confidence interval size}
}
\references{
Moore, B. R. (1980). A Modification of the Rayleigh Test for Vector Data. Biometrika, 67(1), 175–180. https://doi.org/10.2307/2335330

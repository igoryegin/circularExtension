\name{Circular Uniform}
\alias{qcircularuniform}
\alias{pcircularuniform}
\title{
  Circular Uniform Distribution
}
\description{
CDF and Quantile Function for the Circular Uniform distribution
}
\usage{
pcircularuniform(theta)
qcircularuniform(u)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{a vector of angles on [0, 2pi)}
  \item{u}{a vector of probabilities}
}
\details{

\code{qcircularuniform} is \eqn{0} when u < 0, and \eqn{2\pi} when u > 1.

}
\value{
  a vector of angles or probabilities
}
\references{
  Pewsey, A., Neuhaeuser M., & Ruxton, G. D. (2013). Circular Statistics in R. Oxford: OUP Oxford.
}
\author{
  Igor Yegin
}

\seealso{
  \code{\link{rcircularuniform}}, \code{\link{dcircularuniform}}
}
\examples{
probs <- pcircularuniform(c(0.1, pi/9, pi/7, pi/6, pi))
probs
angles <- qcircularuniform(c(0.05, 0.08, 0.34, 0.37, 0.4, 0.99))
angles
}

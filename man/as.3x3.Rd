\name{as.3x3}
\alias{as.3x3}
\title{
Turns a vector into an object of class 3x3
}
\description{
}
\usage{
as.3x3(x)
}
\arguments{
  \item{x}{a numeric vector of length 9. If a longer vector is used, only the first 9
  elements will be considered.
}
}
\details{
 The first 8 numbers are assigned to the "outer" element ob an object. The ninth
 number is assigned to the "zero" element.
}
\references{
}
\author{
}
\note{
}

\examples{
x <- c(16, 28, 26, 27, 17, 22, 19, 19, 45)
as.3x3(x)
}

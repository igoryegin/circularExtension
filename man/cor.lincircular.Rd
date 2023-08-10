\name{cor.lincircular}
\alias{cor.lincircular}

\title{
Linear-Circular Correlation
}
\description{
Performs correlation analysis for two vectors, one circular and one linear. 
}
\usage{
cor.lincircular(lvar, cvar, method = c("jwm", "rank"), p.value = c("asymptotic", "simulated"))
}

\arguments{
  \item{lvar}{
  vector with a linear variable
}
  \item{cvar}{
  vector with a circular variable
}
  \item{method}{
  test type. By default, Johnson-Wehrly-Mardia test of circular-linear association is performed. If \code{"rank"} is selected, the function will perform Mardia's rank correlation test.
}
  \item{p.value}{
  method of calculating the p-value. If \code{simulated}, a Monte Carlo p-value with 9,999 simulations is provided.
}
}
\details{
Johnson-Wehrly-Mardia correlation coefficient \eqn{R^2} is always between 0 and 1. For statistical inference, we use the fact that, under \eqn{H_0}, the test statistic \eqn{F = \frac{n - 3}{2}\times{\frac{R^2}{1 - R^2}}} is approximated by F-distribution with 2 and n - 3 degrees of freedom.

Likewise, Mardia's test statistic \eqn{U = \frac{24(T^2_c + T^2_s)}{n^2(n+1)}} has Chi-square distribution with 2 degrees of freedom.
}
\value{

}
\references{
Johnson, R. A., & Wehrly, T. (1977). Measures and Models for Angular Correlation and Angular-Linear Correlation. Journal of the Royal Statistical Society. Series B (Methodological), 39(2), 222–229.

Mardia, K. V., & Jupp, P. E. (2000). Directional Statistics. Wiley-Blackwell.

Pewsey, A., Neuhaeuser M., & Ruxton, G. D. (2013). Circular Statistics in R. Oxford: OUP Oxford.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
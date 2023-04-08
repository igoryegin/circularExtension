\name{Atan2bunif}
\alias{datan2bunif}
\alias{patan2bunif}
\title{
Atan2 transformation of bivariate Uniform(-1, 1) distribution
}
\description{
PDF and CDF of atan2 function applied to bivariate uniform distribution defined as \eqn{f_{xy}(x, y) = \frac{1}{4}} for \eqn{x,y \in [-1, 1]}
}
\usage{
datan2bunif(x)
patan2bunif(q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector of arguments (\code{datan2bunif}) or quantiles (\code{patan2bunif})}
}
\details{
  This is an auxiliary function whose primary purpose is to transform uniform distribtion on a square to uniform distribution on a circle.
}

\value{
a vector of density function values \eqn{f(x)} or distribution function values \eqn{F(x)}
}
\references{
%% ~put references to the literature/web site here ~
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
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.

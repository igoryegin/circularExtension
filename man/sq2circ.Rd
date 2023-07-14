\name{sq2circ}
\alias{sq2circ}
\title{
Mapping of square uniform distribution onto disc uniform distribution.
}
\description{
The points on a square bounded by \eqn{x,y \in [-1, 1]} are mapped onto a circle with zero centre and unit radius. Meanwhile, the uniformity is preserved, i.e if the points are uniformly distributed on a square, they will still have uniform distribution on a circle/disc.
}
\usage{
sq2circ(x, y = NULL)
}
\arguments{
  \item{x}{a vector of x-axis coordinates, or a list of length 2 containing both x-coordinates and y-coordinates
}
  \item{y}{a vector of y-axis coordinates; ignored if x is a list
}
}
\details{
Certain models in social and political sciences (e.g. several kinds of political spectrum models) are represented on a square. It is somewhat inconvenient to apply directional statistics to such data because, even if the data is uniformly spread on a square, certain angles are still more probable. As a matter of fact, under bivariate uniformity, the angles between points and the x-axis follow \link[circularExtension]{Atan2bunif} distribution.

This function serves to tackle this issue by first transforming cartesian coordinates into polar coordinates, and then applying probability integral transform to make angles truly uniform. Finally, the angles are rescales from [0, 1] to [-pi, pi].
}
\value{
a list of cartesian coordinates of transformed data
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

\examples{
  coord.square <- list(x = c(0.98, 0.36, -0.51, 0.03, -0.85, 0.45, 1), y = c(0.88, -0.75, 0.44, 0.06, -1, -0.9, -0.62))
  sq2circ(coord.square)

  coord.x <- c(0.98, 0.36, -0.51, 0.03, -0.85, 0.45, 1)
  coord.y <- c(0.88, -0.75, 0.44, 0.06, -1, -0.9, -0.62))
  sq2circ(coord.x, coord.y)
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

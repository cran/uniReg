\name{equiknots}
\alias{equiknots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Determine the knot sequence.
%%  ~~function to do ... ~~
}
\description{Determines \emph{g+2k+2} knots for the spline basis. The inner knots lie equidistant in \eqn{[a,b]}. If \code{coinc=T}, \emph{k} knots are equal to each \emph{a} and \emph{b}, otherwise the outer knots are also equidistant beyond \eqn{[a,b]}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
equiknots(a, b, g, k, coinc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{Left boundary of the spline interval.
%%     ~~Describe \code{a} here~~
}
  \item{b}{Right boundary of the spline interval.
%%     ~~Describe \code{b} here~~
}
  \item{g}{Number of inner knots.
%%     ~~Describe \code{g} here~~
}
  \item{k}{Spline degree.
%%     ~~Describe \code{k} here~~
}
  \item{coinc}{Should the outer knots be coincident or not? If \code{coinc=TRUE}, there are \emph{k} coincident outer knots at each \emph{a} and \emph{b}.
%%     ~~Describe \code{coinc} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{A numeric vector of length \emph{g+2k+2} with knot locations.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{Claudia Koellmann
%%  ~~who you are~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{unimat}}, \code{\link{unireg}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
equiknots(0,5,3,3,TRUE)
equiknots(0,5,3,3,FALSE)
%
%## The function is currently defined as
%function (a, b, g, k, coinc) 
%{
%    inner <- seq(a, b, length.out = g + 2)
%    if (coinc == T) {
%        outera <- rep(a, k)
%        outerb <- rep(b, k)
%    }
%    else {
%        width <- (b - a)/(g + 1)
%        outera <- a - (k:1 * width)
%        outerb <- b + (1:k * width)
%    }
%    return(c(outera, inner, outerb))
%  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ }
%\keyword{ }% __ONLY ONE__ keyword per line

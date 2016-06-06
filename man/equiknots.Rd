\name{equiknots}
\alias{equiknots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Determine the knot sequence.
%%  ~~function to do ... ~~
}
\description{Determines \emph{g+2k+2} knots for the spline basis of degree \emph{k} on the interval \eqn{[a,b]}. The \emph{g} inner knots lie equidistant in \eqn{[a,b]}. If \code{coinc=TRUE}, the outer knots (\emph{k} on each side of the interval) are placed coincident with \emph{a} and \emph{b}, otherwise the outer knots are also equidistant beyond \eqn{[a,b]}. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
equiknots(a, b, g, k, coinc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{Left numeric boundary of the spline interval.
%%     ~~Describe \code{a} here~~
}
  \item{b}{Right numeric boundary of the spline interval.
%%     ~~Describe \code{b} here~~
}
  \item{g}{A non-negative integer giving the number of inner knots.
%%     ~~Describe \code{g} here~~
}
  \item{k}{A non-negative integer specifying the degree of the spline basis.
%%     ~~Describe \code{k} here~~
}
  \item{coinc}{Logical indicating, if the outer knots should be coincident with the boundary knots or not. If \code{coinc=TRUE}, there are \emph{k} coincident outer knots at \emph{a} as well as at \emph{b}.
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ }
%\keyword{ }% __ONLY ONE__ keyword per line

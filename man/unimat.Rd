\name{unimat}
\alias{unimat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create the matrix of unimodality constraints.
%%  ~~function to do ... ~~
}
\description{Returns a matrix \emph{C} that can be used to specify linear constraints \eqn{C b \ge 0}{C*b >= 0} to impose unimodality with mode at the \emph{m}-th element on a numeric vector \emph{b} of length \emph{p}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
unimat(p, m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{Integer (>=2) giving the length of the vector \emph{b}.
%%     ~~Describe \code{p} here~~
}
  \item{m}{Location of the mode within the vector \emph{b}. Should be in integer between 1 and p.
%%     ~~Describe \code{m} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{Matrix \emph{C} with coefficients for the linear constraints.
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

\seealso{\code{\link{equiknots}}, \code{\link{unireg}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
unimat(4,2)
unimat(5,3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ }
%\keyword{ }% __ONLY ONE__ keyword per line

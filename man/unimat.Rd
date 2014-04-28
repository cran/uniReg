\name{unimat}
\alias{unimat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create the matrix of unimodality constraints.
%%  ~~function to do ... ~~
}
\description{Returns a matrix \emph{A} that can be used to specify linear constraints \eqn{A b \ge 0}{A*b >= 0} to impose unimodality with mode at the \emph{m}-th element on a numeric vector \emph{b} of length \emph{p}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
unimat(p, m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{Length (>=2) of the vector, for which the constraints will be used.
%%     ~~Describe \code{p} here~~
}
  \item{m}{Location of the mode within the vector. Should be \eqn{1 \le m \le p}{ 1 <= m <= p}.
%%     ~~Describe \code{m} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{Matrix with coefficients for the linear constraints.
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
%## The function is currently defined as
%function (p, m) 
%{
%    if (p < 1) {
%        warning("p must be >= 1")
%    }
%    if (p < m) {
%        warning("m must be < = p")
%    }
%    Amon <- function(q) {
%        return(cbind(diag(-1, q - 1), rep(0, q - 1)) + cbind(rep(0, 
%            q - 1), diag(1, q - 1)))
%    }
%    return(rbind(cbind(Amon(m), matrix(0, nrow = m - 1, ncol = p - 
%        m)), cbind(matrix(0, nrow = p - m, ncol = m - 1), -Amon(p - 
%        m + 1))))
%  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ }
%\keyword{ }% __ONLY ONE__ keyword per line

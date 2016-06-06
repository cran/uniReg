\name{print.unireg}
\alias{print.unireg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Print method for \code{unireg} objects.
%%  ~~function to do ... ~~
}
\description{Prints unimodal regression objects.
%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{\method{print}{unireg}(x, ...)
}
\arguments{
  \item{x}{Object of class \code{'unireg'}, a result of function \code{\link{unireg}}.}
  \item{...}{Further arguments (currently not used).}
}
  
\details{
Prints a short overview of a fitted unimodal regression object to the console, namely, the type of the fitted model (including degree of the spline and type of constraint and penalty), the coefficients and their mode location, the tuning parameter and the variance estimate.
}
\value{Invisibly returns the input \code{x}.
}

\author{Claudia Koellmann}
%\note{}

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{unireg}},\code{\link{plot.unireg}},\code{\link{points.unireg}}}

\examples{
x <- sort(rep(0:5,20))         
n <- length(x)               
set.seed(41333)
func <- function(mu){rnorm(1,mu,0.05)}
y <- sapply(dchisq(x,3),func)

# plot of data
plot(jitter(x), y, xlab="x (jittered)")

# fit with default settings
fit <- unireg(x, y, g=5)
# short overview of the fitted spline
fit
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
\keyword{nonparametric}
%\keyword{ }% __ONLY ONE__ keyword per line

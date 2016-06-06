\name{predict.unireg}
\alias{predict.unireg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Predict method for \code{unireg} objects.
%%  ~~function to do ... ~~
}
\description{Predicted values based on a unimodal regression object.
%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{\method{predict}{unireg}(object, newdata, ...)
}
\arguments{
  \item{object}{Object of class \code{'unireg'}, a result of \code{\link{unireg}}.}
  \item{newdata}{A numeric vector of values at which to evaluate the fitted spline function.}
  \item{...}{Further arguments (currently not used).}
}
\details{\code{predict.unireg} produces predicted values by evaluating the fitted regression spline function at the values in \code{newdata}.
}
\value{Returns a numeric vector of predicted function values.
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

# prediction at interim values
predict(fit, c(1.5,2.5,3.5,4.5))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
\keyword{nonparametric}
%\keyword{ }% __ONLY ONE__ keyword per line

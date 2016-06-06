\name{plot.unireg}
\alias{plot.unireg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot method for \code{unireg} objects.
%%  ~~function to do ... ~~
}
\description{Plotting a unimodal regression object.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{\method{plot}{unireg}(x, onlySpline=FALSE, type="l", xlab="x", ylab=NULL, col="black", ...)
}
\arguments{
  \item{x}{Object of class \code{"unireg"}, a result of \code{\link{unireg}}.}
  \item{onlySpline}{Logical indicating whether only the fitted spline or also the original data points should be plotted. Defaults to FALSE (plotting both).}
  \item{type}{Per default plotting type  \code{"l"} is used for the fitted spline.}
  \item{xlab}{Per default the x-axis is labelled with "x".}
  \item{ylab}{If the user does not specify a label for the y-axis, that is when \code{ylab=NULL}(default), prespecified labels like "Fitted unimodal spline function" (depending on the constraint) are used.}
  \item{col}{Colour of the spline function to be plotted (default: black).}
  \item{...}{other parameters to be passed through to the generic \code{plot} function.}
}
\details{This is a plot method for unimodal regression objects. The spline function is plotted using a grid of x-values equally spaced across the interval on which the spline is defined. The distance between the grid values is given by the minimal distance of the observed x-values (used for fitting) divided by 10.}

%\value{}

\author{Claudia Koellmann}
%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{unireg}},\code{\link{predict.unireg}},\code{\link{points.unireg}}}

\examples{
x <- sort(rep(0:5,20))         
n <- length(x)               
set.seed(41333)
func <- function(mu){rnorm(1,mu,0.05)}
y <- sapply(dchisq(x,3),func)

# fit with default settings
fit <- unireg(x, y, g=5)
# short overview of the fitted spline
fit

# plot of fitted spline with and without data
plot(fit, col="orange")
plot(fit, onlySpline=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
\keyword{nonparametric}
%\keyword{ }% __ONLY ONE__ keyword per line

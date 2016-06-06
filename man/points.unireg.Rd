\name{points.unireg}
\alias{points.unireg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Points method for \code{unireg} objects.
%%  ~~function to do ... ~~
}
\description{Plotting a unimodal regression object into an existing plot.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{\method{points}{unireg}(x, type="l", ...)
}
\arguments{
	\item{x}{Object of class \code{"unireg"}, a result of \code{\link{unireg}}.}
  \item{type}{Per default plotting type \code{"l"} is used for the fitted spline.}
  \item{...}{other parameters to be passed through to the generic \code{points} functions.}
}
\details{This is a points method for unimodal regression objects. The spline function is plotted using a grid of x-values equally spaced across the interval on which the spline is defined. The distance between the grid values is given by the minimal distance of the observed x-values (used for fitting) divided by 10.
}

%\value{}

\author{Claudia Koellmann}
%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{unireg}},\code{\link{predict.unireg}},\code{\link{plot.unireg}}}

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

# plot of true and fitted functions
plot(jitter(x), y, xlab="x (jittered)")
curve(dchisq(x,3), 0, 5, type="l", col="grey", lwd=2, add=TRUE)
points(fit, lwd=2, col="orange")
legend("bottomright", legend = c("true mean function", 
        "difference penalized unimodal fit"),
    col=c("grey","orange"),lwd=c(2,2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
\keyword{nonparametric}
%\keyword{ }% __ONLY ONE__ keyword per line

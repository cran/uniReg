\name{unireg}
\alias{unireg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fitting a unimodal penalized spline regression.
%%  ~~function to do ... ~~
}
\description{Function for fitting spline regressions to data. The fit can be constrained to be unimodal, inverse-unimodal, isotonic or antitonic and an arbitrary penalty on the B-spline coefficients can be used.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
unireg(x, y, w=NULL, sigma=NULL, a=min(x), b=max(x), g=10, k=3, 
    constr=c("unimodal","none","invuni","isotonic","antitonic"), 
    penalty=c("diff", "none", "sigEmax", "self", "diag"), Vinv=NULL, 
    beta0=NULL, coinc=NULL, tuning=TRUE, abstol=0.01, vari=5, ordpen=2)
}
\arguments{
  \item{x}{A numeric vector of \emph{x}-values, length \emph{n}. Contains \eqn{p \le n}{p <= n} distinct values.}

  \item{y}{A numeric vector of observed \emph{y}-values of length \emph{n}.}
  
  \item{w}{A numeric vector of positive weights of length \emph{n}, which will be used as the diagonal of the covariance matrix sigma. The weights do not have to sum to \emph{n}, but will be transformed to do so internally. If \code{sigma} is given, \code{w} should be \code{NULL}.}
  
  \item{sigma}{Estimates of the residual (co-)variances. Can be a positive numeric vector of length \emph{p}, giving estimates for the variance at each of the distinct \emph{x}-values. If it is a vector of length 1, a vector of length \emph{p} will be created from it by repetition. It is also possible to provide a whole covariance matrix, if correlated data are present.\\
  If \code{sigma=NULL}, each \emph{x}-value has to be appear at least twice and a global variance (same for all distinct \emph{x}-values) is estimated internally.\\
  If \code{sigma} is given, \code{w} should be \code{NULL}.}

  \item{a}{The left boundary of the interval, on which the spline is defined. If \code{coinc=TRUE}, the spline is zero to the left of this value. By default \code{a} is equal to the minimal \emph{x}-value.}

  \item{b}{The right boundary of the interval, on which the spline is defined. If \code{coinc=TRUE}, the spline is zero to the right of this value. By default \code{b} is equal to the maximal \emph{x}-value.}

  \item{g}{The number of inner knots of the spline, the default is 10.}

  \item{k}{The degree of the spline. Default is 3 (cubic spline).}

  \item{constr}{Shape constraint for the fit. Can be one of "unimodal" (default), "none", "invuni" (inverse-unimodal), "isotonic" or "antitonic".}

  \item{penalty}{Which penalty on the B-spline coefficients should be used? Possible choices are \code{"none"} (no penalty should be used), \code{"diff"} (default) for the differences penalty of order \code{ordpen} and \code{"sigEmax"} for the sigmoid Emax penalty and "self" for a self-defined penalty.  In the latter case, \code{Vinv} and \code{beta0} have to be provided. In the other cases they should be NULL.}
  
  \item{Vinv}{If a self-defined penalty on the B-spline coefficients should be used, \code{Vinv} is the penalty matrix. \code{Vinv} is a \eqn{(g+k+1)\times(g+k+1)}{(g+k+1)x(g+k+1)} matrix. Default: \code{NULL}.}

 \item{beta0}{If a self-defined penalty on the B-spline coefficients should be used, \code{beta0} is the penalty vector of length \emph{g+k+1}. Default: \code{NULL}.}

  \item{coinc}{Should the outer knots of the knot sequence be coincident with the boundary knots or not? Default is \code{NULL} and altering has no effect, if a pre-defined penalty is used. If \code{penalty="self"}, it has to be specified as \code{TRUE} or \code{FALSE}.}

  \item{tuning}{Should the tuning paramater lambda be optimized with or without consideration of the shape constraint? Default: \code{TRUE}. Changing \code{tuning} has no effect, when \code{constr="none"} or \code{penalty="none"}.}

  \item{abstol}{If \code{sigma=NULL}, the iteration between estimation of the residual standard deviation and the B-spline coefficients is stopped, when the absolute change in the variance estimate is less than \code{abstol}. Default: 0.01. Changing \code{abstol} has no effect, when \code{sigma} is provided.}

  \item{vari}{Variance parameter \eqn{sigma_v^2}{sigma_v^2} in the full-rank precision matrix of the prior for beta. By default 5.}

  \item{ordpen}{Order of the difference penalty. Only has an effect, if \code{penalty="diff"}. By default 2.}

}
\details{This function combines implementations of the spline methods described in Koellmann et al. Given paired data \eqn{(x_1,y_1),...,(x_n,y_n)}{(x_1,y_1),...,(x_n,y_n)} it is possible to fit regression splines using the B-spline basis and the maximum likelihood approach. If the spline is unrestricted, the problem reduces to a simple linear regression problem. If the spline is restricted to be unimodal, inverse unimodal, isotonic or antitonic, this leads to a quadratic programming problem. If a penalty is used on the spline coefficients, the tuning parameter is chosen via restricted maximum likelihood (REML).

The data should contain repeated measurements at certain points on the \emph{x}-axis (at least 2 for each point), so that a start estimate of the residual variance can be calculated. Then the function iterates between estimation of the spline coefficients and of the variance. Both estimates will be weighted, if weights are given.
If there is only one measurement (e.g. an already aggregated (mean) measurement), the function expects an input in \code{sigma}, a vector with variances at each \emph{x}-value. The reciprocal variances are used as weights in this case (no further weights can be used).

If no penalty is used, the number of estimable B-spline coefficients, which is \emph{g+k+1}, equals the number of distinct \emph{x}-values. \emph{g} and \emph{k} have to be chosen accordingly.}

\value{A list of
\item{coef}{The vector of estimated B-spline coefficients.}
\item{x}{The original vector of \emph{x}-values.}
\item{fitted.values}{The fitted values at each \emph{x}-value.}
\item{unimod.func}{A function that evaluates the estimated spline at arbitrary \emph{x}-values.}
\item{lambdaopt}{The optimal tuning parameter found via REML.}
\item{sigma}{The estimated residual variance. If the input for \code{sigma} was not \code{NULL}, \code{sigma} equals the input.}
\item{degree}{The degree \code{k} of the spline.}
\item{knotsequence}{The sequence of knots (length \emph{g+2k+2}) used for spline fitting.}
\item{g}{The number \code{g} of inner knots.}
\item{a}{The left boundary of the domain [a,b].}
\item{b}{The right boundary of the domain [a,b].}
\item{variter}{The number of iterations used to estimate the spline coefficients and the variance.}
}

\references{Koellmann, C., Bornkamp, B., and Ickstadt, K., \emph{Unimodal regression using Bernstein-Schoenberg-splines
and penalties}, Technical report 6. URL: http://sfb876.tu-dortmund.de/FORSCHUNG/techreports.html: SFB 876, TU Dortmund, 2012.}

\author{Claudia Koellmann}
%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{unimat}}, \code{\link{equiknots}}}

\examples{
# generate some data
x <- sort(rep(0:5,20))         
n <- length(x)               
set.seed(41333)
func <- function(mu){rnorm(1,mu,0.05)}
y <- sapply(dchisq(x,3),func)

# plot of data
plot(jitter(x),y,xlab="x (jittered)")

# fit with default settings
fit <- unireg(x,y,g=5)

# fit without penalty (we can use at most g=2 inner knots if k=3)
fit2 <- unireg(x,y, penalty="none",g=2)

# fit without penalty and without constraint 
# (does not differ from fit2 with constraint in this case)
fit3 <- unireg(x,y, penalty="none",g=2,constr="none")

# plot of true and fitted functions
z <- seq(0,5,by=0.01)
plot(jitter(x),y,xlab="x (jittered)")#,ylim=c(-0.15,0.35))
points(z,dchisq(z,3),type="l",col="grey",lwd=2)
points(z,fit$unimod.func(z),type="l",lwd=2)
points(z,fit2$unimod.func(z),type="l",col="blue",lwd=2)
points(z,fit3$unimod.func(z),type="l",col="red",lwd=2)
legend("bottomright", legend = c("true mean function", 
        "difference penalized unimodal fit", 
        "unpenalized fit (with and without constraint)"),
    col=c("grey","black","red"),lwd=c(2,2,2))

# estimated variance
fit$sigma
fit2$sigma

\dontrun{
# fit with isotonic, antitonic and inverse-unimodal constraint (just for completeness)
fit4 <- unireg(x,y,constr="antitonic",g=5)
fit5 <- unireg(x,y,constr="isotonic",g=5)
fit6 <- unireg(x,y,constr="invuni",g=5)

points(z,fit4$unimod.func(z),type="l",col="orange",lwd=2)
points(z,fit5$unimod.func(z),type="l",col="brown",lwd=2)
points(z,fit6$unimod.func(z),type="l",col="yellow",lwd=2)

# suppose only aggregated data had been given
means <- c(mean(y[1:20]), mean(y[21:40]), mean(y[41:60]), mean(y[61:80]), 
           mean(y[81:100]), mean(y[101:120]))
sigma <- c(sd(y[1:20]),sd(y[21:40]),sd(y[41:60]),sd(y[61:80]),sd(y[81:100]),sd(y[101:120]))^2

# unimodal fit with differences penalty
fit7 <- unireg(x=unique(x),y=means,g=5,w=NULL,sigma=sigma,abstol=NULL)
plot(unique(x),means,pch=19,ylim=range(y))
points(z,dchisq(z,3),type="l",col="grey",lwd=2)
points(z,fit7$unimod.func(z),type="l",col="green",lwd=2)
legend("bottomright", legend = c("true mean function", "observed mean values", 
    "diff. penalized unimodal fit for means"),
  col=c("grey","black","green"),lty=c(1,NA,1),lwd=c(2,0,2),pch=c(NA,19,NA))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
\keyword{nonparametric}
%\keyword{ }% __ONLY ONE__ keyword per line

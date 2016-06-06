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
unireg(x, y, w=NULL, sigmasq=NULL, a=min(x), b=max(x), g=10, k=3,
	constr=c("unimodal","none","invuni","isotonic","antitonic"),
	penalty=c("diff", "none", "sigEmax", "self", "diag"), Om=NULL,
	beta0=NULL, coinc=NULL, tuning=TRUE, abstol=0.01,vari=5,ordpen=2,
	m=1:(g+k+1), allfits=FALSE, nCores=1)
}
\arguments{
  \item{x}{A numeric vector of \emph{x}-values, length \emph{n}. Contains at least \eqn{d=g+k+1 \le n}{d = g+k+1 <= n} distinct values.}

  \item{y}{A numeric vector of observed \emph{y}-values of length \emph{n}.}
  
  \item{w}{A positive numeric weight vector of length \emph{n}. The weights do not have to sum to \emph{n}, but will be transformed to do so internally. If \code{sigmasq} is given, \code{w} should be \code{NULL} (default).}
    
  \item{sigmasq}{Estimate(s) of the residual (co-)variance(s). Can be a positive numeric vector of length \emph{n}, giving estimates for the variance at each of the \emph{x}-values. If it is a vector of length 1, equal varainces across all \emph{x}-values are assumed.\cr
  If \code{sigmasq=NULL}, each \emph{x}-value has to be appear at least twice and a global variance (same for all \emph{x}-values) is estimated internally.\cr
  If \code{sigmasq} is given, \code{w} should be \code{NULL}.}


  \item{a}{The left numeric boundary of the interval, on which the spline is defined. If \code{coinc=TRUE}, the spline is zero to the left of this value. By default \code{a} is equal to the minimal \emph{x}-value.}

  \item{b}{The right numieric boundary of the interval, on which the spline is defined. If \code{coinc=TRUE}, the spline is zero to the right of this value. By default \code{b} is equal to the maximal \emph{x}-value.}

  \item{g}{A non-negative integer giving the number of inner knots of the spline (default: 10).}

  \item{k}{A non-negative integer specifying the degree of the spline. By default a cubic spline (k = 3) is fitted.}

  \item{constr}{A character string specifying the shape constraint for the fit. Can be one of "unimodal" (default), "none", "invuni" (inverse-unimodal), "isotonic" or "antitonic".}

  \item{penalty}{A character string specifying, which penalty on the B-spline coefficients should be used. Possible choices are \code{"diff"} (default) for the differences penalty of order \code{ordpen}, \code{"none"} for no penalty, \code{"sigEmax"} for the sigmoid Emax penalty, \code{"self"} for a self-defined penalty and \code{"diag"} for a ridge penalty.  For a self-defined penalty, \code{Om} and \code{beta0} have to be provided.}

  \item{Om}{If a self-defined penalty on the B-spline coefficients should be used, \code{Om} is the penalty matrix of dimension \eqn{d \times d}{d x d} and full rank $d$. Otherwise, \code{Om} should be \code{NULL} (default).}

 \item{beta0}{If a self-defined penalty on the B-spline coefficients should be used, \code{beta0} is the penalty vector of length \emph{d}. Otherwise, \code{beta0} should be \code{NULL} (default).}

  \item{coinc}{Logical indicating, if the outer knots of the knot sequence should be coincident with the boundary knots or not? Default is \code{NULL} and altering has no effect, if a pre-defined penalty is used. If \code{penalty="self"}, it has to be specified.}

  \item{tuning}{Logical indicating, if the tuning parameter lambda should be optimized with (\code{tuning=TRUE}, default, computationally expensive) or without (\code{tuning=FALSE}) consideration of the shape constraint. Changing \code{tuning} has no effect, when \code{constr="none"} or \code{penalty="none"}.}

  \item{abstol}{The iterative estimation of the residual variance \eqn{\sigma^2}{sigmasq} and the coefficient vector stops after iteration \eqn{\varsigma}{v}, when \eqn{|\hat{\sigma}^{(\varsigma)} - \hat{\sigma}^{(\varsigma-1)}|}{|hat{sigma}^{(v)} - hat{sigma}^{(v-1)}|} is less than a positive numeric value \code{abstol} (default: 0.01) or when \eqn{\varsigma=10}{v=10}. If \code{sigmasq} is not \code{NULL}, the supplied value is used as starting value in this iteration scheme. There is no iterative estimation, if \code{abstol} is set to \code{NULL}.}

  \item{vari}{Variance parameter \eqn{sigma_v^2 > 0}{sigma_v^2 > 0} in the full-rank precision matrix of the prior for beta. By default 5.}

  \item{ordpen}{Order of the difference penalty (integer \eqn{\ge 0}{ >= 0}, default 2). Only effective, if \code{penalty="diff"}.}
  
  \item{m}{An integer vector specifying the modes of the coefficient vector which should be used for fitting, in explicit, a subset of \emph{\{1,...,d\}}. This argument only has an effect if \code{constr="unimodal"} or \code{"invuni"}.}
  
  \item{allfits}{ Logical indicating if the estimated coefficient vectors for all modes in \code{m} should be returned (\code{TRUE}) or only the one with minimal residual sum of squares (\code{FALSE}).}
  \item{nCores}{The integer number of cores used for parallelization. If \code{nCores=1}, there is no parallelization (default).}
}
	
\details{This function combines implementations of the spline methods described in Koellmann et al. Given paired data \eqn{(x_1,y_1),...,(x_n,y_n)}{(x_1,y_1),...,(x_n,y_n)} it is possible to fit regression splines using the B-spline basis and the maximum likelihood approach. If the spline is unrestricted, the problem reduces to a simple linear regression problem. If the spline is restricted to be unimodal, inverse unimodal, isotonic or antitonic, this leads to a quadratic programming problem. If a penalty is used on the spline coefficients, the tuning parameter is chosen via restricted maximum likelihood (REML).

The data should contain repeated measurements at certain points on the \emph{x}-axis (at least 2 for each point), so that a start estimate of the residual variance can be calculated. Then the function iterates between estimation of the spline coefficients and of the variance. Both estimates will be weighted, if weights are given.
If there is only one measurement per \emph{x}-value, the function expects an input in \code{sigmasq}, an estimate of the variance(s) used for weighted estimation (no further weights can be used).

If no penalty is used, the number of estimable B-spline coefficients, which is \emph{d=g+k+1}, equals the number of distinct \emph{x}-values. \emph{g} and \emph{k} have to be chosen accordingly.}

\value{Returns an object of class "unireg", that is, a list containing the following components: 
\item{x}{The (sorted) vector of \emph{x}-values.}
\item{y}{The input vector of \emph{y}-values (sorted according to \emph{x}).}
\item{w}{The vector of weights used for fitting (sorted according to \emph{x}).}
\item{a}{The left boundary of the domain [a,b].}
\item{b}{The right boundary of the domain [a,b].}
\item{g}{The number \code{g} of inner knots.}
\item{degree}{The degree \code{k} of the spline.}
\item{knotsequence}{The sequence of knots (length \emph{g+2k+2}) used for spline fitting.}
\item{constr}{The constraint on the coefficients.}
\item{penalty}{The type of penalty used.}
\item{Om}{The penalty matrix.}
\item{beta0}{The penalty vector.}
\item{coinc}{The input parameter \code{coinc}.}
\item{tuning}{The input parameter \code{tuning}.}
\item{abstol}{The input value of \code{abstol}.}
\item{vari}{The input variance parameter \code{vari}.}
\item{ordpen}{The order of the difference penalty.}
\item{coef}{The vector of estimated B-spline coefficients (corresponding to the mode with minimal RSS).}
\item{fitted.values}{The fitted values at each \emph{x}-value (corresponding to the mode with minimal RSS).}
\item{lambdaopt}{The optimal tuning parameter found via REML (corresponding to the mode with minimal RSS).}
\item{sigmasq}{The estimated residual variance. If the input for \code{abstol} was \code{NULL}, \code{sigmasq} equals its input value.}
\item{variter}{The number \eqn{\varsigma}{v} of iterations used to estimate the spline coefficients and the variance.}
\item{ed}{The effective degrees of freedom (corresponding to the mode with minimal RSS).}
\item{modes}{The input vector \code{m} of modes.}
\item{allcoefs}{A matrix of coefficient vectors (corresponding to the modes specified in \code{m}) or \code{NULL} (if \code{allfits=FALSE})}
}

\references{Koellmann, C., Bornkamp, B., and Ickstadt, K. (2104), \emph{Unimodal regression using Bernstein-Schoenberg splines and penalties}, Biometrics 70 (4), 783-793.}
\author{Claudia Koellmann}
%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{unimat}}, \code{\link{equiknots}}, \code{\link{plot.unireg}}, \code{\link{points.unireg}}, \code{\link{print.unireg}}, \code{\link{predict.unireg}},}

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

# fit without penalty (we can use at most g=2 inner knots if k=3)
fit2 <- unireg(x, y, penalty="none", g=2)

# plot of fitted spline with or without data
plot(fit2)
plot(fit2, onlySpline=TRUE)

# fit without penalty and without constraint 
# (does not differ from fit2 with constraint in this case)
fit3 <- unireg(x, y, penalty="none", g=2, constr="none")

# plot of true and fitted functions
plot(jitter(x), y, xlab="x (jittered)")
curve(dchisq(x,3), 0, 5, type="l", col="grey", lwd=2, add=TRUE)
points(fit, lwd=2)
points(fit2, col="blue", lwd=2)
points(fit3, col="red", lwd=2)
legend("bottomright", legend = c("true mean function", 
        "difference penalized unimodal fit", 
        "unpenalized fit (with and without constraint)"),
    col=c("grey","black","red"),lwd=c(2,2,2))

# estimated variance
fit$sigmasq
fit2$sigmasq

\dontrun{
# fit with isotonic, antitonic and inverse-unimodal constraint (just for completeness)
fit4 <- unireg(x,y,constr="antitonic",g=5)
fit5 <- unireg(x,y,constr="isotonic",g=5)
fit6 <- unireg(x,y,constr="invuni",g=5)

points(fit4,col="orange",lwd=2)
points(fit5,col="brown",lwd=2)
points(fit6,col="yellow",lwd=2)

# suppose only aggregated data had been given
means <- c(mean(y[1:20]), mean(y[21:40]), mean(y[41:60]), mean(y[61:80]), 
           mean(y[81:100]), mean(y[101:120]))
sigmasq <- c(sd(y[1:20]),sd(y[21:40]),sd(y[41:60]),sd(y[61:80]),sd(y[81:100]),sd(y[101:120]))^2

# unimodal fit with differences penalty
fit7 <- unireg(x=unique(x), y=means, g=5, w=NULL, sigmasq=sigmasq, abstol=NULL)
plot(unique(x), means, pch=19, ylim=range(y))
curve(dchisq(x,3), 0, 5, type="l", col="grey", lwd=2, add=TRUE)
points(fit7, type="l", col="green", lwd=2)
legend("bottomright", legend = c("true mean function", "observed mean values", 
    "diff. penalized unimodal fit for means"),
  col=c("grey","black","green"), lty=c(1,NA,1), lwd=c(2,0,2), pch=c(NA,19,NA))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
\keyword{nonparametric}
%\keyword{ }% __ONLY ONE__ keyword per line

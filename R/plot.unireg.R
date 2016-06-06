plot.unireg <- function(x, onlySpline=FALSE, type="l", xlab="x", ylab=NULL, col="black", ...){
    constr <- x$constr
    if(is.null(ylab)){if(constr=="none"){constr <- "unconstrained"}
                      if(constr=="invuni"){constr <- "inverse unimodal"}
                      ylab <- paste("Fitted", constr, "spline function")
    }
    if(!onlySpline){
        plot(x$x,x$y, xlab=xlab, ylab=ylab, ...)
        points(x, type=type, col=col, ...)
    }else{minx <- min(diff(unique(x$x)))
        z <- seq(x$a, x$b, by=minx/10)
        plot(z, predict.unireg(x,z), type=type, xlab=xlab, ylab=ylab, col=col, ...)
    }
}

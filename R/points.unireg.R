points.unireg <- function(x, type="l", ...){
    minx <- min(diff(unique(x$x)))
    z <- seq(x$a, x$b, by=minx/10)
    points(z, predict.unireg(x,z), type=type, ...)
}

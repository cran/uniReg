print.unireg <- function(x, ...){
    modes <- which(x$coef==max(x$coef))
    constr <- x$constr
    if(constr=="none"){constr <- "unconstrained"}
    if(constr=="invuni"){constr <- "inverse unimodal"
            modes <- which(x$coef==min(x$coef))
    }
    penalty <- x$penalty
    if(penalty=="none"){penalty <- "no penalty"}
    if(penalty=="diff"){penalty <- paste("difference penalty of order", x$ordpen)}
    if(penalty=="sigEmax"){penalty <- "sigmoid Emax penalty"}
    if(penalty=="self"){penalty <- "self-defined penalty"}
    if(penalty=="diag"){penalty <- "ridge penalty"}

    cat("Fitted", constr, "spline of degree", x$degree,"with", penalty, "\n\n")
    cat("Coefficients          ")
    cat(round(x$coef, 2), "\n")
    cat("Mode of coefficients  ")
    cat(modes, "\n")
    cat("Tuning parameter      ")
    cat(round(x$lambdaopt, 2), "\n")
    cat("Variance estimate     ")
    cat(round(x$sigmasq, 2), "\n")
    invisible(x)
}

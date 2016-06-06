equiknots <- function(a,b,g,k,coinc){
###########################################################################################
# Determines g+2k+2 knots for the spline basis. The inner knots lie equidistant in [a,b]. #
# If coinc=T, k knots are equal to each a and b, otherwise the outer knots are also equi- #
# distant beyond [a,b].                                                                   #
###########################################################################################
    if(coinc==T){inner <- seq(a,b,length.out=g+2); outera <- rep(a,k); outerb <- rep(b,k); return(c(outera,inner,outerb))
    }else{width <- (b-a)/(g+1); return(seq(a - k*width, b + k*width, by = width))}
}

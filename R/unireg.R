negloglikFREQ <- function(lambda,tBSB,tySB,sigmaest,tbetaO,tbetaObeta,Dtilde,beta0,Om,rangO,Cm,constr){
################################################################################
## Calculates the value of the negative (restricted) log-likelihood of lambda. #
################################################################################              
            Einv <- tBSB + exp(lambda)*Om
            E <- solve(Einv)
            e <- as.numeric((tySB + exp(lambda)*tbetaO)%*%E) 
            if(constr!="none"){
                prob2 <- pmvnorm(lower=0, mean=as.numeric(Cm%*%e), sigma=Cm%*%E%*%t(Cm))[[1]]
                term2 <- -log(prob2)
                Omtilde <- Dtilde + exp(lambda)*Om
                prob3 <- pmvnorm(lower=0, mean=as.numeric(Cm%*%beta0), sigma=Cm%*%solve(Omtilde)%*%t(Cm))[[1]]
                term3 <- log(prob3)
            }else{term2 <- 0
                term3 <- 0
            }
            term1 <- 0.5*log(det(mean(sigmaest)*Einv)) - (rangO/2)*lambda - 0.5*t(e)%*%Einv%*%e + 0.5*exp(lambda)*tbetaObeta #
            
            negloglik <- term1+term2+term3
            
            if(term3==-Inf){negloglik <- .Machine$double.xmax}
            if(negloglik==Inf){negloglik <- .Machine$double.xmax}
            if(negloglik==-Inf){negloglik <- -.Machine$double.xmax}
    return(negloglik)
}


unisplinem <- function(m,tBSB,tySB,sigmaest,tbetaO,tbetaObeta,Dtilde,rangO,B,beta0,Om,constr,inverse,penalty,tuning){
######################################################################################
# Optimizes the tuning factor lambda and calculates the constrained (with mode m) or #
# unconstrained estimate of the coefficient vector.                                  #
######################################################################################
    d <- dim(B)[2]
    Cm <- inverse*unimat(d,m)          # Cm is the (transposed) matrix of constraints on the beta coefficients if the mode is m
    
    if(penalty=="none"){lambdaopt <- 0
    }else{
        # if tuning=TRUE, lambda ist optimized with constr=constr, otherwise with constr="none"
        if(!tuning){
            opt <- optimize(negloglikFREQ,interval=c(3,10),tBSB=tBSB,tySB=tySB,sigmaest=sigmaest,tbetaO=tbetaO,tbetaObeta=tbetaObeta,Dtilde=Dtilde,beta0=beta0,Om=Om,rangO=rangO,Cm=Cm,constr="none")
            lambdaopt <- exp(opt$minimum)       
        }else{
            opt <- optimize(negloglikFREQ,interval=c(3,10),tBSB=tBSB,tySB=tySB,sigmaest=sigmaest,tbetaO=tbetaO,tbetaObeta=tbetaObeta,Dtilde=Dtilde,beta0=beta0,Om=Om,rangO=rangO,Cm=Cm,constr=constr)
            lambdaopt <- exp(opt$minimum)
        }
    }
        
    # constructing E_{lambda}^{-1}, E_{lambda} and e_{lambda} with the optimal lambda    
    Einv <- tBSB+ lambdaopt*Om
    tySBlBetaO <- as.numeric(tySB + lambdaopt*tbetaO)
    e <- solve(Einv,tySBlBetaO,system=Einv)
    
    # matrices for solve.QP.compact
    Aneu <- unimatind(d,m)
    Amat <- inverse*Aneu$Amat
    Aind <- Aneu$Aind
    
    sc <- norm(Einv,"2")
    # constrained coefficients are estimated with solve.QP.compact; otherwise the vector e is the unconstrained solution
    if(constr!="none"){coeff <- try(solve.QP.compact(Dmat=Einv/sc,dvec=tySBlBetaO/sc,Amat=Amat,Aind=Aind)$solution, silent=TRUE)
            if(class(coeff)=="try-error" || any(is.nan(coeff))){
            coeff <- solve.QP.compact(Dmat=t(solve(chol(Einv/sc))),dvec=tySBlBetaO/sc,Amat=Amat,Aind=Aind, factorized=TRUE)$solution}
        }else{coeff <- e}

    fit <- as.numeric(B%*%coeff)   
    return(list(coef=coeff,fitted.values=fit,lambdaopt=lambdaopt))
}

unireg <- function(x, y, w=NULL, sigmasq=NULL, a=min(x), b=max(x), g=10, k=3, constr=c("unimodal","none","invuni","isotonic","antitonic"), 
    penalty=c("diff", "none", "sigEmax", "self", "diag"), Om=NULL, beta0=NULL, coinc=NULL, tuning=TRUE, abstol=0.01,vari=5,ordpen=2,m=1:(g+k+1),
    allfits=FALSE, nCores=1){
#######################################################
# Applies the specified spline regression to x and y. #
#######################################################
    n <- length(x)
    if(!is.vector(x,mode="numeric")){stop("x should be a numeric vector.")}
    if(!is.vector(y,mode="numeric")){stop("y should be a numeric vector.")}
    if(!is.null(w) && !is.vector(w,mode="numeric")){stop("w should be NULL or a numeric vector.")}
    if(!is.null(sigmasq) && !is.vector(sigmasq,mode="numeric")){stop("sigmasq should be NULL or a numeric vector.")}
    if(!is.null(sigmasq) && length(sigmasq)!=1 && length(sigmasq)!=n){stop("sigmasq a numeric vector of length 1 or same length as x.")}
    if(!(is.vector(a,mode="numeric") && length(a)==1 && is.finite(a))){stop("a should be a finite numeric vector of length 1.")}
    if(!(is.vector(b,mode="numeric") && length(b)==1 && is.finite(b))){stop("b should be a finite numeric vector of length 1.")}
    if(!(g%%1==0 && g>=0 && is.finite(g))){stop("g should be a finite whole number >=0.")}
    if(!(k%%1==0 && k>=0 && is.finite(k))){stop("k should be a finite whole number >=0.")}
    constr <- match.arg(constr)
    penalty <- match.arg(penalty)
    if(!is.null(Om) && !(is.matrix(Om) && isTRUE(all.equal(dim(Om),c(g+k+1,g+k+1))))){stop("Om should be NULL or a (g+k+1)x(g+k+1) matrix.")}
    if(!is.null(beta0) && !is.vector(beta0,mode="numeric")){stop("beta0 should be NULL or a numeric vector.")}
    if(!is.null(coinc) && !is.logical(coinc)){stop("coinc should be NULL, TRUE or FALSE.")}
    if(!is.logical(tuning)){stop("tuning should be TRUE or FALSE.")}
    if(!(is.vector(abstol,mode="numeric") && length(abstol)==1 && is.finite(abstol)) && !is.null(abstol)){stop("abstol should be NULL or a finite numeric vector of length 1.")}
    
    if(!is.null(sigmasq) && !(length(sigmasq)==1 || all(sigmasq==sigmasq[1])) && !is.null(abstol)){stop("Cannot estimate sigmasq in case of heteroscedasticity. Set abstol to NULL.")}
    if(is.null(sigmasq) && is.null(abstol)){stop("Either sigmasq or abstol have to be non-NULL.")}
    if(is.null(sigmasq) && !is.null(abstol) && !all(table(x)>=2)){stop("Provide a startvalue for sigmasq or at least 2 repeated measurements for each x-value.")}
    if(b<=a){stop("[a,b] is not a proper interval.")}
    if(penalty=="none" && (g+k+1)>length(unique(x))){warning("Parameters not estimable. Reduce g+k or increase number ob observation points.")}
    
    if(penalty=="self" && (is.null(Om) || is.null(beta0))){warning("Om and beta0 have to be specified, if penalty='self'.")}
    if(penalty!="self" && (!is.null(Om) || !is.null(beta0))){warning("Om and beta0 will be ignored, if penalty!='self'.")}
    if(penalty!="self" && is.logical(coinc)){warning("coinc has no influence, if penalty!='self'.")}
    if(penalty=="self" && is.null(coinc)){warning("coinc has to be TRUE or FALSE, if penalty='self'.")}

    if(constr=="none" && !tuning){warning("tuning has no influence, if constr='none'.")}
    if(penalty=="none" && !tuning){warning("tuning has no influence, if penalty='none'.")}
        
    if(!is.null(w)){
        if(length(w)!=n){stop("w should be a vector of length n.")
            }else{w <- n*w/sum(w)}
        }else{w <- rep(1,n)}
    
    orderx <- order(x)     
    y <- y[orderx]
    w <- w[orderx]
    x <- sort(x)
    dose <- unique(x)
    nod <- length(dose)
    nd <- as.numeric(table(x))
    d <- g+k+1
    
    constraint <- constr
    inverse <- 1
    if(constr=="invuni"){
        constr <- "unimodal"
        inverse <- -1
    }
    if(constr=="antitonic"){
        constr <- "isotonic"
        inverse <- -1
    }
       
    if(is.null(sigmasq)){
        wssd <- numeric()
        for(u in seq_along(dose)){
            yd <- y[x==dose[u]]
            wd <- w[x==dose[u]]
            wmd <- sum(wd*yd)/sum(wd)
            wssd[u] <- sum(wd*(yd-wmd)^2)
        }
        varest <- sum(wssd)/sum(w)
    }else{varest <- sigmasq}

    yold <- y
    scalfactor <- 0.5*(max(yold) - min(yold))
    shiftfactor <- min(yold) + scalfactor
    y <- (yold-shiftfactor)/scalfactor
    
    if(length(varest)!=1){
        varest <- varest[orderx]
        if(!all(varest==varest[1])){
            w <- 1/varest
            sigmaest <- 1/(scalfactor)^2
        }else{sigmaest <- varest[1]/(scalfactor)^2}
    }else{sigmaest <- varest/(scalfactor)^2}

    if(penalty!="self"){
        if(penalty=="none"){coinc <- TRUE; ordpen=2; Dtilde <- diag(0,d); beta0 <- rep(0,d)
        }else if(penalty=="diff"){coinc <- FALSE; ordpen=ordpen; Dtilde <- diag(1/vari,d); beta0 <- rep(0,d)
        }else if(penalty=="sigEmax"){coinc <- TRUE; ordpen=1; Dtilde <- diag(1/vari,d)
        }else if(penalty=="diag"){coinc <- TRUE; ordpen=0; Dtilde <- diag(0,d); beta0 <- rep(0,d)
        }
        # penalty matrix
        Dq <- diag(d)        
        if(ordpen>=1){Dq <- diff(Dq,difference=ordpen)
            }else{Dtilde <- diag(0,d)} # wie bei "diag"
        Om <- t(Dq)%*%Dq
        
        # determination of knots depending on the interval [a,b], number g of inner knots and degree k of the spline
        # if coinc=T there are k coincident knots at the boundaries a and b
        knotseq <- equiknots(a,b,g,k,coinc) 

        # calculating beta0 for sigEmax
        # we fit a sigmoidEmax model and evaluate it at the knot averages
        if(penalty=="sigEmax"){
            bounds <- defBnds(mD = 8)$sigEmax
            bounds[2,] <- c(1,10)      
            sigE <- fitMod(x, y, model="sigEmax", bnds=bounds)
            knotloc <- knotave(knotseq,d=k)
            beta0 <- predict(sigE, newdata=data.frame(x = knotloc), predType="full-model")
            }
    }else{Dtilde <- diag(0,d)
        beta0 <- (beta0-shiftfactor)/scalfactor
        knotseq <- equiknots(a,b,g,k,coinc)}
    
    # B-spline design matrix
    B <- splineDesign(knots=knotseq, x=x, ord = k+1, outer.ok = T)
        
    tbetaObeta <- t(beta0)%*%Om%*%beta0
    tbetaO <- t(beta0)%*%Om
    rangO <- qr(Om)$rank 
    tBB <- t(B)%*%diag(w)%*%B
    tyB <- t(y)%*%diag(w)%*%B
    
    variter <- 0
    repeat{
        tBSB <- tBB/sigmaest
        tySB <- tyB/sigmaest
        
        if(constr=="unimodal"){    
            vals <- numeric(length(m))
            if(nCores>1){
                cl <- makeCluster(rep("localhost", nCores), type = "SOCK")
                regs <- parLapply(cl=cl, X=m, fun=unisplinem, tBSB=tBSB,tySB=tySB,sigmaest=sigmaest,tbetaO=tbetaO,tbetaObeta=tbetaObeta,Dtilde=Dtilde,rangO=rangO,B=B,beta0=beta0,Om=Om,constr=constr,inverse=inverse,penalty=penalty,tuning=tuning)
                stopCluster(cl)
            }else{
                regs <- lapply(X=m,FUN=unisplinem,tBSB=tBSB,tySB=tySB,sigmaest=sigmaest,tbetaO=tbetaO,tbetaObeta=tbetaObeta,Dtilde=Dtilde,rangO=rangO,B=B,beta0=beta0,Om=Om,constr=constr,inverse=inverse,penalty=penalty,tuning=tuning)
                }
            for (i in seq_along(m)){vals[i] <- sum(w*(y - regs[[i]]$fitted.values)^2)}   
            mini <- which.min(vals)    
            coef.unimod <- regs[[mini]]$coef
            fit.unimod <- regs[[mini]]$fitted.values
            lambdaopt <- regs[[mini]]$lambdaopt
        }else{ # constr = isotonic or none
            erg <- unisplinem(m=d,tBSB=tBSB,tySB=tySB,sigmaest=sigmaest,tbetaO=tbetaO,tbetaObeta=tbetaObeta,Dtilde=Dtilde,rangO=rangO,B=B,beta0=beta0,Om=Om,constr=constr,inverse=inverse,penalty=penalty,tuning=tuning)
            coef.unimod <- erg$coef
            fit.unimod <- erg$fitted.values
            lambdaopt <- erg$lambdaopt
        }
        Hperm <- tBSB%*%solve(tBSB+ lambdaopt*Om)
        ed <- sum(diag(Hperm))
        
        if(!is.null(sigmasq) && is.null(abstol)){break}
        
        resd <- y-fit.unimod
        varest <- sum(w*(resd-sum(resd*w)/sum(w))^2)/sum(w)
        variter <- variter+1
        if(abs(varest - sigmaest)<abstol || variter>=10){
            sigmaest <- varest
            break}
        sigmaest <- varest
    }
    if(allfits){
        allcoefs <- matrix(0,nrow=0,ncol=d)
        for(i in seq_along(m)){
            allcoefs <- rbind(allcoefs, scalfactor*regs[[i]]$coef + shiftfactor)
            }
    }else{
        allcoefs <- NULL
    }
    
    res <- list(x=x,y=yold,w=w,a=a,b=b,g=g,degree=k,knotsequence=knotseq,constr=constraint,penalty=penalty,Om=Om,beta0=beta0,coinc=coinc,tuning=tuning,abstol=abstol,vari=vari,ordpen=ordpen,
            coef=scalfactor*coef.unimod+shiftfactor,fitted.values=scalfactor*fit.unimod+shiftfactor,lambdaopt=lambdaopt,sigmasq=sigmaest*scalfactor^2,variter=variter,ed=ed,modes=m,allcoefs=allcoefs)
    class(res) <- "unireg" 
    res
}

unimat <- function(p,m){
#######################################################################################
# Generates the matrix A of unimodality-constraints with mode m (for p coefficients). #
#######################################################################################
    if(p<1){warning("p must be >= 1")}
    if(p<m){warning("m must be < = p")}
    Amon <- function(q){return(cbind(diag(-1,q-1),rep(0,q-1)) + cbind(rep(0,q-1),diag(1,q-1)))}
    return(rbind(cbind(Amon(m),matrix(0,nrow=m-1,ncol=p-m)),
                 cbind(matrix(0,nrow=p-m,ncol=m-1),-Amon(p-m+1))))
}


unimatind <- function(p,m){
    Amat <- rbind(c(rep(-1,m-1),rep(1,p-m)), c(rep(1,m-1),rep(-1,p-m)))
    Aind <- rbind(rep(2,p-1), 1:(p-1), 2:p)
    return(list(Amat=Amat,Aind=Aind))
}

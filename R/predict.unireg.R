predict.unireg <- function(object, newdata, ...){
      Bneu <- splineDesign(knots=object$knotsequence, newdata, ord = object$degree+1, outer.ok = T)
      as.numeric(Bneu%*%object$coef)
}

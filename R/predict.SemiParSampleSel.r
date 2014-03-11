predict.SemiParSampleSel <- function(object, eq, ...){

   if(eq==1){ ss.pred <- object$gam1; ind <- 1:length(coef(ss.pred)) }
   else     { ss.pred <- object$gam2; ind <- (length(coef(object$gam1))+1):(length(coef(object$gam1)) + length(coef(ss.pred))) }
   
   ss.pred$coefficients <- object$coefficients[ind]
   ss.pred$Vp <- object$Vb[ind,ind]
   
   predict.gam(ss.pred, ...)
   

}


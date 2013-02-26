plot.SemiParSampleSel <- function(x, eq, ...){
   
    if(eq==1){ 
          ss.plot <- x$gam1
          ind <- 1:length(coef(ss.plot)) 
    } 
    else{ 
          ss.plot <- x$gam2
          fir <- length(coef(x$gam1)) + 1 
          sec <- fir - 1 + length(coef(ss.plot))
          ind <- fir:sec 
    }
                                      
    ss.plot$coefficients <- x$coefficients[ind]
    ss.plot$Vp  <- x$Vb[ind,ind]
    ss.plot$edf <- diag(x$F)[ind]
   
    plot.gam(ss.plot, ...)
          
}
   
   
   
   
   
   
   
   
st.theta.star <- function(start.theta, co, BivD) {


  if(BivD=="N") { 
        if(is.null(start.theta)) {
             co <- ifelse( abs(co) > 0.9, sign(co)*0.9, co); 
             a.theta <- atanh(co)
        }
        else { if(abs(start.theta)<=1) a.theta <- atanh(start.theta) else stop("Wrong initial value of theta/rho.") }
  }
  

  if(BivD=="C" | BivD=="rC") {
        if(is.null(start.theta)) a.theta <- log(3) 
        else { if(start.theta>0) a.theta <- log(start.theta) else stop("Wrong initial value of theta.") } 
  } 
 
  
  if(BivD=="J" | BivD=="rJ") {
        if(is.null(start.theta)) a.theta <- log(3)
        else { if(start.theta>1) a.theta <- log(start.theta-1) else stop("Wrong initial value of theta.") }
  }  
  
  
  if(BivD=="F") {
        if(is.null(start.theta)) { a.theta <- iRho(frankCopula(),co); a.theta <- ifelse(abs(a.theta)<15,a.theta,sign(a.theta)*15) }
        else a.theta <- start.theta
  }  
  

  if(BivD=="FGM") {
        if(is.null(start.theta)) {
              co <- ifelse( abs(co) < 0.233, co, sign(co)*0.233 )   # so that abs(a.theta)<0.7
              a.theta <- iRho(fgmCopula(),co); a.theta <- atanh(a.theta)
        }
        else { if(abs(start.theta)<=1) a.theta <- atanh(start.theta) else stop("Wrong initial value of theta.") }
  } 
  

  if(BivD=="AMH") {
        if(is.null(start.theta)) {
              co <- ifelse( co > -0.18, co, -0.18 ); co.m <- ifelse( co < 0.33, co, 0.33 )   
              a.theta <- iTau(amhCopula(),co); a.theta <- ifelse( abs(a.theta) < 0.9, a.theta, sign(a.theta)*0.9 ); a.theta <- atanh(a.theta)
        }
        else { if(abs(start.theta)<=1) a.theta <- atanh(start.theta) else stop("Wrong initial value of theta.") }
  }  
  

  if(BivD=="G" | BivD=="rG") {
        if(is.null(start.theta)) a.theta <- log(3)
        else { if(start.theta>1) a.theta <- log(start.theta-1) else stop("Wrong initial value of theta.") }
  }  



return(a.theta)


}















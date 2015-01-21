working.comp <- function(x, VC=VC, myf=myf){

  der.1 <- x$dl.dbe1
  der.2 <- x$dl.dbe2
  e.par <- x$argument
  
  ll <- length(e.par) 

  if(VC$margins[2] %in% c("N","G","NB","PIG")) expa <- c((ll-1) : ll)
  if(VC$margins[2] == "P")                     expa <- c(ll)
  if(VC$margins[2] %in% c("D","S"))            expa <- c((ll-2) : ll)
   
  e.par <- e.par[-expa]
  

  sqq <- seq(1,(2*VC$n-1),by=2)

  X <- rW.X <- Matrix(0,2*VC$n,(VC$X1.d2+VC$X2.d2))
  D <- rW.Z <- Matrix(0,2*VC$n,1)

  X[sqq,   1:VC$X1.d2]                      <- VC$X1
  X[sqq+1,(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)] <- VC$X2

  D[sqq,1]   <- x$dl.dbe1
  D[sqq+1,1] <- x$dl.dbe2

  be1be2 <- cbind(x$d2l.be1.be1, x$d2l.be1.be2, x$d2l.be2.be2)

  L.W <- unlist( apply(be1be2, 1, myf ) , recursive = FALSE) 

  W.inv <- c.W <- list()

    for(i in 1:VC$n){
    
      W.eig <- eigen(L.W[[i]], symmetric=TRUE)  
     
      s.ch <- sum(as.numeric(as.vector(L.W[[i]])==0))
   
      if(s.ch==3){ c.W[[i]] <- W.inv[[i]] <- matrix(0,2,2)
                   c.W[[i]][1,1] <- sqrt(L.W[[i]][1,1]) 
                   W.inv[[i]][1,1] <- 1/L.W[[i]][1,1]
                 }
      if(s.ch!=3){ if(min(W.eig$values) < sqrt(.Machine$double.eps)) { pep <- which(W.eig$values < sqrt(.Machine$double.eps))
      	                                                               W.eig$values[pep] <- 0.0000001}     
                                                                                                              
                 }      
      
 if( s.ch!=3 ){
      c.W[[i]]   <- W.eig$vec%*%tcrossprod(diag(sqrt(W.eig$val)),W.eig$vec) 
      W.inv[[i]] <- W.eig$vec%*%tcrossprod(diag(1/W.eig$val    ),W.eig$vec)
              }      

                 } # end loop
                  

      c.W   <- bdiag(c.W)
      W.inv <- bdiag(W.inv)

      rW.X <- as.matrix(c.W%*%X) 
      rW.Z <- as.matrix(c.W%*%( X%*%e.par + W.inv%*%D )) 
           
rm(c.W, X, W.inv, D)
gc()

   
 list( rW.X=rW.X , rW.Z=rW.Z )

}


      
                           























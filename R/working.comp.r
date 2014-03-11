working.comp <- function(x,X1=X1,X2=X2,X1.d2=X1.d2,X2.d2=X2.d2,n=n,myf=myf){

  der.1 <- x$dl.dbe1
  der.2 <- x$dl.dbe2
  e.par <- x$argument

  sqq <- seq(1,(2*n-1),by=2)

  X <- rW.X <- Matrix(0,2*n,(X1.d2+X2.d2))
  D <- rW.Z <- Matrix(0,2*n,1)

  X[sqq,   1:X1.d2]                <- X1
  X[sqq+1,(X1.d2+1):(X1.d2+X2.d2)] <- X2

  D[sqq,1]   <- x$dl.dbe1
  D[sqq+1,1] <- x$dl.dbe2

  be1be2 <- cbind(x$d2l.be1.be1, x$d2l.be1.be2, x$d2l.be2.be2)

  L.W <- unlist( apply(be1be2, 1, myf ) , recursive = FALSE) 

  W.inv <- c.W <- list()

    for(i in 1:n) {

      W.eig <- eigen(L.W[[i]], symmetric=TRUE)  

      if(min(W.eig$values) < .Machine$double.eps){ L.W[[i]] <- nearPD( L.W[[i]], ensureSymmetry = FALSE )$mat
                                                   W.eig <- eigen(L.W[[i]], symmetric=TRUE) 
                                                 }
      
      c.W[[i]]   <- W.eig$vec%*%tcrossprod(diag(sqrt(W.eig$val)),W.eig$vec) 
      W.inv[[i]] <- W.eig$vec%*%tcrossprod(diag(1/W.eig$val    ),W.eig$vec)       

                  }

      c.W <- bdiag(c.W)
      W.inv <- bdiag(W.inv)

      rW.X <- as.matrix(c.W%*%X) 
      rW.Z <- as.matrix(c.W%*%( X%*%e.par + W.inv%*%D )) 
              
 list( rW.X=rW.X , rW.Z=rW.Z )

}



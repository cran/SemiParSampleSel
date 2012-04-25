working.comp <- function(x,X1=X1,X2=X2,X1.d2=X1.d2,X2.d2=X2.d2,n=n){
  #x=fit
  der.1 <- x$dl.dbe1
  der.2 <- x$dl.dbe2
  e.par <- x$argument
  W <- c.W <- i.W <- matrix(0,2,2)
  X <- rW.X <- matrix(0,2*n,(X1.d2+X2.d2))
  D <- iW.D <- Z <- rW.Z <- matrix(0,2*n,1)
  j <- 1
    for(i in seq(1,(2*n-2),by=2)) {
      X[i,1:X1.d2]                   <- X1[j,]
      X[i+1,(X1.d2+1):(X1.d2+X2.d2)] <- X2[j,]
      D[i,1]   <- der.1[j]
      D[i+1,1] <- der.2[j]
      W <- matrix(c( x$d2l.be1.be1[j],x$d2l.be1.be2[j],     
                     x$d2l.be1.be2[j],x$d2l.be2.be2[j] ) , 2 , 2 ) 
      W.eig <- eigen(W)
      c.W   <- W.eig$vec%*%diag(sqrt(pmax(W.eig$val,sqrt(.Machine$double.eps))))%*%t(W.eig$vec) 
      W.inv <- W.eig$vec%*%diag(1/pmax(W.eig$val,sqrt(.Machine$double.eps)))%*%t(W.eig$vec) 

      rW.X[i:(i+1),1:(X1.d2+X2.d2)] <- c.W%*%X[i:(i+1),1:(X1.d2+X2.d2)]
      rW.Z[i:(i+1),1] <- c.W%*%( X[i:(i+1),]%*%e.par + W.inv%*%D[i:(i+1),1] )
      j <- j + 1
    }
 list( rW.X=rW.X , rW.Z=rW.Z )
}



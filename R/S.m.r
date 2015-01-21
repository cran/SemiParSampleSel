S.m <- function(gam1,gam2,l.sp1,l.sp2){

  Ss <- list()
  
  off <- rank <- pc1 <- pc2 <- 0 
  
  i1 <- i2 <- 1
  
	for( j in 1:(l.sp1 + l.sp2) ){
	
		if(j<=l.sp1){
		
		             if(length(gam1$paraPen)!=0 && pc1==0 ){ pc1 <- pc1 + 1 
		             Ss[[j]] <- gam1$paraPen$S[[1]]
            		     rank[j] <- gam1$paraPen$rank
                              off[j] <- gam1$paraPen$off  
                                            }else{                                          
		             Ss[[j]] <- gam1$smooth[[i1]]$S[[1]]
            		     rank[j] <- gam1$smooth[[i1]]$rank
                              off[j] <- gam1$smooth[[i1]]$first.para 
                                  i1 <- i1 + 1  
                                                 }
                             }                    
                else{
		             if(length(gam2$paraPen)!=0 && pc2==0 ){ pc2 <- pc2 + 1
		             Ss[[j]] <- gam2$paraPen$S[[1]]
            		     rank[j] <- gam2$paraPen$rank
                              off[j] <- length(coef(gam1)) + gam2$paraPen$off                  
                                            }else{                  
                     	     Ss[[j]] <- gam2$smooth[[i2]]$S[[1]]  
                             off[j]  <- length(coef(gam1)) + gam2$smooth[[i2]]$first.para 
                             rank[j] <- gam2$smooth[[i2]]$rank
                                  i2 <- i2 + 1 
                                                  }
                     }                                            
       }
       
       
       
list(rank=rank,off=off,Ss=Ss)



}





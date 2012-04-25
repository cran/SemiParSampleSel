spS <- function(sp,gam1,gam2){
  di <- NA
  j  <- 1
	for(i in 1:length(gam1$sp)) di[i] <- dim(gam1$smooth[[i]]$S[[1]])[1]
	for(i in (length(gam1$sp)+1):(length(gam1$sp)+length(gam2$sp))) {di[i] <- dim(gam2$smooth[[j]]$S[[1]])[1]; j <- j + 1}
  S <- matrix(0,sum(di),sum(di))
  i <- 1
  Cum <- c(0,cumsum(di))
	for(j in 0:(length(gam1$sp)+length(gam2$sp)-1) ){
		if(j<length(gam1$sp))S[ (Cum[j+1]+1):Cum[j+2] , (Cum[j+1]+1):Cum[j+2] ] <- sp[j+1]*gam1$smooth[[j+1]]$S[[1]] 
		else{S[ (Cum[j+1]+1):Cum[j+2] , (Cum[j+1]+1):Cum[j+2] ] <- sp[j+1]*gam2$smooth[[i]]$S[[1]]
                 i <- i + 1
            } 
      }
S
}

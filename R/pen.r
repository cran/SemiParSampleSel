pen <- function(params, qu.mag, sp, VC){

    dimP1 <- dimP2 <- 0  
    
    if(VC$margins[2] %in% c("N", "NB", "PIG", "G")) dzer <- diag(0,2) 
    if(VC$margins[2] %in% c("P"))                   dzer <- diag(0,1) 
    if(VC$margins[2] %in% c("D", "S"))              dzer <- diag(0,3)  

      
    S1 <- S2 <- matrix(0,1,1)   

    S <- mapply("*", qu.mag$Ss, sp, SIMPLIFY=FALSE)
    S <- do.call(adiag, lapply(S, unlist))

    ma1 <- matrix(0,VC$gp1,VC$gp1) 
    ma2 <- matrix(0,VC$gp2,VC$gp2)

    if(length(VC$pPen1)!=0){ indP1 <- qu.mag$off[1]:(qu.mag$off[1]+qu.mag$rank[1]-1)
                          dimP1 <- length(indP1)
                          ma1[indP1,indP1] <- S[1:dimP1,1:dimP1]
                                } 
    if(length(VC$pPen2)!=0){ 
                          indP2 <- (qu.mag$off[VC$l.sp1+1]-VC$X1.d2):(-VC$X1.d2+qu.mag$off[VC$l.sp1+1]+qu.mag$rank[VC$l.sp1+1]-1)
                          dimP2 <- length(indP2)
                          ma2[indP2,indP2] <- S[(dimP1+1):(length(indP2)+dimP1),(dimP1+1):(length(indP2)+dimP1)]
                                }                                 
    
    lP1 <- length(VC$pPen1); lP2 <- length(VC$pPen2) 
    
    if((lP1!=0 && VC$l.sp1>1) || (lP1==0 && VC$l.sp1>0)) S1 <- S[(dimP1+1):(dimP1+VC$X1.d2-VC$gp1),(dimP1+1):(dimP1+VC$X1.d2-VC$gp1)]
    if((lP2!=0 && VC$l.sp2>1) || (lP2==0 && VC$l.sp2>0)){dS1 <- dim(S1)[2]; if(dS1==1) dS1 <- 0; 
                                                   S2 <- S[(dimP1+dimP2+dS1+1):dim(S)[2],(dimP1+dimP2+dS1+1):dim(S)[2]]}
    
    lS1 <- length(S1); lS2 <- length(S2) 
    
    if(lS1==1 && lS2==1) S.h <- adiag(ma1, ma2, dzer)
    if(lS1 >1 && lS2==1) S.h <- adiag(ma1, S1, ma2, dzer)
    if(lS1==1 && lS2 >1) S.h <- adiag(ma1, ma2, S2, dzer)
    if(lS1 >1 && lS2 >1) S.h <- adiag(ma1, S1, ma2, S2, dzer)
        
   
   S.h1 <- 0.5*crossprod(params,S.h)%*%params
   S.h2 <- S.h%*%params
   
   list(S.h = S.h, S.h1 = S.h1, S.h2 = S.h2)
   
   
   
         }






















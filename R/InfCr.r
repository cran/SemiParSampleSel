InfCr <- function(x,cr="AIC"){

n <- x$n; edf <- x$t.edf; l <- x$fit$l   

if(cr=="AIC")  value <- 2*l + 2*edf 
#if(cr=="AICc") value <- 2*l + 2*edf + (2*edf*(edf+1))/(n-edf-1)  
if(cr=="AICc") value <- 2*l + 2*edf*(1+(edf+1)/(n-edf-1))  
if(cr=="BIC")  value <- 2*l + log(n)*edf 

value

}










#' Function: Compute the bivariate empirical joint probability
#' 
#' @param X The vector of a monthly hydro-climatic variable of n years.
#' @param Y The vector of a monthly hydro-climatic variable of n years.  
#' 
#' @return The empirical joint probability time scale 
#' @export
#' 
#' @references Hao, Z., et al. (2017). "An integrated package for drought monitoring, prediction and analysis to aid drought modeling and assessment." Environmental Modelling & Software 91: 199-209.
#' 
#' @examples
#' X=runif(120, min = 0, max = 100)
#' Y=runif(120, min = 0, max = 100)
#' fit<-BiEmp(X,Y) 
BiEmp<-function (X,Y)
{
  X=matrix(X,ncol=1)
  Y=matrix(Y,ncol=1)
  
  n=length(X)
  
  Z=matrix(NA,nrow=length(X),ncol=1)
  
  for(k in 1:n)
  {  
    Z[k]=sum((X<=X[k])&(Y<=Y[k])) 
    
    Z[k]=(Z[k]-0.44)/(n+0.12)
    
  } 
  
  return(Z)
  
}
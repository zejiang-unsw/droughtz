#' Function: Compute the univariate empirical joint probability (EMP)
#' 
#' @param X The vector of a monthly hydro-climatic variable of n years.
#' @param dist is the funciton for the plotting position formula (Gringorten or Weibull).
#' 
#' @return The univariate EMP 
#' @export
#' 
#' @examples
#' X=runif(120, min = 0, max = 100)
#' fit<-UniEmp(X,dist = "Gringorten") 
UniEmp<-function (X,dist = "Gringorten")
{
  
  if (dist!= "Weibull" & dist != "Gringorten") 
  {
    stop("Please select either Weibull or Gringorten
         Other distribution form  will be updated soon.")
  }
  
  X=matrix(X,ncol=1)
  n=length(X)
  Z=matrix(NA,nrow=length(X),ncol=1)
  
  for  (k in 1:n)
  {  
    Z[k]=sum(X<=X[k]) 
    
    if (dist == "Weibull") 
    {
      Z[k]=Z[k]/(n+1)
    }
    else if (dist== "Gringorten") 
    {
      Z[k]=(Z[k]-0.44)/(n+0.12)
    }
  } 
  
  return(Z)
  
}
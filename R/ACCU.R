#' Function: Obtain the accumulation of monthly hydro-climatic variables
#' 
#' @param X The vector of monthly hydro-climatic variables of n years. ts is the accumulated time scale.
#' @param ts  The accumulated time scale  
#' 
#' @return The accumulated series
#' @export
#' 
#' @examples
#' X=runif(120, min = 0, max = 100)
#' Y<-ACCU(X,ts=3) # Compute the 3 month  accumulated series
ACCU<-function (X,ts=6) 
{
  
  X=matrix(X,ncol=1)
  
  X1=matrix(0,nrow=length(X)-ts+1,ncol=ts)
  
  #  define the accumulation of the time scale ts for each month
  
  X2=matrix(NA,nrow=length(X),ncol=1)
  
  for (i in 1:ts)
  {
    X1[,i]=X[i:(length(X)-ts+i)]   
  }
  
  #  Get the accumulate of the scale ts of each month
  
  X2[ts:length(X),1]=matrix(rowSums(X1),ncol=1)
  
  #  reshape the accumulation to the matrix with 12 column
  
  Y=matrix(X2,ncol=12,byrow=T)
  
  return(Y)
}

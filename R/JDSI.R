#' Function: Compute the multivariate drought index with joint distribution 
#' 
#' @param X  is The vector of a monthly hydro-climatic variable of n years. 
#' @param Y  is The vector of a monthly hydro-climatic variable of n years.
#' @param ts is the accumulated time scale. 
#' @param type is the method used to compute the JDSI (1 is Joint distribution and 2 is the Kendall function). 
#' 
#' @return The multivariate drought index of different time scales from the marginal probability (or percentile) 
#' @export
#' 
#' @references Hao, Z., et al. (2017). "An integrated package for drought monitoring, prediction and analysis to aid drought modeling and assessment." Environmental Modelling & Software 91: 199-209.
#' 
#' @examples
#' X=runif(120, min = 0, max = 100)
#' Y=runif(120, min = 0, max = 100)
#' fit<-JDSI(X,Y,ts=6)  
#' z=matrix(t(fit$JDSI),ncol=1)
#' plot(z, type="l", col=1, lwd=2, lty=1, xlim=c(0,120),xlab="Time",ylab="JDSI")
JDSI<-function (X,Y,ts=6,type=1) 
{
  X=matrix(X,ncol=1)
  Y=matrix(Y,ncol=1)
  
  if ((length(X)/12 != round(length(X)/12)))
  {
    return("ERROR: The input data should be the vector of the monthly data for a few years")  
  }
  
  # Obtain the accumulation of the variable for a specific time scale (ts) 
  AX=ACCU(X,ts)
  AY=ACCU(Y,ts)

  #  define the SI=Y0 for each month
  #  define the empirical distribution
  JDSI_A=matrix(NA,nrow=length(X)/12,ncol=12,byrow=T)
  JDSI_K=matrix(NA,nrow=length(X)/12,ncol=12,byrow=T)
  
  for (k in 1:12)
  {
    mdx=AX[,k]
    mdy=AY[,k]
    
    # if the first number is NA, exclude the first number and then compute SPI
    if (k<ts) 
    {
      xd=mdx[2:length(mdx)]
      yd=mdy[2:length(mdy)]
      
      p=BiEmp(xd,yd)
      u=copula::pobs(cbind(xd,yd))
      
      p1=p
      p2=copula::Kn(p,u)*length(xd)/(length(xd)+1)
      
      JDSI_A[2:length(mdx),k]=stats::qnorm(p1)
      JDSI_K[2:length(mdx),k]=stats::qnorm(p2)
    
    } else { 
      # if the first number is not NA, take the whole month to compute SPI
      xd=mdx
      yd=mdy    
      
      p=BiEmp(xd,yd)
      u=copula::pobs(cbind(xd,yd))
      
      p1=p
      p2=copula::Kn(p,u)*length(xd)/(length(xd)+1)

      JDSI_A[,k]=stats::qnorm(p1)
      JDSI_K[,k]=stats::qnorm(p2)
    }
  }
  
  if (type==1) 
  {
    result<-list(JDSI=JDSI_A)
  } else if (type==2) {
    result<-list(JDSI=JDSI_K)
  } else {
    stop("Only the type =1 and 2 are included. Others will be included in the future ")
  }
  
  return(result)
}

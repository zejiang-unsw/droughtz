#' Function: Compute the multivariate drought index
#' 
#' Based on the vector of a monthly hydro-climatic variable, 
#' the multivariate drought index is computed based on the marginal (or univariate) probability (or percentile) from the function SDI.
#' 
#' @param X  is The vector of a monthly hydro-climatic variable of n years. 
#' @param Y  is The vector of a monthly hydro-climatic variable of n years.
#' @param ts is the accumulated time scale. 
#' 
#' @return The multivariate drought index of different time scales from the marginal probability (or percentile) 
#' @export
#' 
#' @references Hao, Z., et al. (2017). "An integrated package for drought monitoring, prediction and analysis to aid drought modeling and assessment." Environmental Modelling & Software 91: 199-209.
#' 
#' @examples
#' 
#' X=runif(120, min = 0, max = 100)
#' Y=runif(120, min = 0, max = 100)
#' 
#' fit<-MDI(X,Y,ts=6) # Compute the 6 month drought index
#' fit$ProbEmp2 #Get the empirival drought index
MDI<-function (X,Y,ts=6) 
{
  X=matrix(X,ncol=1)
  Y=matrix(Y,ncol=1)
  
  # Obtain the accumulation of the variable for a specific time scale (ts) 
  
  XA=ACCU(X,ts)
  YA=ACCU(Y,ts)
      
      #  define the SI=Y0 for each month
      #  define the empirical distribution
      
  P_emp=matrix(NA,nrow=length(X)/12,ncol=12,byrow=T)
  P_gam=matrix(NA,nrow=length(X)/12,ncol=12,byrow=T)

  Emp2=matrix(NA,nrow=length(X)/12,ncol=12,byrow=T)
  
  
  for  (k in 1:12)
  {
  
    mdx=XA[,k]
    mdy=YA[,k]
    
    # if the first number is NA, exclude the first number and then compute SPI
    
    if (is.na(mdx[1])==TRUE) 
    {
      xd=mdx[2:length(mdx)]
      yd=mdy[2:length(mdy)]
    

      Emp2[2:length(mdx),k]=BiEmp(xd,yd)
      
    }
    
    else
      # if the first number is not NA, take the whole month to compute SPI
    {  
      
      xd=mdx
      yd=mdy    
       
      Emp2[,k]=BiEmp(xd,yd)
      
    }
    
  }
  
  #result<-list(Probemp=P_emp,Probgam=P_gam,ProbEmp2=Emp2)
   
r=stats::qnorm(Emp2)
 result<-list(ProbEmp2=r)
#result<-list(MDI=r)
  return(result)
}

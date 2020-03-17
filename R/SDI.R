#' Function: Compute the standardized drought index  
#' 
#' Based on the vector of monthly variables, the standardized drought 
#' index is computed. Note here the standardized precipitation index (SPI) is used as 
#' the example of the drought index in the univariate case. 
#' It also represents other drought indices computed in the similar way as SPI.
#' 
#' Apart from the standardized drought index, the percentile (probability) is also provided,
#' 
#' @param X The vector of a monthly hydro-climatic variable of n years. 
#' @param ts is the accumulated time scale.
#' @param dist is distribution funciton.
#' 
#' @return The univariate and multivariate drought index of different time scale from both the empirical and gamma distribution
#' @export
#' 
#' @references Hao, Z., et al. (2017). "An integrated package for drought monitoring, prediction and analysis to aid drought modeling and assessment." Environmental Modelling & Software 91: 199-209.
#' 
#' @examples
#' X=runif(120, min = 0, max = 100)
#' fit<-SDI(X,ts=3) # Compute the 3 month drought index
#' fit$SDI # Get the empirical drought index 
#' z=matrix(t(fit$SDI),ncol=1)
#' plot(z, type="l", col=1, lwd=2, lty=1, xlim=c(0,120),xlab="Time",ylab="SDI")
SDI<-function (X,ts=6,dist="EmpGrin") 
{
  
  X=matrix(X,ncol=1)
  
  if ((length(X)/12 != round(length(X)/12)))
  {
    return("ERROR: The input data should be the vector of the monthly data for a few years")  
  }
  
  
  if (dist!= "EmpGrin" & dist!= "EmpWeib" &dist != "Gamma"& dist != "Lognormal") 
  {
    stop("Only the empirical and gamma distributions are included for the current version! 
         Other distribution form  will be updated soon.
         Please select EmpGrin/EmpWeib/Gamma/Lognormal")
  }
  
  # Obtain the accumulation of the variable for a specific time scale (ts) 
  
  XA=ACCU(X,ts)
  
  #  define the SI=Y0 for each month
  #  define the empirical distribution
  
  SI_emp1=matrix(NA,nrow=length(X)/12,ncol=12,byrow=T)
  SI_emp2=matrix(NA,nrow=length(X)/12,ncol=12,byrow=T)
  SI_gam=matrix(NA,nrow=length(X)/12,ncol=12,byrow=T)
  SI_log=matrix(NA,nrow=length(X)/12,ncol=12,byrow=T)
  
  
  for  (k in 1:12)
  {
    md=XA[,k]
    
    # if the first number is NA, exclude the first number and then compute SPI
    
    if (is.na(md[1])==TRUE) 
    {
      xd=md[2:length(md)]
      
      if (dist == "EmpGrin")
      {
      emp1_cdf=UniEmp(xd,dist = "Gringorten")
      SI01 <- stats::qnorm(emp1_cdf)
      SI_emp1[2:length(md),k]=matrix(SI01,ncol=1)
          
      }
      else if(dist == "EmpWeib")
      {
      emp2_cdf=UniEmp(xd,dist = "Weibull")
      SI02 <- stats::qnorm(emp2_cdf)
      SI_emp2[2:length(md),k]=matrix(SI02,ncol=1)
      }
      else if(dist == "Gamma") 
      {
      par<-MASS::fitdistr(xd,"gamma")
      gam_cdf=stats::pgamma(xd,par$estimate[1],par$estimate[2])
      SI1 <- stats::qnorm(gam_cdf)
      SI_gam[2:length(md),k]=matrix(SI1,ncol=1)
      }
      else  (dist == "Lognormal")  
      {
      par<-MASS::fitdistr(xd,"lognormal")
      log_cdf=stats::plnorm(xd,par$estimate[1],par$estimate[2])
      SI2 <- stats::qnorm(log_cdf)
      SI_log[2:length(md),k]=matrix(SI2,ncol=1)
      }
      
      

    }
    
    else
      # if the first number is not NA, take the whole month to compute SPI
    {  xd=md
       
    if (dist == "EmpGrin")
    {
      emp1_cdf=UniEmp(xd,dist = "Gringorten")
      SI01 <- stats::qnorm(emp1_cdf)
      SI_emp1[,k]=matrix(SI01,ncol=1)
      
    }
    else if(dist == "EmpWeib")
    {
      emp2_cdf=UniEmp(xd,dist = "Weibull")
      SI02 <- stats::qnorm(emp2_cdf)
      SI_emp2[,k]=matrix(SI02,ncol=1)
    }
    else if(dist == "Gamma") 
    {
      par<-MASS::fitdistr(xd,"gamma")
      gam_cdf=stats::pgamma(xd,par$estimate[1],par$estimate[2])
      SI1 <- stats::qnorm(gam_cdf)
      SI_gam[,k]=matrix(SI1,ncol=1)
    }
    else if (dist == "Lognormal")  
    {
    par<-MASS::fitdistr(xd, "lognormal")
    log_cdf=stats::plnorm(xd,par$estimate[1],par$estimate[2])
    SI2 <- stats::qnorm(log_cdf)
    SI_log[,k]=matrix(SI2,ncol=1)
    }
    
    }
    
  }
  
  # result<-list(Probemp=P_emp,SDIemp=SI_emp,Probgam=P_gam,SDIgam=SI_gam)
  
  if (dist=="EmpGrin") 
    
  {result<-list(SDI=SI_emp1) }
  
  else if (dist=="EmpWeib")
  {
 {result<-list(SDI=SI_emp2) }
  }
  else if (dist =="Gamma")
    
  {result<-list(SDI=SI_gam) }
  
  else if (dist =="Lognormal")
    
  {result<-list(SDI=SI_log) }
  
  else
  {
    stop("Only the empirical/gamma/weibull distributions are included! Please select either emp or gam ")
  }
  
  return(result)
}

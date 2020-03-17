#' Function: Univariate and multivariate return period (Gumbel copula)
#' 
#' @param X   the drought properties or indices
#' @param Y   the drought properties or indices
#' @param EL  the average reocurrence time, default EL=1 (annual maxima)
#' 
#' @return The univariate and multivariate return period
#' @export
#' 
#' @references Hao, Z., et al. (2017). "An integrated package for drought monitoring, prediction and analysis to aid drought modeling and assessment." Environmental Modelling & Software 91: 199-209.
#' @references Shiau, J.-T. and R. Modarres (2009). "Copula‐based drought severity‐duration‐frequency analysis in Iran." Meteorological Applications: A journal of forecasting, practical applications, training techniques and modelling 16(4): 481-489.
#' 
#' @examples
#' 
UMFreq<-function (X,Y,EL=1) 
{
  n=length(X)
  if(sum(X<=0)|sum(Y<=0)) stop("duration or serverity smaller than 0!")
  
  # Compute the univariate return period
  UT1<-matrix(NA, nrow=n,ncol=1)
  UT2<-matrix(NA, nrow=n,ncol=1)
  
  pa<-MASS::fitdistr(X,"exponential")
  P1=stats::pexp(X,pa$estimate)
  UT1=1/(1-P1)*EL               # Exceedance return period
  
  pg<-MASS::fitdistr(Y,"gamma")
  P2=stats::pgamma(Y,shape=pg$estimate[1],rate=pg$estimate[2])
  UT2=1/(1-P2)*EL               # Exceedance return period
  
  # Compute the joint return period using "copula" package
  d=2                           # the dimension of the copula
  u=copula::pobs(cbind(X,Y))
  
  ## maximum pseudo-likelihood, mpl
  fit <-copula::fitCopula(copula::gumbelCopula(), u, method="mpl", start=1)
  theta<-stats::coef(fit)
 
  cop <- copula::gumbelCopula(theta, dim=d)
  
  copk <- copula::onacopulaL("Gumbel", list(theta, 1:2))
  
  MTA<-matrix(NA, nrow=n,ncol=1) #joint return period for AND case
  MTO<-matrix(NA, nrow=n,ncol=1) #joint return period for OR case
  MTK<-matrix(NA, nrow=n,ncol=1) #joint return period based on Kendall distribution
  
  for (i in 1:n){
    P=copula::pCopula(c(P1[i],P2[i]),cop) 
      
    MTA[i]= 1/(1-P1[i]-P2[i]+P)*EL  
    MTO[i]= 1/(1-P)*EL  
    
    Kg <- copula::pK(P, cop=copk@copula, 2)
    MTK[i]= 1/(1-Kg)*EL
  }
  
  result<-list(pa=pa,pg=pg,cop=cop,
               P1=P1, P2=P2,
               UT1=UT1,UT2=UT2,
               MTA=c(MTA),MTO=c(MTO),MTK=c(MTK))
  return(result)
  
}

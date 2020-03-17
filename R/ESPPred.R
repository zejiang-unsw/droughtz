#' Function: Drought prediction with ESP method
#' 
#' @param X is the monthly variables.
#' @param Y is the monthly variables.
#' @param L is the lead time. 
#' @param m is the monthly variables. 
#' @param ts is the monthly variables.
#' 
#' @return The prediction of univariate and multivariate drought index
#' @export
#' 
#' @examples
#' X=runif(120, min = 0, max = 100)
#' Y=runif(120, min = 0, max = 100)
#' ESPPred(X,Y,L=1,m=7,ts=6)
ESPPred<-function (X,Y,L=1,m=7,ts=6) 
{
  
  # use 6-month time scale as an example
  # Note L+m<=12 (i.e., do not predict the next year)
  
  # Load the monthly data
  xd=as.vector(X)
  yd=as.vector(Y)
  
  ny=length(xd)/12 
  
  mdx=matrix(xd,nrow=ny,ncol=12,byrow=T)
  mdy=matrix(yd,nrow=ny,ncol=12,byrow=T)
  
  #Accumulation
  m_ax=ACCU(xd,ts)
  m_ay=ACCU(yd,ts)
  
  
  #Resampling
  rdx=matrix(NA,nrow=ny-1,ncol=L)
  rdy=matrix(NA,nrow=ny-1,ncol=L)
  
  if (L==1)
  {
    mk=m+1

    rdx[,1]=mdx[1:ny-1,mk]
    rdy[,1]=mdy[1:ny-1,mk]
  }

  if (L==2)
  {
    mk=matrix(c(m+1,m+2),nrow=1,ncol=2)
    
    rdx[,1:2]=mdx[1:ny-1,mk[,1]:mk[,2]]
    rdy[,1:2]=mdy[1:ny-1,mk[,1]:mk[,2]]
  }
  
  #Obs.
  n0=length(xd)-12+m
  
  odx=sum(xd[(n0-(ts-L-1)):n0])
  ody=sum(yd[(n0-(ts-L-1)):n0])
  
  #Blend pred.
  adx=matrix(NA,nrow=ny-1,ncol=1)
  ady=matrix(NA,nrow=ny-1,ncol=1)
  
  for  (i in 1:ny-1)
  {
    s1=rowSums(rdx)
    s2=rowSums(rdy)
    
    adx[i,1]=odx+s1[i];
    ady[i,1]=ody+s2[i];
  }
  
  #Historical
  mo=m+L
  
  if (m>=ts){
    ohx=matrix(NA,nrow=ny-1,ncol=1)
    ohy=matrix(NA,nrow=ny-1,ncol=1)
    
    ohx[,1]=m_ax[1:ny-1,mo]
    ohy[,1]=m_ay[1:ny-1,mo]
  } else {
    
    ohx=matrix(NA,nrow=ny-2,ncol=1)
    ohy=matrix(NA,nrow=ny-2,ncol=1)
    
    ohx[,1]=m_ax[2:ny-1,mo]
    ohy[,1]=m_ay[2:ny-1,mo]
  }

  #Pred. of univariate index
  SDIx=matrix(NA,nrow=ny-1,ncol=1)
  SDIy=matrix(NA,nrow=ny-1,ncol=1)
  DIxy=matrix(NA,nrow=ny-1,ncol=1)
  
  for(i in 1:length(ohx))
  {
    cdx=as.vector(c(ohx,adx[i]))
    cdy=as.vector(c(ohy,ady[i]))
    
    #Pred. of univariate index SDIx
    empcdf <- stats::ecdf(cdx) 
    cdf0 <- sapply(cdx, empcdf)
    px<-cdf0*length(cdx)/(length(cdx)+1)
    DIx<-stats::qnorm(px)
    
    SDIx[i,1]<-DIx[length(DIx)]
    
    #Pred. of univariate index SDIy
    empcdf <- stats::ecdf(cdy) 
    cdf0 <- sapply(cdy, empcdf)
    py<-cdf0*length(cdy)/(length(cdy)+1)
    DIy<-stats::qnorm(py)
    
    SDIy[i,1]<-DIy[length(DIy)]
    
    #Pred. of multivarite index
    p2=BiEmp(cdx,cdy)
    DI2<-stats::qnorm(p2)
    
    DIxy[i,1]<-DI2[length(DI2)]
    
  }
  
  xp=stats::median(SDIx)
  yp=stats::median(SDIy)
  xyp=stats::median(DIxy)
  
  
  result<-list(SDIx=SDIx,SDIy=SDIy,DIxy=DIxy,xp=xp,yp=yp,xyp=xyp)
  return(result)

}

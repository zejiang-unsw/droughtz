#' Function: Compute fixed drought threshold
#'
#' @param data        A vector of time series data
#' @param percentile  The percentile of FDC (exceedance), range from (0,1)
#'
#' @return
#' @export
#' 
#' @references Van Loon, A. F. (2015). "Hydrological drought explained." Wiley Interdisciplinary Reviews: Water 2(4): 359-392.
#'
#' @examples
#' 
FixThres <- function(data, percentile){
  perc = 1-percentile
  fthreshold = array()
  fthreshold = quantile(data, probs=(perc), na.rm=T)
  
  return(list(FixThres=fthreshold))
}

#' Function: Compute varying drought threshold
#'
#' @param data        A vector of time series data
#' @param percentile  The percentile of FDC (exceedance), range from (0,1)
#' @param dates       A vector of date variable
#' @param by          Time step of the input data set
#'
#' @return
#' @export
#' 
#' @references Van Loon, A. F. (2015). "Hydrological drought explained." Wiley Interdisciplinary Reviews: Water 2(4): 359-392.
#'
#' @examples
#' 
VarThres <- function(data, percentile, dates, by="month"){
  perc = 1-percentile
  vthreshold = array()
  thres = array()
  
  if(by=="day"){
    # Julian days
    Jdays = strptime(dates,format="%Y-%m-%d")$yday+1 
  
    # data filtering by 30 days moving window
    datafilter = filter(data, rep(1/30,30),sides=2)
    for(d in 1:366){
      sel <- which(Jdays == d)
      vthreshold[d] = quantile(datafilter[sel], probs=(perc), na.rm=T)
    }
    
    for(i in 1:length(Jdays)){
      thres[i] = vthreshold[Jdays[i]]
    }
    
  } else if(by=="month"){
    # the number of month
    month = strptime(dates,format="%Y-%m-%d")$mon+1
    
    for(d in 1:12){
      sel <- which(month == d)
      vthreshold[d] = quantile(data[sel], probs=(perc), na.rm=T)
    }
    
    for(i in 1:length(month)){
      thres[i] = vthreshold[month[i]]
    }
  }
  
  return(list(VarThres=thres))
}

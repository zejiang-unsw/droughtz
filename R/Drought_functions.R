#' Function: Compute drought characteristics
#'
#' @param X       A vector of time series data
#' @param start   Start of time series 
#' @param end     End of time series
#' @param theta   Threshold of a drought event
#'
#' @return
#' @export
#' @import zoo
#' 
#' @references Van Loon, A. F. (2015). "Hydrological drought explained." Wiley Interdisciplinary Reviews: Water 2(4): 359-392.
#'
#' @examples
#' 
drought <- function(X, start, end, theta){

  data = window(X, start=start, end = end)
  dates = as.Date(as.yearmon(time(data)))
  if(length(theta)==1) theta = rep(theta,length(data))
    
  dur = array() #duration
  def = array() #deficit volume, one measure of severity
  sev = array() #severity
  st  = array() #start of the event
  end = array() #end of the event
  p = 0         #number of events
  
  sel = which(data < theta)
  int = data - theta            #intensity: lowest SPI value of the drought event
  last = -999

  for(i in sel){
    if((i-1) != last){ #if still the same drought event
      p = p + 1
      st[p] = as.character(dates[i])
      dur[p] = 0
      def[p] = 0
      sev[p] = 0
      end[p] = 0
    }
    last = i
    dur[p] = dur[p] + 1
    def[p] = def[p] + data[i] - theta[i]  #severity by the accumulated value of the event
    sev[p] = min(int[i],sev[p])           #severity by the min value of the event
    end[p] = as.character(dates[i])
  }
  
  return(list(start=st, end=end, 
              dur=dur, 
              def=signif(-def, digits=4), 
              sev=signif(-sev, digits=4)))
    
}
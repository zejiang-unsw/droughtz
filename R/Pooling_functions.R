#' Function: Pooling droughts
#'
#' @param drought   drought characteristics data variable
#' @param IT        interevent time period 
#' @param by        time step
#'
#' @return Pooled droughts by combining drought events within interevent time period
#' @export
#'
#' @examples
#' 
Pooling <- function(drought, IT=1, by="month"){
  
  st = drought$start
  end = drought$end
  dur = drought$dur
  def = drought$def
  maxint = drought$sev
  
  if(!is.na(st)){
  
  start.date = as.Date(st,format="%Y-%m-%d")
  end.date = as.Date(end,format="%Y-%m-%d")

  elapsed_months <- function(end_date, start_date) {
    ed <- as.POSIXlt(end_date)
    sd <- as.POSIXlt(start_date)
    12 * (ed$year - sd$year) + (ed$mon - sd$mon)
  }

  intereventtime = array()
  for(i in 1:length(dur)){
    if(by=="day"){
      intereventtime[i] = (start.date[i+1])-(end.date[i]+1)
    } else if(by=="month"){
      intereventtime[i] = elapsed_months(start.date[i+1],end.date[i])-1
    }
  }
  intereventtime[length(dur)] = 999 
  dur.pooled = array()
  def.pooled = array()
  maxint.pooled = array()
  st.pooled = array()
  end.pooled = array()
  duration = dur[1]
  deficit = def[1]
  maxintensity = maxint[1]
  begin = start.date[1]
  j=1
  for(i in 1:length(dur)){
    if (intereventtime[i] < IT){
      duration = duration + intereventtime[j] + dur[j+1]
      deficit = deficit + def[j+1]
      maxintensity = max(maxintensity, maxint[j+1])
      if (j == 1 || intereventtime[j-1] > IT) {
        begin = start.date[j]
      }
      j = j+1
    } else { 
      dur.pooled[i] = duration
      def.pooled[i] = deficit
      maxint.pooled[i] = maxintensity
      st.pooled[i] = begin
      j=i+1
      duration = dur[i+1]
      deficit = def[i+1]
      maxintensity = maxint[i+1]
      begin = start.date[i+1]
    }
  }
  
  st.pooled = as.Date(st.pooled,format="%Y-%m-%d",origin="1970-01-01")
  
  end.pooled = sapply(which(!is.na(st.pooled)), function(i) seq(st.pooled[i], by = by, length = dur.pooled[i])[dur.pooled[i]])
  end.pooled = as.Date(end.pooled,format="%Y-%m-%d",origin="1970-01-01")
  
  # droughts.pooled = na.omit(data.frame(st.pooled,end.pooled,dur.pooled,def.pooled,maxint.pooled))
  # row.names(droughts.pooled) <- 1:nrow(droughts.pooled)
  # return(as.list(droughts.pooled))

  return(list(start=na.omit(st.pooled), end=end.pooled, 
               dur=dur.pooled[which(!is.na(st.pooled))], 
               def=def.pooled[which(!is.na(st.pooled))], 
               sev=maxint.pooled[which(!is.na(st.pooled))]))
 } else {
  return(drought)
 }
}

#' Function: Remove minor droughts 
#'
#' @param file drought characteristics data variable
#' @param dur  drought duration variable
#' @param MIN  minimum drought duration
#'
#' @return A subset of droughts with duration larger than x days or months
#' @export
#'
#' @examples
#' 
Minor <- function(file, dur, MIN){
  
  #dur=file$dur
  droughts.sub <- lapply(file, function(i) subset(i, dur > MIN)) 

  return(droughts.sub)
}


#' Function: Get linear regression equations
#'
#' @param df drought characteristics data variable
#'
#' @return A string of equation, including intercept, coefficients and statistics.
#' @export
#' 
#' @examples
#' 
lm_eqn <- function(df){
  y <- df$Sev; x <- df$Dur
  
  if(!is.na(y[1])){
    m <- lm(y ~ x, df)
    # eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
    #      list(a = format(unname(coef(m)[1]), digits = 2),
    #           b = format(unname(coef(m)[2]), digits = 2),
    #          r2 = format(summary(m)$r.squared, digits = 3)))
    
    eq <- substitute(~~R^2~"="~r2, 
                     list(r2 = format(round(summary(m)$r.squared,3), digits = 3, nsmall = 3, justify="left")))
    return(as.character(as.expression(eq)))
  } else {
    return(NA)
  }
}
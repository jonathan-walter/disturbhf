#' Calculate the decimal day of year
#'
#' \code{decimal_doy} computes the decimal day of year
#'
#' @param datetime a date-time object in POSIX format
#'
#' @return \code{decimal_doy} returns a numeric vector with the decimal day of year
#'
#' @description The decimal day of year uses hours, minutes, and seconds to compute the fractional
#' day, expressed as a decimal.
#'
#' @author Jonathan Walter, \email{jaw3es@@virginia.edu}
#'
#' @export

decimal_doy<-function(datetime){
  dd<-as.numeric(strftime(datetime, format="%j"))
  hh<-as.numeric(strftime(datetime, format="%H"))
  mm<-as.numeric(strftime(datetime, format="%M"))
  ss<-as.numeric(strftime(datetime, format="%S"))

  totsec<-24*60*60
  dec<-(hh*60*60 + mm*60 + ss)/totsec

  return(dd + dec)

}

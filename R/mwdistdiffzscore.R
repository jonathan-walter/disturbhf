
#' Calculate differences in continuous distribution functions between reference and moving window observations
#'
#' \code{mwdistdiffz} computes differences between the continuous distribution functions (cdf) for observations within
#' a moving window and a reference distribution, and the z-scores of differences relative to samples of the reference distribution.
#' It is used for identifying recovery times from disturbance in time series data.
#'
#' @param testy a data frame representing the "test" period, i.e., the period over which to search for disturbance and recovery.
#' It must contain the columns \code{tt} and \code{yy}.
#' \code{tt} contains numeric times assumed to be in units of days, or POSIXct (date) format corresponding to observations of \code{yy}.
#' If \code{is.numeric(tt)}, tt is taken to be in units of day-of-year, with January 1 = doy 1.
#' \code{yy} is a numeric time series.
#' @param refy a data frame representiong the refrence period that \code{testy} is compared to.
#' It must contain the columns \code{tt} and \code{yy}, as in \code{testy}.
#' @param wwidth the moving window width, in number of time steps
#' @param refwidth the width of the rolling reference window, in number of time steps; if \code{NULL} do not use rolling reference window
#' @param dx increment between values at which to evaluate differences between the cdf for \code{yy} and \code{refy}
#' @param stride number of time steps by which the moving window advances
#' @param dmin the fraction of data that must be present (i.e., non-NA) in test (and, if applicable, adaptive reference) moving windows to procede with computations.
#' @param ddiff_method statistic quantifying difference between test and reference distributions; default is "dist" for maximum absolute difference between the ecdf's (bounded [0,1]), "integral" computes the integral of the absolute value of the differences between the ecdf's over the range of observed value (bounded by the distance between the distribution if they are completely non-overlapping)
#'
#' @return \code{mwdistdiff} returns a data frame containing the columns:
#' \item{tleft}{the time corresponding to the beginning (left edge) of the moving window}
#' \item{tright}{the time corresponding to the end (right edge) of the moving window}
#' \item{ddiff}{differences from the reference distritbution}
#' \item{zz}{Z-scores representing the strength of excursion from reference}
#'
#' @details The value suppplied to \code{refwidth} determines whether a rolling reference window will be used,
#' or if all of refy will be used as the reference period.
#' \code{testy} and \code{refy} must be regularly spaced, ordered time series. If there are gaps in the time series, fill
#' them with NAs prior to applying this function. If there are insufficient data then output will contain \code{NA}s.
#' \code{dmin} takes into account NAs in time series as well as reference windows overlapping the edges or gaps in the time series.
#' Currently, the year is assumed to be 365 days long (i.e., there is no internal handling of leap years).
#'
#' @author Jonathan Walter, \email{jaw3es@@virginia.edu}
#'
#' @examples
#' #need to add some
#'
#' @export

## create DOY column to index off of, accept tt as numeric (1,Inf), or date
mwdistdiffz<-function(testy, refy, wwidth, refwidth=NULL, dx=0.01, stride=1, dmin=0.5, ddiff_method = "dist"){

  #a little error handling
  if(!is.data.frame(testy) | !"tt" %in% colnames(testy) | !"yy" %in% colnames(testy)){
    stop("testy must be a data frame containing the columns 'tt' and 'yy'")
  }
  if(!is.data.frame(refy) | !"tt" %in% colnames(refy) | !"yy" %in% colnames(refy)){
    stop("refy must be a data frame containing the columns 'tt' and 'yy'")
  }
  if(any(class(testy$tt) != class(refy$tt))){
    stop("testy$tt and refy$tt must have same format")
  }

  #TODO: Check that testy is regularly spaced
  #warnings
  if(nrow(testy)<50){
    warning("the test period is very short, results may not be robust")
  }
  if(nrow(refy)<50){
    warning("the referece period is very short, results may not be robust")
  }

  #Objects that stay the same regardless of inputs
  xx<-seq(min(c(refy$yy,testy$yy),na.rm=T)-dx, max(c(refy$yy,testy$yy),na.rm=T)+dx, dx)

  #------------------------------------------------------------------------------------------------
  # Using ALL of refy as the reference period (no seasonal adaptive referecne)

  if(is.null(refwidth)){
    if(is.numeric(testy$tt)){
      #Get mu and sd for excursions during the reference period
      wind0<-seq(from=ceiling(wwidth/2)+1, to=length(refy$tt)-ceiling(wwidth/2), by=stride)
      ddiff0<-rep(NA, length(wind0))
      dt<-diff(testy$tt)[1]
      for(ww in 1:length(wind0)){
        pd<-(wind0[ww]-(wwidth/2)):(wind0[ww]+wwidth/2-1) #get period for window
        if(mean(!is.na(refy$yy[pd])) < dmin){
          ddiff0[ww]<-NA
        }
        # else if(any(diff(diff(refy$tt[pd]))>1e-6)){
        #   #returns NA if any times are skipped
        #   ddiff0[ww]<-NA
        # }
        else{
          tmprefdist<-ecdf(refy$yy[-pd])
          tmpwdist<-ecdf(refy$yy[pd])
          if(ddiff_method == "dist"){
            ddiff0[ww]<-max(abs(tmprefdist(xx)-tmpwdist(xx)))
          }else if(ddiff_method == "integral"){
            ddiff0[ww]<-dx*sum(abs(tmprefdist(xx)-tmpwdist(xx)))
          }

        }
      }
      #rm(tmprefdist,tmpwdist)
      mu.ref<-mean(ddiff0, na.rm=T)
      sd.ref<-sd(ddiff0, na.rm=T)

      #Compute Z-scores for excursions during the test period vs. reference
      refdist<-ecdf(refy$yy) #make ecdf for reference period

      wind<-seq(from=ceiling(wwidth/2)+1, to=length(testy$tt)-ceiling(wwidth/2), by=stride)
      ddiff<-rep(NA, length(wind))
      zz<-rep(NA, length(wind))

      for(ww in 1:length(wind)){
        pd<-(wind[ww]-(wwidth/2)):(wind[ww]+wwidth/2-1) #get period for window
        if(mean(!is.na(testy$yy[pd])) < dmin){
          ddiff[ww]<-NA
        }
        else{
          wdist<-ecdf(testy$yy[pd])
          if(ddiff_method == "dist"){
            ddiff[ww]<-max(abs(refdist(xx)-wdist(xx)))
          }else if(ddiff_method == "integral"){
            ddiff[ww]<-dx*sum(abs(refdist(xx)-wdist(xx)))
          }
          zz[ww]<-(ddiff[ww]-mu.ref)/sd.ref
        }
      }
      wleft<-testy$tt[wind]-wwidth*dt/2
      wright<-testy$tt[wind]+wwidth*dt/2
    }

    if(any(grepl("POSIX",class(refy$tt)))){

      dt<-diff(testy$tt)[1]

      #Get mu and sd for excursions during the reference period
      wind0<-seq(from=ceiling(wwidth/2)+1, to=length(refy$tt)-ceiling(wwidth/2), by=stride)
      ddiff0<-rep(NA, length(wind0))

      for(ww in 1:length(wind0)){
        pd<-(wind0[ww]-(wwidth/2)):(wind0[ww]+wwidth/2-1) #get period for window
        if(mean(!is.na(refy$yy[pd])) < dmin){
          ddiff0[ww]<-NA
        }
        # else if(any(diff(diff(refy$tt[pd]))>1e-6)){
        #   #returns NA if any times are skipped
        #   ddiff0[ww]<-NA
        # }
        else{
          tmprefdist<-ecdf(refy$yy[-pd])
          tmpwdist<-ecdf(refy$yy[pd])
          if(ddiff_method == "dist"){
            ddiff0[ww]<-max(abs(tmprefdist(xx)-tmpwdist(xx)))
          }else if(ddiff_method == "integral"){
            ddiff0[ww]<-dx*sum(abs(tmprefdist(xx)-tmpwdist(xx)))
          }

        }
      }
      #rm(tmprefdist,tmpwdist)
      mu.ref<-mean(ddiff0, na.rm=T)
      sd.ref<-sd(ddiff0, na.rm=T)

      #Compute Z-scores for excursions during the test period vs. reference
      refdist<-ecdf(refy$yy) #make ecdf for reference period

      wind<-seq(from=ceiling(wwidth/2)+1, to=length(testy$tt)-ceiling(wwidth/2), by=stride)
      # tleft<-refy$tt[wind-wwidth/2]
      # tright<-refy$tt[wind+wwidth/2]
      ddiff<-rep(NA, length(wind))
      zz<-rep(NA, length(wind))

      for(ww in 1:length(wind)){
        pd<-(wind[ww]-(wwidth/2)):(wind[ww]+wwidth/2-1) #get period for window
        if(mean(!is.na(testy$yy[pd])) < dmin){
          ddiff[ww]<-NA
        }
        else{
          wdist<-ecdf(testy$yy[pd])
          if(ddiff_method == "dist"){
            ddiff[ww]<-max(abs(refdist(xx)-wdist(xx)))
          }else if(ddiff_method == "integral"){
            ddiff[ww]<-dx*sum(abs(refdist(xx)-wdist(xx)))
          }
          zz[ww]<-(ddiff[ww]-mu.ref)/sd.ref
        }
      }
      wleft<-testy$tt[wind]-wwidth*dt/2#*24*60*60
      wright<-testy$tt[wind]+wwidth*dt/2#*24*60*60
    }

  }

  # -----------------------------------------------------------------------------------------------
  # This is with the seasonal adaptive reference window

  if(!is.null(refwidth)){

    #for numeric tt ...
    if(is.numeric(testy$tt)){

      #Compute day of year corresponding to tt, Jan 1=doy 1
      testy$doy <- testy$tt %% 365
      testy$doy[testy$doy==0]<-365
      refy$doy <- refy$tt %% 365
      refy$doy[refy$doy==0]<-365
      dt<-diff(testy$doy)[1]
      tmin<-min(refy$tt)
      tmax<-max(refy$tt)



      #Compute excursions in ref period and get mean and sd
      wind<-seq(from=tmin, to=tmax, by=stride*dt)
      ddiff<-rep(NA, length(wind))

      for(ww in 1:length(wind)){

        wind.ww<-wind[ww] %% 365
        if(wind.ww==0){wind.ww<-365}

        if(wind.ww < min(testy$doy) | wind.ww > max(testy$doy)){next}

        ltest<-wind[ww]-wwidth*dt/2 #left side of "test" window
        if(ltest < tmin){next} #skip indices where window overhangs beginning of time series
        rtest<-wind[ww]+wwidth*dt/2 #right side of "test" window
        if(rtest > tmax){break} #stop computation when window overhangs end of time series
        if(rtest > ltest){
          tpd<-refy$tt > ltest & refy$tt <= rtest
        }
        if(rtest < ltest){
          tpd<-refy$tt > ltest | refy$tt <= rtest
        }

        lref<-min(refy$doy[abs(refy$tt-wind[ww]) < dt/10]-refwidth*dt/2) %% 365 #left side of reference window
        if(lref==0){lref==365}
        rref<-max(refy$doy[abs(refy$tt-wind[ww]) < dt/10]+refwidth*dt/2) %% 365 #right side of reference window
        if(rref==0){rref==365}
        if(rref > lref){
          rpd<-refy$doy > lref & refy$doy <= rref
        }
        if(rref < lref){
          rpd<-refy$doy > lref | refy$doy <= rref
        }

        #check if sufficient non-missing values in reference and test periods
        if(sum(!is.na(refy$yy[tpd]))/wwidth < dmin | sum(!is.na(refy$yy[rpd]))/refwidth < dmin){
          next
        }
        refdist<-ecdf(refy$yy[rpd])
        wdist<-ecdf(refy$yy[tpd])
        if(ddiff_method == "dist"){
          ddiff[ww]<-max(abs(refdist(xx)-wdist(xx)))
        }else if(ddiff_method == "integral"){
          ddiff[ww]<-dx*sum(abs(refdist(xx)-wdist(xx)))
        }
      }
      mu.ref<-mean(ddiff, na.rm=T)
      sd.ref<-sd(ddiff, na.rm=T)

      #Compute excursions in test period and get z-score
      tmin<-min(testy$tt)
      tmax<-max(testy$tt)
      wind<-seq(from=tmin, to=tmax, by=stride*dt)
      ddiff<-rep(NA, length(wind))
      zz<-rep(NA, length(wind))

      for(ww in 1:length(wind)){

        ltest<-wind[ww]-wwidth*dt/2 #left side of test window
        if(length(lref)==0){next}#enables skipping indices when time series are gappy
        if(ltest < tmin){next} #skip indices where window overhangs beginning of time series
        rtest<-wind[ww]+wwidth*dt/2 #right side of test window
        if(rtest > tmax){break} #stop computation when window overhands end of time series
        if(rtest > ltest){
          tpd<-testy$tt > ltest & testy$tt <= rtest
        }
        if(rtest < ltest){
          tpd<-testy$tt > ltest | testy$tt <= rtest
        }

        lref<-min(testy$doy[abs(testy$tt-wind[ww]) < dt/10]-refwidth*dt/2) %% 365 #left side of reference window
        if(length(lref)==0){next}#enables skipping indices when time series are gappy
        if(lref==0){lref==365}
        rref<-max(testy$doy[abs(testy$tt-wind[ww]) < dt/10]+refwidth*dt/2) %% 365 #right side of reference window
        if(rref==0){rref==365}
        if(rref > lref){
          rpd<-refy$doy > lref & refy$doy <= rref
        }
        if(rref < lref){
          rpd<-refy$doy > lref | refy$doy <= rref
        }

        #check if sufficient non-missing values in reference and test periods
        if(sum(!is.na(testy$yy[tpd]))/wwidth/dt < dmin | sum(!is.na(refy$yy[rpd]))/refwidth/dt < dmin){
          next
        }

        refdist<-ecdf(refy$yy[rpd])
        wdist<-ecdf(testy$yy[tpd])
        if(ddiff_method == "dist"){
          ddiff[ww]<-max(abs(refdist(xx)-wdist(xx)))
        }else if(ddiff_method == "integral"){
          ddiff[ww]<-dx*sum(abs(refdist(xx)-wdist(xx)))
        }
        zz[ww]<-(ddiff[ww]-mu.ref)/sd.ref

      }
      wleft<-wind-wwidth*dt/2
      wright<-wind+wwidth*dt/2
    }

    #for tt as a date
    if(any(grepl("POSIX", class(testy$tt)))){

      testy$doy <- decimal_doy(testy$tt) #compute decimal day of year corresponding to tt
      refy$doy <- decimal_doy(refy$tt)
      dt<-diff(testy$doy)[1] #time difference in decimal days
      dtt<-diff(testy$tt)[1] #time difference in time units
      tmin<-min(refy$tt) # minimum time in test period
      tmax<-max(refy$tt) # maximum time in test period

      #Compute excursions in ref period and get mean and sd
      wind<-seq(from=tmin, to=tmax, by=stride*dtt)
      ddiff<-rep(NA, length(wind))

      for(ww in 1:length(wind)){

        if(floor(decimal_doy(wind[ww])) < min(testy$doy) | floor(decimal_doy(wind[ww])) > max(testy$doy)){next} #skip if doy of window is not part of test period

        ltest<-wind[ww]-wwidth*dtt/2 #left side of "test" window
        if(ltest < tmin){next} #skip indices where window overhangs beginning of time series
        rtest<-wind[ww]+wwidth*dtt/2 #right side of "test" window
        if(rtest > tmax){break} #stop computation when window overhangs end of time series
        if(rtest > ltest){
          tpd<-refy$tt > ltest & refy$tt <= rtest
        }
        if(rtest < ltest){
          tpd<-refy$tt > ltest | refy$tt <= rtest
        }

        lref<-min(refy$doy[abs(refy$tt-wind[ww]) < dtt/10]-refwidth*dt/2) %% 365 #left side of reference window
        if(length(lref)==0){next} #enables skipping indices when time series are gappy
        if(lref==0){lref==365}
        rref<-max(refy$doy[abs(refy$tt-wind[ww]) < dtt/10]+refwidth*dt/2) %% 365 #right side of reference window
        if(rref==0){rref==365}
        if(rref > lref){
          rpd<-refy$doy > lref & refy$doy <= rref
        }
        if(rref < lref){
          rpd<-refy$doy > lref | refy$doy <= rref
        }

        #check if sufficient non-missing values in reference and test periods
        if(sum(!is.na(refy$yy[tpd]))/wwidth < dmin | sum(!is.na(refy$yy[rpd]))/refwidth < dmin){
          next
        }
        refdist<-ecdf(refy$yy[rpd])
        wdist<-ecdf(refy$yy[tpd])
        if(ddiff_method == "dist"){
          ddiff[ww]<-max(abs(refdist(xx)-wdist(xx)))
        }else if(ddiff_method == "integral"){
          ddiff[ww]<-dx*sum(abs(refdist(xx)-wdist(xx)))
        }

      }
      mu.ref<-mean(ddiff, na.rm=T)
      sd.ref<-sd(ddiff, na.rm=T)

      #Compute excursions in test period and get z-score
      tmin<-min(testy$tt) # minimum time in test period
      tmax<-max(testy$tt) # maximum time in test period
      wind<-seq(from=tmin, to=tmax, by=stride*dtt)
      ddiff<-rep(NA, length(wind))
      zz<-rep(NA, length(wind))


      for(ww in 1:length(wind)){

        ltest<-wind[ww]-wwidth*dtt/2 #left side of test window
        if(ltest < tmin){next} #skip indices where window overhangs beginning of time series
        rtest<-wind[ww]+wwidth*dtt/2 #right side of test window
        if(rtest > tmax){break} #stop computation when window overhands end of time series
        if(rtest > ltest){
          tpd<-testy$tt > ltest & testy$tt <= rtest
        }
        if(rtest < ltest){
          tpd<-testy$tt > ltest | testy$tt <= rtest
        }

        lref<-min(testy$doy[abs(testy$tt-wind[ww]) < dtt/10]-refwidth*dt/2) %% 365 #left side of reference window
        if(length(lref)==0){next}#enables skipping indices when time series are gappy
        if(lref==0){lref==365}
        rref<-max(testy$doy[abs(testy$tt-wind[ww]) < dtt/10]+refwidth*dt/2) %% 365 #right side of reference window
        if(rref==0){rref==365}
        if(rref > lref){
          rpd<-refy$doy > lref & refy$doy <= rref
        }
        if(rref < lref){
          rpd<-refy$doy > lref | refy$doy <= rref
        }

        #check if sufficient non-missing values in reference and test periods
        if(sum(!is.na(testy$yy[tpd]))/wwidth < dmin | sum(!is.na(refy$yy[rpd]))/refwidth < dmin){
          next
        }

        refdist<-ecdf(refy$yy[rpd])
        wdist<-ecdf(testy$yy[tpd])
        if(ddiff_method == "dist"){
          ddiff[ww]<-max(abs(refdist(xx)-wdist(xx)))
        }else if(ddiff_method == "integral"){
          ddiff[ww]<-dx*sum(abs(refdist(xx)-wdist(xx)))
        }
        zz[ww]<-(ddiff[ww]-mu.ref)/sd.ref

      }
      wleft<-wind-wwidth*dtt/2
      wright<-wind+wwidth*dtt/2
    } #close if statement for date-formatted data

  } #close if statement for adaptive reference window

  return(data.frame(wleft=wleft,wright=wright,ddiff=ddiff,zz=zz))
}

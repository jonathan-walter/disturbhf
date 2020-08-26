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
mwdistdiffz<-function(testy, refy, wwidth, refwidth=NULL, dx=0.01, stride=1, dmin=0.5){

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
  # Using ALL of refy as the reference period (no seasonal adaptive referece)

  if(is.null(refwidth)){
    ## TODO: Reconcile this section to make indexing and syntax as close as possible to bottom
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
        else if(any(diff(diff(refy$tt[pd]))>1e-6)){
          #returns NA if any times are skipped
          ddiff0[ww]<-NA
        }
        else{
          tmprefdist<-ecdf(refy$yy[-pd])
          tmpwdist<-ecdf(refy$yy[pd])
          ddiff0[ww]<-mean((tmprefdist(xx)-tmpwdist(xx))^2)
        }
      }
      rm(tmprefdist,tmpwdist)
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
          ddiff[ww]<-mean((refdist(xx)-wdist(xx))^2)
          zz[ww]<-(ddiff[ww]-mu.ref)/sd.ref
        }
      }
      wleft<-testy$tt[wind]-wwidth*dt/2
      wright<-testy$tt[wind]+wwidth*dt/2
    }

    if(any(grepl("POSIX",class(refy$tt)))){

      #Get mu and sd for excursions during the reference period
      wind0<-seq(from=ceiling(wwidth/2)+1, to=length(refy$tt)-ceiling(wwidth/2), by=stride)
      ddiff0<-rep(NA, length(wind0))
      for(ww in 1:length(wind0)){
        pd<-(wind0[ww]-(wwidth/2)):(wind0[ww]+wwidth/2-1) #get period for window
        if(mean(!is.na(refy$yy[pd])) < dmin){
          ddiff0[ww]<-NA
        }
        else if(any(diff(diff(refy$tt[pd]))>1e-6)){
          #returns NA if any times are skipped
          ddiff0[ww]<-NA
        }
        else{
          tmprefdist<-ecdf(refy$yy[-pd])
          tmpwdist<-ecdf(refy$yy[pd])
          ddiff0[ww]<-mean((tmprefdist(xx)-tmpwdist(xx))^2)
        }
      }
      rm(tmprefdist,tmpwdist)
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
          ddiff[ww]<-mean((refdist(xx)-wdist(xx))^2)
          zz[ww]<-(ddiff[ww]-mu.ref)/sd.ref
        }
      }
      wleft<-testy$tt[wind]-wwidth*dt/2*24*60*60
      wright<-testy$tt[wind]+wwidth*dt/2*24*60*60
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

      #Compute excursions in ref period and get mean and sd
      dt<-diff(testy$doy)[1]
      tmin<-min(refy$tt)
      tmax<-max(refy$tt)
      wind<-seq(from=tmin, to=tmax, by=stride*dt)
      ddiff<-rep(NA, length(wind))

      for(ww in 1:length(wind)){

        ltest<-wind[ww]-wwidth*dt/2 #left side of "test" window
        if(ltest < tmin){next} #skip indices where window overhangs beginning of time series
        rtest<-wind[ww]+wwidth*dt/2 #right side of "test" window
        if(rtest > tmax){break} #stop computation when window overhangs end of time series
        tpd<-refy$tt > ltest & refy$tt <= rtest

        lref<-refy$doy[abs(refy$tt-wind[ww]) < dt/10]-refwidth*dt/2 %% 365 #left side of reference window
        if(lref==0){lref==365}
        rref<-refy$doy[abs(refy$tt-wind[ww]) < dt/10]+refwidth*dt/2 %% 365 #right side of reference window
        if(rref==0){rref==365}
        rpd<-refy$doy > lref & refy$doy <= rref

        #check if sufficient non-missing values in reference and test periods
        if(mean(!is.na(refy$yy[tpd])) < dmin | mean(!is.na(refy$yy[rpd])) < dmin){
          ddiff[ww]<-NA
          next
        }
        #subref<-refy$yy[rpd]
        refdist<-ecdf(refy$yy[rpd])
        wdist<-ecdf(testy$yy[tpd])
        ddiff[ww]<-mean((refdist(xx)-wdist(xx))^2)
      }
      mu.ref<-mean(ddiff, na.rm=T)
      sd.ref<-sd(ddiff, na.rm=T)

      #Compute excursions in test period and get z-score
      dt<-diff(testy$doy)[1]
      tmin<-min(testy$tt)
      tmax<-max(testy$tt)
      wind<-seq(from=tmin, to=tmax, by=stride)
      ddiff<-rep(NA, length(wind))
      zz<-rep(NA, length(wind))

      for(ww in 1:length(wind)){

        ltest<-wind[ww]-wwidth*dt/2 #left side of test window
        if(ltest < tmin){next} #skip indices where window overhangs beginning of time series
        rtest<-wind[ww]+wwidth*dt/2 #right side of test window
        if(rtest > tmax){break} #stop computation when window overhands end of time series
        tpd<-testy$tt > ltest & testy$tt <= rtest

        lref<-testy$doy[abs(testy$tt-wind[ww]) < dt/10]-refwidth*dt/2 %% 365 #left side of reference window
        if(lref==0){lref==365}
        rref<-testy$doy[abs(testy$tt-wind[ww]) < dt/10]+refwidth*dt/2 %% 365 #right side of reference window
        if(rref==0){rref==365}
        rpd<-refy$doy > lref & refy$doy <= rref

        #check if sufficient non-missing values in reference and test periods
        if(mean(!is.na(testy$yy[tpd])) < dmin | mean(!is.na(refy$yy[rpd])) < dmin){
          ddiff[ww]<-NA
          next
        }

        refdist<-ecdf(refy$yy[rpd])
        wdist<-ecdf(testy$yy[tpd])
        ddiff[ww]<-mean((refdist(xx)-wdist(xx))^2)
        zz[ww]<-(ddiff[ww]-mu.ref)/sd.ref

      }
      wleft<-wind-wwidth*dt/2
      wright<-wind+wwidth*dt/2
    }

    #for tt as a date
    if(any(grepl("POSIX", class(testy$tt)))){

      #compute decimal day of year corresponding to ttt
      testy$doy <- decimal_doy(testy$tt)
      refy$doy <- decimal_doy(refy$tt)

      #Compute excursions in ref period and get mean and sd
      #many functions consider time in seconds, need to make sure to convert to days
      dt<-diff(testy$doy)[1]
      tmin<-min(refy$tt)
      tmax<-max(refy$tt)
      wind<-seq(from=tmin, to=tmax, by=stride*dt*24*60*60)
      ddiff<-rep(NA, length(wind))

      for(ww in 1:length(wind)){

        ltest<-wind[ww]-refwidth*dt/2*24*60*60 #left side of "test" window
        if(ltest < tmin){next} #skip indices where window overhangs beginning of time series
        rtest<-wind[ww]+refwidth*dt/2*24*60*60 #right side of "test" window
        if(rtest > tmax){break} #stop computation when window overhangs end of time series
        tpd<-refy$tt > ltest & refy$tt <= rtest

        lref<-testy$doy[testy$tt==wind[ww]]-refwidth*dt/2 %% 365 #left side of reference window
        if(lref==0){lref==365}
        rref<-testy$doy[testy$tt==wind[ww]]+refwidth*dt/2 %% 365 #right side of reference window
        if(rref==0){rref==365}
        rpd<-refy$doy > lref & refy$doy <= rref

        #check if sufficient non-missing values in reference and test periods
        if(mean(!is.na(refy$yy[tpd])) < dmin | mean(!is.na(refy$yy[rpd])) < dmin){
          ddiff0[ww]<-NA
          next
        }
        refdist<-ecdf(refy$yy[rpd])
        wdist<-ecdf(testy$yy[tpd])
        ddiff[ww]<-mean((refdist(xx)-wdist(xx))^2)
      }
      mu.ref<-mean(ddiff, na.rm=T)
      sd.ref<-sd(ddiff, na.rm=T)

      #Compute excursions in test period and get z-score
      dt<-diff(testy$doy)[1]
      tmin<-min(testy$tt)
      tmax<-max(testy$tt)
      wind<-seq(from=tmin, to=tmax, by=stride*dt*24*60*60)
      ddiff<-rep(NA, length(wind))
      zz<-rep(NA, length(wind))

      for(ww in 1:length(wind)){

        ltest<-wind[ww]-refwidth*dt/2*24*60*60 #left side of test window
        if(ltest < tmin){next} #skip indices where window overhangs beginning of time series
        rtest<-wind[ww]+refwidth*dt/2*24*60*60 #right side of test window
        if(rtest > tmax){break} #stop computation when window overhands end of time series
        tpd<-testy$tt > ltest & testy$tt <= rtest

        lref<-testy$doy[testy$tt==wind[ww]]-refwidth*dt/2 %% 365 #left side of reference window
        if(lref==0){lref==365}
        rref<-testy$doy[testy$tt==wind[ww]]+refwidth*dt/2 %% 365 #right side of reference window
        if(rref==0){rref==365}
        rpd<-refy$doy > lref & refy$doy <= rref

        #check if sufficient non-missing values in reference and test periods
        if(mean(!is.na(testy$yy[tpd])) < dmin | mean(!is.na(refy$yy[rpd])) < dmin){
          ddiff0[ww]<-NA
          next
        }

        refdist<-ecdf(refy$yy[rpd])
        wdist<-ecdf(testy$yy[tpd])
        ddiff[ww]<-mean((refdist(xx)-wdist(xx))^2)
        zz[ww]<-(ddiff[ww]-mu.ref)/sd.ref

      }
      wleft<-wind-wwidth*dt/2*24*60*60
      wright<-wind+wwidth*dt/2*24*60*60
    }

  }
  return(data.frame(wleft=wleft,wright=wright,ddiff=ddiff,zz=zz))
}

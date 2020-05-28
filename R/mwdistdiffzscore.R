#' Calculate differences in continuous distribution functions between reference and moving window observations
#'
#' \code{mwdistdiffz} computes differences between the continuous distribution functions (cdf) for observations within
#' a moving window and a reference distribution, and the z-scores of differences relative to samples of the reference distribution.
#' It is used for identifying recovery times from disturbance in time series data.
#'
#' @param testy a data frame representing the "test" period, i.e., the period over which to search for disturbance and recovery.
#' It must contain the columns \code{tt} and \code{yy}.
#' \code{tt} is time in day-of-year format corresponding to observations of \code{yy}.
#' \code{yy} is a numeric time series.
#' @param refy a data frame representiong the refrence period that \code{testy} is compared to.
#' It must contain the columns \code{tt} and \code{yy}, as in \code{testy}.
#' @param wwidth the moving window width, in number of time steps
#' @param refwidth the width of the rolling reference window, in number of time steps; if \code{NULL} do not use rolling reference window
#' @param dx increment between values at which to evaluate differences between the cdf for \code{yy} and \code{refy}
#' @param stride number of time steps by which the moving window advances
#' @param stat the summary statistic describing differences in the cdf between \code{yy} and \code{refy}. See details.
#' @param dmin the fraction of data that must be present (i.e., non-NA) in test (and, if applicable, adaptive reference) moving windows to procede with computations.
#'
#' @return \code{mwdistdiff} returns a data frame containing the columns:
#' \item{wstart}{the time corresponding to the beginning (left edge) of the moving window}
#' \item{wend}{the time corresponding to the end (right edge) of the moving window}
#' \item{ddiff}{differences from the reference distritbution}
#' \item{zz}{Z-scores representing the strength of excursion from reference}
#'
#' @details The value suppplied to \code{refwidth} determines whether a rolling reference window will be used,
#' or if all of refy will be used as the reference period.
#' Currently accepted values for \code{stat} are "mean" and "sum".
#' \code{testy} must be a regularly spaced time series, but \code{refy} can be irregularly spaced; however, if
#' there is insufficient data then output will contain \code{NA}s.
#' \code{dmin} takes into account NAs in time series as well as reference windows overlapping the edges or gaps in the time series.
#'
#' @author Jonathan Walter, \email{jaw3es@@virginia.edu}
#'
#' @examples
#' #need to add some
#'
#' @export
#'

#TODO: write this so that the function can accept data data? this would make it such that testy can span multiple years.

mwdistdiffz<-function(testy, refy, wwidth, refwidth=NULL, dx=0.01, stride=1, stat="mean", dmin=0.5){

  #a little error handling
  if(!is.data.frame(testy) | !"tt" %in% colnames(testy) | !"yy" %in% colnames(testy)){
    stop("testy must be a data frame containing the columns 'tt' and 'yy'")
  }
  if(!is.data.frame(refy) | !"tt" %in% colnames(refy) | !"yy" %in% colnames(refy)){
    stop("refy must be a data frame containing the columns 'tt' and 'yy'")
  }
  if(!stat %in% c("mean","sum")){
    stop("available values of 'stat' are mean and sum")
  }
  #TODO: Check that testy is regularly spaced
  #warnings
  if(nrow(testy)<50){
    warning("the test period is very short, results may not be robust")
  }
  if(nrow(refy)<50){
    warning("the referece period is very short, results may not be robust")
  }

  #Begin function
  xx<-seq(min(c(refy$yy,testy$yy),na.rm=T)-dx, max(c(refy$yy,testy$yy),na.rm=T)+dx, dx)

  #TODO: change referencing in this section to midpoint of moving window to be consistent with the adaptive window
  if(is.null(refwidth)){
    #get mean and SD of excursions based on reference period
    wind0<-seq(1, to=nrow(refy)-wwidth, by=stride)
    ddiff0<-rep(NA, length(wind0))

    for(ww in 1:length(wind0)){
      if(mean(!is.na(refy$yy[wind0[ww]:(wind0[ww]+wwidth)])) < dmin){
        ddiff0[ww]<-NA
      }
      else if(any(diff(diff(refy$tt[wind0[ww]:(wind0[ww]+wwidth)]))>1e-6)){
        ddiff0[ww]<-NA
      }
      else{
        tmprefdist<-ecdf(refy$yy[-(wind0[ww]:(wind0[ww]+wwidth))])
        tmpwdist<-ecdf(refy$yy[wind0[ww]:(wind0[ww]+wwidth)])
        if(stat=="sum"){
          ddiff0[ww]<-sum((tmprefdist(xx)-tmpwdist(xx))^2)
        }
        if(stat=="mean"){
          ddiff0[ww]<-mean((tmprefdist(xx)-tmpwdist(xx))^2)
        }
      }
    }
    rm(tmprefdist,tmpwdist)
    mu.ref<-mean(ddiff0, na.rm=T)
    sd.ref<-sd(ddiff0, na.rm=T)

    #measure excursions during the study period
    refdist<-ecdf(refy$yy) #make ecdf for reference period

    wind<-seq(from=1,to=length(testy$tt)-wwidth,by=stride)
    wstart<-refy$tt[wind]
    wend<-refy$tt[wind+wwidth]
    ddiff<-rep(NA, length(wind))
    zz<-rep(NA, length(wind))

    for(ww in 1:length(wind)){
      if(mean(!is.na(testy$yy[wind[ww]:(wind[ww]+wwidth)])) < dmin){
        ddiff[ww]<-NA
      }
      else{
        wdist<-ecdf(testy$yy[wind[ww]:(wind[ww]+wwidth)])
        if(stat=="sum"){
          ddiff[ww]<-sum((refdist(xx)-wdist(xx))^2)
        }
        if(stat=="mean"){
          ddiff[ww]<-mean((refdist(xx)-wdist(xx))^2)
        }
        zz[ww]<-(ddiff[ww]-mu.ref)/sd.ref
      }
    }
  }

  # This is with the adaptive window

  #TODO: Reference all subsetting based on tt!
  else{
    dt<-mean(diff(testy$tt))
    tmin<-min(testy$tt)
    tmax<-max(testy$tt)
    #get mean and SD of excursions based on reference period
    wind0<-seq(from=ceiling(tmin+refwidth*dt/2), to=floor(tmax-refwidth*dt/2-1), by=stride)
    #wind0<-seq(from=ceiling(refwidth/2), to=floor(length(refy)-(refwidth/2-1)), by=stride)
    ddiff0<-rep(NA, length(wind0))

    for(ww in 1:length(wind0)){
      # if(mean(!is.na(refy$yy[wind0[ww]:(wind0[ww]+wwidth)])) < dmin){
      #   ddiff0[ww]<-NA
      # }
      #return NA if too few observations
      if(mean(!is.na(refy$yy[refy$tt >= wind0[ww]-refwidth*dt/2 & refy$tt < wind0[ww]+refwidth*dt/2-1])) < dmin){
        ddiff0[ww]<-NA
      }
      else{
        tmpref<-refy$yy[refy$tt >= wind0[ww]-refwidth*dt/2 & refy$tt < wind0[ww]+refwidth*dt/2-1]
        if(length(tmpref)/refwidth < dmin){ #if insufficient data in reference window, return NA
          ddiff0[ww]<-NA
        }
        else{
          tmpref<-tmpref[-((length(tmpref)/2-wwidth/2):(length(tmpref)/2+wwidth/2-1))]
          tmprefdist<-ecdf(tmpref)
          tmpwdist<-ecdf(refy$yy[refy$tt >= wind0[ww]-refwidth*dt/2 & refy$tt < wind0[ww]+refwidth*dt/2-1])
          if(stat=="sum"){
            ddiff0[ww]<-sum((tmprefdist(xx)-tmpwdist(xx))^2)
          }
          if(stat=="mean"){
            ddiff0[ww]<-mean((tmprefdist(xx)-tmpwdist(xx))^2)
          }
        }
      }
    }
    rm(tmprefdist,tmpwdist)
    mu.ref<-mean(ddiff0, na.rm=T)
    sd.ref<-sd(ddiff0, na.rm=T)

    #measure excursions during the study period
    wind<-seq(from=ceiling(tmin+refwidth*dt/2), to=floor(tmax-refwidth*dt/2-1), by=stride)
    wstart<-testy$tt[wind]-refwidth*dt/2
    wend<-testy$tt[wind]+refwidth*dt/2 #changing referencing so wind gives the center of the moving window
    ddiff<-rep(NA, length(wind))
    zz<-rep(NA, length(wind))

    for(ww in 1:length(wind)){
      if(mean(!is.na(testy$yy[testy$tt >= wind0[ww]-wwidth*dt/2 & testy$tt < wind0[ww]+wwidth*dt/2-1])) < dmin
         | mean(!is.na(refy$yy[refy$tt >= wind0[ww]-refwidth*dt/2 & refy$tt < wind0[ww]+refwidth*dt/2-1])) < dmin){
        ddiff[ww]<-NA
      }
      else{
        if(wind0[ww]-refwidth*dt/2 < 1){
          tmpref<-refy$yy[refy$tt >= wind0[ww]-refwidth*dt/2 %% 365 + 1 & testy$tt < wind0[ww]+refwidth*dt/2-1]
        }
        else{
          tmpref<-refy$yy[refy$tt >= wind0[ww]-refwidth*dt/2 & testy$tt < wind0[ww]+refwidth*dt/2-1]
        }
        tmpref<-tmpref[-((length(tmpref)/2-wwidth/2):(length(tmpref)/2+wwidth/2-1))]
        refdist<-ecdf(tmpref)
        wdist<-ecdf(testy$yy[testy$tt >= wind0[ww]-wwidth*dt/2 & testy$tt < wind0[ww]+wwidth*dt/2-1])
        if(stat=="sum"){
          ddiff[ww]<-sum((refdist(xx)-wdist(xx))^2)
        }
        if(stat=="mean"){
          ddiff[ww]<-mean((refdist(xx)-wdist(xx))^2)
        }
        zz[ww]<-(ddiff[ww]-mu.ref)/sd.ref
      }
    }
  }

  return(data.frame(wstart=wstart,wend=wend,ddiff=ddiff,zz=zz))
}

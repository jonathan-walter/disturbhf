#' Sound an alarm for disturbances and determine recovery time
#'
#' \code{recovalarm} interprets output from \code{mwdistdiffz} to detect disturbances and determine recovery time
#'
#' @param mwdistdiffz.obj output from \code{mwdistdiffz}
#' @param dthresh threshold z-score corresponding to a disturbance alarm
#' @param rthresh threshold z-score corresponding to recovery from disturbance
#'
#' @return \code{recovalarm} returns a data frame containing the columns:
#' \item{dist.date}{time a disturbance was detected}
#' \item{recov.date}{time a recovery was detected}
#' \item{tdiff}{difference between \code{recov} and \code{disturb}; the recovery time}
#' \item{peakz}{peak z-score recorded during disturbance}
#' \item{peak.date}{date peak z-score was recorded}
#'
#' @details Values of \code{disturb} and \code{recov} are the times corresponding to the middle of the moving window.
#'
#' @author Jonathan Walter, \email{jaw3es@@virginia.edu}
#'
#' @examples
#' #need to add some
#'
#' @export
#'

#TODO: Add computations for peakz and peak.date

recovalarm<-function(mwdistdiffz.obj, dthresh=c(2), rthresh=c(0.5)){

  #some basic error handling
  if(!all(c("zz","wstart","wend") %in% colnames(mwdistdiffz.obj))){
    stop("mwdistdiffzz.obj must be output from function mwdistdiffz")
  }

  #helper function to get midpoints of moving windows
  med.date<-function(d){
    #by design d[1] is the start and d[2] is the end
    return(d[1] + floor((d[2]-d[1])/2))
  }

  zz<-mwdistdiffz.obj$zz
  wstart<-mwdistdiffz.obj$wstart
  wend<-mwdistdiffz.obj$wend

  #are there even disturbances?
  if(max(zz)<dthresh){
    out<-data.frame(dist.date=NA
                    ,recov.date=NA
                    ,tdiff=NA
                    ,peakz=NA
                    ,peak.date=NA)
    print("no disturbances detected with provided values to dthresh")
  }
  else{
    # zz.disturb<-zz < min(dthresh) | zz > max(dthresh)
    zz.recov<-zz < rthresh
    rle.disturb<-rle(zz > dthresh)
    dstarts<-cumsum(rle.disturb$lengths)[rle.disturb$values]-rle.disturb$lengths[rle.disturb$values]+1
    rle.recov<-rle(zz > min(rthresh) & zz < max(rthresh))
    #rstarts<-cumsum(rle.recov$lengths)[rle.recov$values]
    recov.ind<-rep(NA, length(dstarts))
    for(ii in 1:length(dstarts)){
       recov.ind[ii]<-min(which(zz.recov & 1:length(zz.recov)>dstarts[ii]))
    }

    if("POSIXct" %in% class(wstart)){
      out<-data.frame(
        dist.date=apply(cbind(wstart[dstarts],wend[dstarts]),1,med.date),
        recov.date=apply(cbind(wstart[recov.ind],wend[recov.ind]),1,med.date)
      )
    }
    else if(is.numeric(wstart)){
      out<-data.frame(
        dist.date=apply(cbind(wstart[dstarts],wend[dstarts]),1,median),
        recov.date=apply(cbind(wstart[recov.ind],wend[recov.ind]),1,median)
      )
    }
    else{stop("wstart and wend must be numeric or POSIXct (date)")}

  }
  out<-out[!duplicated(out$recov.date),]
  out$tdiff<-out$recov.date-out$dist.date

  peakz<-rep(NA, nrow(out))
  peak.date<-rep(NA, nrow(out))

  return(out)
}

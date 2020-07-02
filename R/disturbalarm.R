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

disturbalarm<-function(mwdistdiffz.obj, dthresh=2, rthresh=0.5){

  #some basic error handling
  if(!all(c("zz","wleft","wright") %in% colnames(mwdistdiffz.obj))){
    stop("mwdistdiffzz.obj must be output from function mwdistdiffz")
  }

  #helper function to get midpoints of moving windows
  med.date<-function(d){
    #by design d[1] is the start and d[2] is the end
    return(d[1] + floor((d[2]-d[1])/2))
  }

  zz<-mwdistdiffz.obj$zz
  wleft<-mwdistdiffz.obj$wleft
  wright<-mwdistdiffz.obj$wright

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
    zz.recov<-zz < rthresh
    rle.disturb<-rle(zz > dthresh)
    dstarts<-cumsum(rle.disturb$lengths)[rle.disturb$values]-rle.disturb$lengths[rle.disturb$values]+1
    rle.recov<-rle(zz > min(rthresh) & zz < max(rthresh))
    recov.ind<-rep(NA, length(dstarts))
    for(ii in 1:length(dstarts)){
       recov.ind[ii]<-min(which(zz.recov & 1:length(zz.recov)>dstarts[ii]))
    }

    if("POSIXct" %in% class(wleft)){
      out<-data.frame(
        dist.date=apply(cbind(wleft[dstarts],wright[dstarts]),1,med.date),
        recov.date=apply(cbind(wleft[recov.ind],wright[recov.ind]),1,med.date)
      )
    }
    else if(is.numeric(wleft)){
      out<-data.frame(
        dist.date=apply(cbind(wleft[dstarts],wright[dstarts]),1,median),
        recov.date=apply(cbind(wleft[recov.ind],wright[recov.ind]),1,median)
      )
    }
    else{stop("wleft and wright must be numeric or POSIXct (date)")}

  }
  out<-out[!duplicated(out$recov.date),]
  out$tdiff<-out$recov.date-out$dist.date

  peakz<-rep(NA, nrow(out))
  peak.date<-rep(NA, nrow(out))
  wmid<-apply(cbind(wleft,wright),1,median)
  dt<-diff(wleft)[1]

  for(ii in 1:nrow(out)){
    tmp<-zz[wmid>=out$dist.date[ii] & wmid<out$recov.date[ii]]
    peakz[ii]<-max(tmp)
    peak.date[ii]<-out$dist.date[ii] + which.max(tmp)*dt
  }

  out$peakz<-peakz
  out$peak.date<-peak.date
  return(out)
}

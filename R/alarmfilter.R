#' Filter alarm signals to reduce errors from transient events
#'
#' \code{alarmfilter} filters output from \code{recovalarm}
#'
#' @param recovalarm.obj output from \code{recovalarm}
#' @param dmin.dist a minimum duration for disturbances, in days
#' @param dmin.recov a minimum duration for recovery events, in days
#'
#' @return \code{alarmfilter} returns a data frame containing the columns:
#' \item{dist.date} time a disturbance was detected
#' \item{recov.date} time a recovery was detected
#' \item{tdiff} difference between \code{recov} and \code{disturb}; the recovery time.
#' This is the same format as output from \code{recovalarm}, but there will be fewer rows.
#'
#' @details \code{alarmfilter} merges disturbances separated by recoveries shorter than \code{dmin.recov}
#' and then removes disturbances shorter than \code{dmin.dist}
#'
#' @export

# TODO Check that this code works if dist.date and recov.date are POSIXct objects. And that it needs to if this is the case.

alarmfilter<-function(recovalarm.obj, dmin.dist, dmin.recov=dmin.dist){

  out<-recovalarm.obj
  if(nrow(out)<=1)return(out)
  else{

    #merge disturbances between transient recoveries
    for(rep in 1:1000){ #this is iterative because multiple merges could be possible
      if(nrow(out)<=1){break}
      drop1<-NULL
      for(ii in 1:(nrow(out)-1)){
        if(out$dist.date[ii+1] - out$recov.date[ii] < dmin.recov){
          out$recov.date[ii]<-out$recov.date[ii+1]
          out$tdiff[ii]<-out$recov.date[ii]-out$dist.date[ii]
          drop1<-c(drop1, ii)
        }
      }
      if(!is.null(drop1)){
        out<-out[-drop1,]
      }
      #break if there are no more that need to be merged
      if(!any(diff(drop1))==1){break}
    }

    #remove disturbances that are too short
    drop2<-NULL
    for(ii in 1:nrow(out)){
      if(!is.na(out$tdiff[ii])){
        if(out$tdiff[ii] < dmin.recov){
          drop2<-c(drop2,ii)
        }
      }
    }
    if(!is.null(drop2)){
      out<-out[-drop2,]
    }
    return(out)
  }

}


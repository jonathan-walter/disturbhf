#' Computes stats used to evaluate disturbance/recovery detection algorithm performance
#' 
#' @param true the \code{key} slot from \code{simDisturb}
#' @param estimated output from \code{recovalarm}
#' @param tol tolerance for error associated with the width of the moving window
#' 
#' @return a data frame with the number of rows = 1 (if no disturbance was detected) or equal to
#' the number of disturbances detected in the time series, and columns corresponding to the following:
#' \item{detect} 0 if a disturbance was detected inside the true disturbance window + tolerance
#' \item{delta.dist} the absolute difference from the day of true disturbance
#' \item{delta.recov} the absolute difference from the day of true recovery
#' \item{delta.tdiff} the absolute difference from the true recover time
#' 
#' @export
#' 
#' TODO: Fix this function to give a true positive if alarm is on when the disturbance hits

evalStats<-function(true,estimated,tol){
  
  if(all(is.na(estimated))){ #if there are no disturbances . . . 
    out<-data.frame(detect=0, delta.dist=NA, delta.recov=NA, delta.tdiff=NA)
  }
  else{
    out<-NULL
    for(ii in 1:nrow(estimated)){
      detect.ii<-ifelse(estimated$dist.date[ii]-tol <= true["dday"]-tol & 
                          estimated$recov.date[ii]+tol >= true["rday"],1,0)
      delta.dist.ii<-abs(true["dday"]-estimated$dist.date[ii])
      delta.recov.ii<-abs(true["rday"]-estimated$recov.date[ii])
      delta.tdiff.ii<-abs(true["recovtime"]-estimated$tdiff[ii])
      out<-rbind(out, c(detect.ii, delta.dist.ii, delta.recov.ii, delta.tdiff.ii))
    }
    out<-as.data.frame(out)
    colnames(out)<-c("detect","delta.dist","delta.recov","delta.tdiff")
  }
  
  return(out)
}
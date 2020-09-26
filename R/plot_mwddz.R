#' Plot output (& optinonally input) from mwdistdiffz()
#'
#' @param mwddz data frame, output from call to mwdistdiffz()
#' @param window_side options are "right"  (default), "left", or "center"; which time index of each window should be used
#' @param diff_ts options are "ddiff" or "zz"; should the actual difference between distributions be plotted ("ddiff"), or the z-score ("zz")
#' @param testy test data frame with columns 'tt' and 'yy' that went into mwdistdiffz()
#' @param refy reference data frame with columns 'tt' and 'yy' that went into mwdistdiffz()
#'
#' @return a plot of the difference between test and reference distribution through time
#' @export
#'
#' @examples
plot_mwddz <- function(mwddz, window_side = "right", diff_ts="ddiff", testy=NULL, refy=NULL){
  if(window_side == "center"){
    mwddz$wcenter = (mwddz$wleft + mwddz$wright) / 2
    timeCol = "wcenter"
  }else{
    timeCol = paste0("w", window_side)
  }
  if(!is.null(testy)){
    par(mfrow=c(2,1))
    Xrange = range(c(testy$tt, refy$tt), na.rm=TRUE)
    Yrange = range(c(testy$yy, refy$yy), na.rm=TRUE)
    plot(refy[, c("tt", "yy")], xlim=Xrange, ylim=Yrange, pch=16)
    points(testy[, c("tt", "yy")], col="blue", pch=16)
    legend(x="topright", legend=c("ref", "test"), col=c("black", "blue"), pch=16)
  }
  plot(mwddz[, c(timeCol, diff_ts)], type="l", lwd=3, col="darkgreen")
}

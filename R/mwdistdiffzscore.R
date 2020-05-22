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
#' @param refwidth the width of the rolling reference window; if \code{NULL} do not use rolling reference window
#' @param dx increment between values at which to evaluate differences between the cdf for \code{yy} and \code{refy}
#' @param stride number of time steps by which the moving window advances
#' @param stat the summary statistic describing differences in the cdf between \code{yy} and \code{refy}. See details.
#' 
#' @return \code{mwdistdiff} returns a data frame containing the columns:
#' \item{wstart} the time corresponding to the beginning (left edge) of the moving window
#' \item{wend} the time corresponding to the end (right edge) of the moving window
#' \item{ddiff} differences from the reference distritbution
#' \item{zz} Z-scores representing the strength of excursion from reference
#' 
#' @details The value suppplied to \code{refwidth} determines whether a rolling reference window will be used,
#' or if all of refy will be used as the reference period. Currently accepted values for \code{stat} are "mean" and "sum".
#' 
#' @author Jonathan Walter, \email{jaw3es@@virginia.edu}
#' 
#' @examples
#' #need to add some
#' 
#' @export

mwdistdiffz<-function(testy, refy, wwidth, refwidth=NULL, dx=0.01, stride=1, stat="mean"){
  
  #a little error handling
  if(!is.data.frame(testy) | !grepl("tt", colnames(testy)) | !grepl("yy", colnames(testy))){
    stop("testy must be a data frame containing the columns 'tt' and 'yy'")
  }
  if(!is.data.frame(refy) | !grepl("tt", colnames(refy)) | !grepl("yy", colnames(refy))){
    stop("refy must be a data frame containing the columns 'tt' and 'yy'")
  }
  if(!stat %in% c("mean","sum")){
    stop("available values of stat are mean and sum")
  }
  #warnings
  if(nrow(testy)<50){
    warning("the test period is very short, results may not be robust")
  }  
  if(nrow(refy)<50){
    warning("the referece period is very short, results may not be robust")
  }  

  xx<-seq(min(c(refy,yy),na.rm=T)-dx, max(c(refy,yy),na.rm=T)+dx, dx)
  
  
  #TODO: Reference all subsetting based on tt!
  
  
  #this uses all of refy as the reference period
  if(is.null(refwidth)){
    #get mean and SD of excursions based on reference period
    wind0<-seq(from=1,to=length(refy)-wwidth,by=stride)
    ddiff0<-rep(NA, length(wind0))
    
    for(ww in 1:length(wind0)){
      if(all(is.na(refy[wind0[ww]:(wind0[ww]+wwidth)]))){ #in future may want to set a minimum number of observations
        ddiff0[ww]<-NA
      }
      else{
        tmprefdist<-ecdf(refy[-(wind0[ww]:(wind0[ww]+wwidth))])
        tmpwdist<-ecdf(refy[wind0[ww]:(wind0[ww]+wwidth)])
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
    refdist<-ecdf(refy) #make ecdf for reference period
    xx<-seq(min(c(refy,yy),na.rm=T)-dx, max(c(refy,yy),na.rm=T)+dx, dx)
    
    wind<-seq(from=1,to=length(tt)-wwidth,by=stride)
    wstart<-tt[wind]
    wend<-tt[wind+wwidth]
    ddiff<-rep(NA, length(wind))
    zz<-rep(NA, length(wind))
    
    for(ww in 1:length(wind)){
      if(all(is.na(yy[wind[ww]:(wind[ww]+wwidth)]))){ #in future may want to set a minimum number of observations
        ddiff[ww]<-NA
      }
      else{
        wdist<-ecdf(yy[wind[ww]:(wind[ww]+wwidth)])
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
  else{
    #get mean and SD of excursions based on reference period
    wind0<-seq(from=ceiling(refwidth/2),to=floor(length(refy)-(refwidth/2-1)),by=stride)
    ddiff0<-rep(NA, length(wind0))
    
    for(ww in 1:length(wind0)){
      if(all(is.na(refy[wind0[ww]:(wind0[ww]+wwidth)]))){ #in future may want to set a minimum number of observations
        ddiff0[ww]<-NA
      }
      if(wind0[ww] < refwidth/2){
        ddiff0[wind0[ww]]<-NA
      }
      else{
        tmpref<-refy[(wind0[ww]-refwidth/2+1):(wind0[ww]+refwidth/2)]
        tmpref<-tmpref[-((length(tmpref)/2-wwidth/2):(length(tmpref)/2+wwidth/2-1))]
        tmprefdist<-ecdf(tmpref)#[-(wind0[ww]:(wind0[ww]+wwidth))])
        tmpwdist<-ecdf(refy[wind0[ww]:(wind0[ww]+wwidth-1)])
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
    wind<-seq(from=ceiling(refwidth/2),to=floor(length(refy)-(refwidth/2-1)),by=stride)
    wstart<-tt[wind]
    wend<-tt[wind+wwidth]
    ddiff<-rep(NA, length(wind))
    zz<-rep(NA, length(wind))
    
    for(ww in 1:length(wind)){
      if(all(is.na(yy[wind[ww]:(wind[ww]+wwidth)]))){ #in future may want to set a minimum number of observations
        ddiff[ww]<-NA
      }
      if(wind[ww] < refwidth/2){
        ddiff[ww]<-NA
      }
      else{
        tmpref<-refy[(wind[ww]-refwidth/2+1):(wind[ww]+refwidth/2)]
        tmpref<-tmpref[-((length(tmpref)/2-wwidth):(length(tmpref)/2+wwidth/2-1))]
        refdist<-ecdf(tmpref)#[-(wind0[ww]:(wind0[ww]+wwidth))])
        #refdist<-ecdf(refy[(wind[ww]-refwidth/2+1):(wind[ww]+refwidth/2)]) #need to add the reference distribution down here relative to window placement.
        
        wdist<-ecdf(yy[wind[ww]:(wind[ww]+wwidth-1)])
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
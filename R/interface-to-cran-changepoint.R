
## ************************************************************************
## 
## 
## 
## (c) Xiaobei Zhao
## 
## Wed Oct 07 09:59:21 EDT 2015 -0400 (Week 40)
## 
## 
## Reference: 
## 
## 
## ************************************************************************


##' .parse_segmfunc
##'
##' 
##' @title .parse_segmfunc
##' @param segtype 
##' @return function
##' @author Xiaobei Zhao
##' @examples
##' NGScopy:::.parse_segmfunc("mean.norm")
.parse_segmfunc <- function(segtype){
  segtype <- .parse_segmtype(segtype)
  segfunc <- switch(
    segtype,
    mean.norm=Xmisc::func(
      function(...) {
        dots <- list(...);
        dots["test.stat"] <- "Normal";
        dots["method"] <- "PELT"; #only allow for the default
        dots["penalty"] <- "SIC"; #only allow for the default    
        do.call(changepoint::cpt.mean,dots);
      },
      'cpt.mean','changepoint'),
    mean.cusum=Xmisc::func(
      function(...) {
        dots <- list(...);
        dots["test.stat"] <- "CUSUM";
        dots["method"] <- "BinSeg"; #only allow for the default
        dots["penalty"] <- "Asymptotic"; #only allow for the default        
        do.call(changepoint::cpt.mean,dots);
      },
      'cpt.mean','changepoint'),
    meanvar.norm=Xmisc::func(
      function(...) {
        dots <- list(...);
        dots["test.stat"] <- "Normal";
        dots["method"] <- "BinSeg"; #only allow for the default
        dots["penalty"] <- "Asymptotic"; #only allow for the default    
        do.call(changepoint::cpt.meanvar,dots);
      },
      'cpt.meanvar','changepoint'),
    var.css=Xmisc::func(
      function(...) {
        dots <- list(...);
        dots["test.stat"] <- "CSS";
        dots["method"] <- "BinSeg"; #only allow for the default
        dots["penalty"] <- "SIC"; #only allow for the default    
        do.call(changepoint::cpt.var,dots);
      },
      'cpt.var','changepoint')
    )
  invisible(segfunc)
}



##' .help_segmtype
##'
##' 
##' @title .help_segmtype
##' @param segtype 
##' @return NULL
##' @author Xiaobei Zhao
##' @examples
##' NGScopy:::.help_segmtype("mean.norm")
.help_segmtype <- function(segtype){
  segfunc <- .parse_segmfunc(segtype)
  message(lprintf('%(package)s::%(name)s',as.environment(attributes(segfunc))),'\n')
  help(attr(segfunc,'name'),attr(segfunc,'package'))
}

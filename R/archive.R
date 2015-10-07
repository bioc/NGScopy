
## ************************************************************************
## 
## 
## 
## (c) Xiaobei Zhao
## 
## Wed Oct 07 10:05:50 EDT 2015 -0400 (Week 40)
## 
## 
## Reference: 
## 
## 
## ************************************************************************



## .parse_segmfunc <- function(segtype){#changepoint_1.1.5
##   segtype <- .parse_segmtype(segtype)
##   segfunc <- switch(
##     segtype,
##     mean.norm=func(changepoint::multiple.mean.norm,'multiple.mean.norm','changepoint'),
##     meanvar.norm=func(changepoint::multiple.meanvar.norm,'multiple.meanvar.norm','changepoint'),
##     mean.cusum=func(changepoint::multiple.mean.cusum,'multiple.mean.cusum','changepoint'),
##     var.css=func(changepoint::multiple.var.css,'multiple.var.css','changepoint')
##     )
##   invisible(segfunc)
## }


## .help_segmtype <- function(segtype){#changepoint_1.1.5
##   segfunc <- .parse_segmfunc(segtype)
##   message(lprintf('%(package)s::%(name)s',as.environment(attributes(segfunc))),'\n')
##   help(attr(segfunc,'name'),attr(segfunc,'package'))
## }

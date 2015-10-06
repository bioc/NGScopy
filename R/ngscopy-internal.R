
## ************************************************************************
## 
## 
## 
## (c) Xiaobei Zhao
## 
## Thu May 22 10:36:58 EDT 2014 -0400 (Week 20)
## 
## 
## Reference: 
## require(rbamtools)
## 
## ************************************************************************


## ------------------------------------------------------------------------
## 
## ------------------------------------------------------------------------


.check.packages <- function(x)
{
  if (length(x)==1){
    if ( sprintf("package:%s",x) %in% search() ){
      return(TRUE)
    }
    ret <-
      tryCatch(
        {library(x,character.only=TRUE,logical.return=TRUE)},
        error=function(e){cat(as.character(e));return(FALSE)}
        )
    return(ret)
  } else{
    ret <- all(sapply(x,check.packages))
    return(ret)
  }
}


## ------------------------------------------------------------------------
## 
## ------------------------------------------------------------------------

.parse_scales <- function(scales){
  if (missing(scales)){
    scales <- c("fixed","free_y","free_x","free")
  } else{
    scales <- match.arg(scales,.parse_scales())
  }
  return(scales)
}

## ------------------------------------------------------------------------
## regions
## ------------------------------------------------------------------------

.check_regions <- function(x){
  if (!is.matrix(x) & !is.data.frame(x) & !is.character(x) & !is.file(x) & !is.connection(x)){
    stop('.check_regions | regions must be either matrix, data.frame, character, file or connection.')
  }
  
  if (is.character(x)) {
    if (x %in% c("","NULL")){
      x <- NULL
    }
  }

  if (!length(x)){
    return(data.frame(chr=character(),start=numeric(),end=numeric()))
  }

  if (is.matrix(x)){
    x <- as.data.frame(x)
  }

  if (is.data.frame(x)){
    if (ncol(x)!=3){
      stop('.check_regions | x must have three columns.')
    }
    .msg <- '.check_regions | x must be data.frame with three columns (chr, start, end).'
    tryCatch(
      {
        x[,1] <- as.character(x[,1]);
        x[,2] <- as.numeric(as.character(x[,2]));
        x[,3] <- as.numeric(as.character(x[,3]));
        colnames(x) <- c('chr','start','end')
      },
      warning=function(w){
        cat(as.character(w));
        stop(.msg)
      },
      error=function(e){
        ##cat(as.character(e));
        stop(.msg)
      }
      )
    if (any(is.na(x))){
      stop('.check_regions | x must not have missing values.')      
    }
    return(x)
  }
  
  if (is.file(x)){
    return(file(x))
  }
  
  if(is.character(x)){
    return(textConnection(x))
  }
}


.sort_regions <- function(x,.ref){
  'sort chromosomes by reference'
  stampmsg(".sort_regions | sort regions by reference")
  ret <- dfsort(x,1,levels=.ref_to_refid(.ref))
  return(ret)
}

.trim_regions <- function(x,.ref){
  'trim regsions if exceeding the reference'
  stampmsg(".trim_regions | trim regsions if exceeding the reference")
  ## ##logme(x)
  ## ##logme(.ref)
  .merge <- merge(x,.ref,by.x='refid',by.y='ID')
  .merge[.merge[,'start']<0,'start'] <- 0
  .merge[.merge[,'end']>.merge[,'LN'],'end'] <- .merge[.merge[,'end']>.merge[,'LN'],'LN']
  return(.merge[,c('refid','start','end')])  
}



## ------------------------------------------------------------------------
## segtype
## ------------------------------------------------------------------------


## segtype <- paste(.parse_segmtype(),sep='',collapse=',')
.check_segmtype <- function(segtype){
  if (!length(segtype)){
    return(character())
  }
  
  if (!is.character(segtype)){
    stop('.check_segmtype | segtype must be character.')
  }

  if (segtype %in% c("","NULL")){
    return(character())
  }
  
  if(length(grep(',',segtype))){
    segtype <- unlist(strsplit(segtype,','))
  }
  return(segtype)
}


.parse_segmtype <- function(segtype){
  if (missing(segtype)){
    segtype <- c("mean.norm","meanvar.norm","mean.cusum","var.css")
  } else{
    segtype <- match.arg(segtype,.parse_segmtype())
  }
  return(segtype)
}

.parse_segmfunc <- function(segtype){
  segtype <- .parse_segmtype(segtype)
  segfunc <- switch(
    segtype,
    mean.norm=func(changepoint:::multiple.mean.norm,'multiple.mean.norm','changepoint'),
    meanvar.norm=func(changepoint:::multiple.meanvar.norm,'multiple.meanvar.norm','changepoint'),
    mean.cusum=func(changepoint:::multiple.mean.cusum,'multiple.mean.cusum','changepoint'),
    var.css=func(changepoint:::multiple.var.css,'multiple.var.css','changepoint')
    )
  invisible(segfunc)
}


.help_segmtype <- function(segtype){
  segfunc <- .parse_segmfunc(segtype)
  message(lprintf('%(package)s::%(name)s',as.environment(attributes(segfunc))),'\n')
  help(attr(segfunc,'name'),attr(segfunc,'package'))
}

## ------------------------------------------------------------------------
## ref
## ------------------------------------------------------------------------


.refid_to_refname <- function(x,.ref,refname='chr'){
  'convert columns from refid (0-based) to refname'
  .refid <- x[,'refid']
  .refname <- .ref_to_refname(.ref)[as.character(.refid)]
  .which <- which(colnames(x)=='refid')
  x[,.which] <- .refname
  colnames(x)[.which] <- refname
  return(x)
}

.refname_to_refid <- function(x,.ref,refname='chr'){
  'convert columns from refname to refid'
  .refname <- x[,refname]
  .refid <- as.numeric(.ref_to_refid(.ref)[as.character(.refname)])
  .which <- which(colnames(x)==refname)
  x[,.which] <- .refid
  colnames(x)[.which] <- 'refid'
  return(x)
}

.parse_ref <- function(inFpath){
  reader <- rbamtools::bamReader(inFpath,idx=TRUE)
  .ref <- rbamtools::getRefData(reader)
  rbamtools::bamClose(reader)
  return(.ref)
}

.ref_to_refname <- function(.ref){
  structure(.ref[,2],names=.ref[,1])
}

.ref_to_refid <- function(.ref){
  structure(.ref[,1],names=.ref[,2])
}

.ref_to_reflength <- function(.ref){
  structure(.ref[,3],names=.ref[,1])
}

.ref_to_reflength2 <- function(.ref){
  structure(.ref[,3],names=.ref[,2])
}

.ref_to_regions <- function(.ref){
  .regions.tmp <- sapply(seq(nrow(.ref)),function(i) c(.ref[i,1],0,.ref[i,3]))
  .regions.mat <- t(.regions.tmp) # to matrix
  ## .regions.list <- split(.regions.tmp,rep(1:ncol(.regions.tmp),each=nrow(.regions.tmp))) # to list
  ret <- as.data.frame(.regions.mat)
  colnames(ret) <- c('refid','start','end')
  return(ret)
}



## ------------------------------------------------------------------------
## 
## ------------------------------------------------------------------------

## .bound_start <- function(x,coords){
##   min(x,coords[2])
## }

.bound_end <- function(x,coords){
  max(x,coords[3])
}


## .bound_start.v0.99.0 <- function(x,range){
##   as.numeric(min(x,rbamtools::getCoords(range)['begin']))
## }
## .bound_end.v0.99.0 <- function(x,range){
##   as.numeric(max(x,rbamtools::getCoords(range)['end']))
## }


##' .. content for \description{} (no empty lines) ..
##'
##' Make windows according to
##' 1) window's mindepth
##' 2) window's minsize
##' @title Make windows
##' @return matrix
##' @author Xiaobei Zhao
##' @param regions 
##' @param inFpath 
##' @param mindepth minimal number of reads overlapping the window
##' reads overlapping multiple windows only count once.
##' @param minsize minimal size of the window
##' @param ds numeric, downsampling factor
##' @param pcThreads
##' 
.make_windows <- function(
  regions=NULL,inFpath,mindepth,minsize,
  ds=1,
  pcThreads=1
  ){

  calc.ext <- FALSE ##TBD

  .unit <- function(reader,coords,inFpath,mindepth,minsize,calc.ext){
    stampmsg('(PID: ', Sys.getpid(), ') Processing coords (refid, start, end): ',vconcat(coords))
    ## ## require(rbamtools)
    ## reader <- rbamtools::bamReader(inFpath,idx=TRUE)
    rbamtools::rewind(reader)
    ##  [bamRange.initialize] coords must be numeric with length=3 (ref,start,stop)!
    ## one row of the region !
    if (length(coords)!=3 | !is.numeric(coords) | !is.vector(coords)){
      stop('.make_windows | Invalid format: coords')
    }
    range <- rbamtools::bamRange(reader,coords)
    ## The chr should be the same as identified by refid
    
    ## refname <- rbamtools::getRefName(range)
    ## rbamtools::rewind(range)

    if (calc.ext){
      ret <- data.frame(
        refid=numeric()
        ,start=numeric()
        ,end=numeric()
        ,size=numeric()
        ,depth=numeric()
        ,pos.mean=numeric()
        ,pos.median=numeric()
        ,pos.sd=numeric()
        ## ,depth.mean=numeric()
        ## ,depth.median=numeric()
        ## ,depth.sd=numeric()
        )
    } else {
      ret <- data.frame(
        refid=numeric()
        ,start=numeric()
        ,end=numeric()
        ,size=numeric()
        ,depth=numeric()
        ) 
    }
    ret <- as.matrix(ret)
    start <- coords[2] ##
    end <- 0 ## 
    depth <- 0
    start2 <- 0
    depth2 <- 0
    counter <- 0
    
    if (calc.ext){
      pos.v <- c()
      pos.v2 <- c()
      pos.v.pre <- c()
    }
    
    read <- NULL
    while(!is.null(read <- getNextAlign(range))){
      
      if (!counter%%ds) { ## downsampling
        refid <- rbamtools::refID(read)
        pos <- position(read) + 1 # convert to 1-based position
        ## stampmsg('(PID: ', Sys.getpid(), ') counter=',counter,', refid=',refid,', pos=',pos)
        ##
        start2 <- if (start==pos) start2 else start                   # store
        depth2 <- if (depth==0) depth2+1 else depth+1                 # store
        if (calc.ext) {
          pos.v2 <- if (!length(pos.v)) c(pos.v2,pos) else c(pos.v,pos) # store
        }
        ##
        depth <- depth+1
        if (calc.ext) {
          pos.v <- c(pos.v,pos)
        }
        if ( (depth>=mindepth & pos-start>=minsize) | pos==end ){ #satisfy rqmt.
          ## when pos==end, correct boundary effect
          end <- pos
          if (calc.ext) {
            pos.v.pre <- pos.v2
          }
          ##
          if (calc.ext) {
            .row <- cbind(
              refid=refid
              ,start=start2
              ,end=end
              ,size=pos-start2
              ,depth=depth2 ## Note: boundary effect
              ,pos.mean=mean(pos.v2)
              ,pos.median=median(pos.v2)
              ,pos.sd=sd(pos.v2)
              ##,depth.mean=depth2/(pos-start2)
              )
          } else {
            .row <- cbind(
              refid=refid
              ,start=start2
              ,end=end
              ,size=pos-start2
              ,depth=depth2 ## Note: boundary effect
              )
          }
          if (depth>=mindepth & pos-start>=minsize){
            ret <- rbind(ret,.row)
          } else {
            ret[nrow(ret),] <- .row
          }
          start <- pos
          depth <- 0
          if (calc.ext) {
            pos.v <- c()
          }
        }
      }
      counter <- counter+1
    }
    ## ##logme(ret)
    ## ##logme(coords)

    ## the remain (fail to meet rqmt.)
    
    .cnames <- colnames(ret)
    .which.start <- which(.cnames=="start")
    .which.end <- which(.cnames=="end")
    .which.size <- which(.cnames=="size")
    .which.depth <- which(.cnames=="depth")
    if (calc.ext) {
      .which.pos.mean <- which(.cnames=="pos.mean")
      .which.pos.median <- which(.cnames=="pos.median")
      .which.pos.sd <- which(.cnames=="pos.sd")
      ##.which.depth.mean <- which(.cnames=="depth.mean")
    }
    ##
    if (!is.null(read)){
      if (!nrow(ret)){
        ## save the remain to a new region if ret is empty
        if (calc.ext) {
          ret <- rbind(ret,cbind(
            refid=refid
            ,start=start
            ,end=pos
            ,size=pos-start
            ,depth=depth
            ,pos.mean=mean(pos.v)
            ,pos.median=median(pos.v)
            ,pos.sd=sd(pos.v)
            ##,depth.mean=depth/(pos-start)
            ))
        } else {
          ret <- rbind(ret,cbind(
            refid=refid
            ,start=start
            ,end=pos
            ,size=pos-start
            ,depth=depth
            ))
        }
      } else {
        ## merge the remain to the last region
        if (calc.ext) {
          .pos.v <- c(pos.v.pre,pos.v)
        }
        ret[nrow(ret),.which.depth] <- ret[nrow(ret),.which.depth]+depth
        if (calc.ext) {
          ret[nrow(ret),.which.pos.mean] <- mean(.pos.v)
          ret[nrow(ret),.which.pos.median] <- median(.pos.v)
          ret[nrow(ret),.which.pos.sd] <- sd(.pos.v)
          ##ret[nrow(ret),.which.depth.mean] <- ret[nrow(ret),.which.depth]/ret[nrow(ret),.which.size]
        }
      }
    }
    
    ## adjust end
    ## ## ret[nrow(ret),.which.end] <- coords[3]
    ret[nrow(ret),.which.end] <- .bound_end(ret[nrow(ret),.which.end],coords)
    ret[nrow(ret),.which.size] <- ret[nrow(ret),.which.end]-ret[nrow(ret),.which.start]
    
    ## rbamtools::bamClose(reader)
    return(ret)
  }
  
  if (pcThreads==1){
    ## ##logme("pcThreads==1")
    reader <- rbamtools::bamReader(inFpath,idx=TRUE)
    if (calc.ext) {
      ret <- data.frame(
        refid=numeric()
        ,start=numeric()
        ,end=numeric()
        ,size=numeric()
        ,depth=numeric()
        ,pos.mean=numeric()
        ,pos.median=numeric()
        ,pos.sd=numeric()
        ## ,depth.mean=numeric()
        )
    } else {
      ret <- data.frame(
        refid=numeric()
        ,start=numeric()
        ,end=numeric()
        ,size=numeric()
        ,depth=numeric()
        ) 
    }
    ret <- as.matrix(ret)
    for (i in seq(nrow(regions))){
      .coords <- as.numeric(regions[i,])
      ret <- rbind(ret,.unit(reader,.coords,inFpath=inFpath,mindepth=mindepth,minsize=minsize,calc.ext=calc.ext))
    }
    rbamtools::bamClose(reader)
    return(ret)
  } else {
    ## ##logme("pcThreads!=1")
    .regions <- dfsplit(regions,'refid',regions$refid)
    ##logme(.regions)
    ##logme(pcThreads)
    tmp <- parallel::mclapply(
      .regions,FUN=.make_windows,
      inFpath=inFpath, mindepth=mindepth,minsize=minsize,ds=ds,
      pcThreads=1,
      mc.cores=pcThreads)
    ret <- do.call(rbind,tmp)
    return(ret)
  }
}



##' .. content for \description{} (no empty lines) ..
##'
##' Count read starts in windows
##' reads overlapping multiple windows only count once.
##' @title Count read starts in windows
##' @param windows data.frame with three columns (refid,start,end)
##' @param inFpath 
##' @param ds numeric, downsampling factor
##' @param pcThreads 
##' @return numeric, the depths.
##' @author Xiaobei Zhao
.count_starts <- function(
  windows,inFpath,
  ds=1,
  pcThreads=1
  )
{
  calc.ext <- FALSE ##TBD
  
  if (pcThreads<=0){
    stop('.count_starts | pcThreads must be no less than 1.')
  }
  if (pcThreads==1){
    if (calc.ext){
      ret <- data.frame(
        depth=numeric()
        ,pos.mean=numeric()
        ,pos.median=numeric()
        ,pos.sd=numeric()
        )
    } else {
      ret <- data.frame(
        depth=numeric()
        )    
    }
    ret <- as.matrix(ret)
    if (nrow(windows)){ 
      reader <- rbamtools::bamReader(inFpath,idx=TRUE)
      ## ##logme(windows)
      counter=0
      .seq <- seq(nrow(windows))
      for (i in .seq){
        coords <- as.numeric(windows[i,])
        ## ##logme(coords)
        ## stampmsg('(PID: ', Sys.getpid(), ') Processing coords (refid, start, end): ',vconcat(coords))
        counter <- counter+1
        if (! counter %% 1000){
          stampmsg('(PID: ', Sys.getpid(), ') Processed ', counter,' windows.')
        }
        range <- rbamtools::bamRange(reader,coords)
        depth <- 0
        if (calc.ext){
          pos.v <- c()
        }
        counter <- 0
        while(!is.null(read <- rbamtools::getNextAlign(range))){
          if (!counter%%ds) { ## downsampling
            pos <- rbamtools::position(read) + 1 # convert to 1-based position
            if (pos>coords[2] & pos<=coords[3]){
              depth <- depth+1
              if (calc.ext){
                pos.v <- c(pos.v,pos)
              }
            }
          }
          counter <- counter+1
        }
        if (calc.ext){
          if (!depth){
            pos.mean <- pos.median <- pos.sd <- NA
          } else {
            pos.mean=mean(pos.v)
            pos.median=median(pos.v)
            pos.sd=sd(pos.v)
          }
          ret <- rbind(ret,cbind(
            depth=depth
            ,pos.mean=pos.mean
            ,pos.median=pos.median
            ,pos.sd=pos.sd
            ))
        } else {
          ret <- rbind(ret,cbind(
            depth=depth
            ))
        }
      }
      rbamtools::bamClose(reader)
    }
    return(ret)
  } else {
    .windows <- dfchunk(windows,n=pcThreads)    
    tmp <- parallel::mclapply(
      .windows,FUN=.count_starts,
      inFpath=inFpath,
      ds=ds,
      pcThreads=1,
      mc.cores=pcThreads)
    ## names(tmp) <- NULL
    ## ret <- do.call(c,tmp)
    ret <- do.call(rbind,tmp)
    return(ret)
  }
}


## ##' .. content for \description{} (no empty lines) ..
## ##'
## ##' Count reads in windows
## ##' reads overlapping multiple windows count multiple times.
## ##' @title Count reads in windows
## ##' @param inFpath 
## ##' @param windows 
## ##' @return numeric, the depths.
## ##' @author Xiaobei Zhao
## .count_reads <- function(
##   windows,inFpath,
##   pcThreads=1
##   )
## {
##   if (pcThreads<=0){
##     stop('.count_reads | pcThreads must be no less than 1.')
##   }
##   if (pcThreads==1){
##     reader <- rbamtools::bamReader(inFpath,idx=TRUE)
##     ret <- c()
##     for (i in seq(nrow(windows))){
##       coords <- as.numeric(windows[i,])
##       range <- rbamtools::bamRange(reader,coords)
##       ret <- c(ret,rbamtools::size(range))
##     }
##     rbamtools::bamClose(reader)
##     return(ret)
##   } else {
##     .windows <- dfchunk(windows,n=pcThreads)    
##     tmp <- parallel::mclapply(
##       .windows,FUN=.count_reads,
##       inFpath=inFpath,
##       mc.cores=pcThreads)
##     ret <- do.call(c,tmp)
##     return(ret)
##   }
## }


.calc_cn <- function(depthN,depthT,libsizeN,libsizeT,dsN,dsT,pseudocount=1,logr=TRUE){
  ##logme(pseudocount,'.calc_cn')
  ##logme(logr,'.calc_cn')
  cnr <- (depthT+pseudocount)/(depthN+pseudocount)/libsizeT*libsizeN*dsT/dsN ## normalize
  if (logr){
    cnr <- log2(cnr)
  }
  ret <- list(cnr=cnr,pseudocount=pseudocount,logr=logr)
  return(ret)
}



##' Compute segmentation
##'
##' 
##' @title Compute segmentation
##' @param data data.frame or matrix with columns of
##' 'chr','cnr'
##' @param segtype 
##' @param .ref 
##' @param .dots 
##' @param pcThreads 
##' @return segmants and plot
##' @author Xiaobei Zhao
.calc_segm <- function(
  data,segtype,.ref,.dots,
  pcThreads=1
  ){
  
  segtype <- match.arg(segtype,.parse_segmtype())
  segfunc <- .parse_segmfunc(segtype)
  ## ##logme(segtype)

  .data <- dfsplit(data,'chr',levels=.ref_to_refname(.ref))
  
  if (!length(.data)){
    stop('.calc_segm | chr must not be empty.')
  }

  .unit <- function(x,segfunc,.dots){
    chr <- x[,'chr'][1]
    ## pos <- x[,'pos']
    cn <- x[,'cnr']
    if (segtype=='mean.cusum'){
      if (!'penalty' %in% names(.dots)){
        .dots$penalty <- 'None'
      }
    }
    .dots$data <- cn
    ## logme(str(.dots))
    obj.cpt <- do.call(segfunc,.dots)
    return(obj.cpt)
  }

  tmp <- parallel::mclapply(
    .data,FUN=.unit,
    segfunc=segfunc,
    .dots=.dots,
    mc.cores=pcThreads)
  ret <- do.call(list,tmp)
  return(ret)
}

##' Process segmentation
##' 
##' Partially adapted from changepoint::plot
##' @title Process segmentation
##' @param obj.cpt list of cpt objects
##' @param data data.frame or matrix with columns of
##' 'chr','pos'
##' @param pcThreads 
##' @return matrix
##' @author Xiaobei Zhao
.proc_segm <- function(obj.cpt,data,pcThreads=1){
  
  ## 
  .unit <- function(.obj.cpt,chr,.data){
    ##.pos <- .data[.data[,'chr']==chr,'pos']
    .start <- .data[.data[,'chr']==chr,'start']
    .end <- .data[.data[,'chr']==chr,'end']
    
    ret <- data.frame(chr=character(),win.from=numeric(),win.to=numeric(),start=numeric(),end=numeric(),mean=numeric(),plottype=numeric())
    ret <- as.matrix(ret)    
    if(changepoint::cpttype(.obj.cpt)=="variance"){
      .plottype <- 2 # vertical segments
      ## .cpts <- changepoint::cpts(.obj.cpt) ## TBF(1.1.5) access directly
      .cpts <- .obj.cpt@cpts
      ncpt <- length(.cpts) # breaks, including end but not start
      .win.from <- NA
      .win.to <- seq_along(changepoint::data.set.ts(.obj.cpt))[.cpts]
      .start <- NA
      .end <- .end[.win.to]
      .means <- NA
    } else if(changepoint::cpttype(.obj.cpt)=="mean"  ||  changepoint::cpttype(.obj.cpt)=="mean and variance"){
      .plottype <- 1 # horizontal segments
      .cpts0 <- .obj.cpt@cpts# breaks, including end but not start
      ncpt <- length(.cpts0) # number of segments, not need to plus 1.
      .cpts <- c(0,.cpts0)   # breaks, including end and start
      if((changepoint::test.stat(.obj.cpt)=="Normal")||(changepoint::test.stat(.obj.cpt)=="CUSUM")){
        .means <- changepoint::param.est(.obj.cpt)$mean
      }else if(changepoint::test.stat(.obj.cpt)=="Gamma"){
        .means <- changepoint::param.est(.obj.cpt)$scale*changepoint::param.est(.obj.cpt)$shape
      }else if(changepoint::test.stat(.obj.cpt)=="Exponential"){
        .means <- 1/changepoint::param.est(.obj.cpt)$rate
      }else if(changepoint::test.stat(.obj.cpt)=="Poisson"){
        .means <- changepoint::param.est(.obj.cpt)$lambda
      }else{
        stop('Invalid Changepoint test statistic')
      }      
      .win.from <- sapply(seq(ncpt),function(i) {
        seq_along(changepoint::data.set.ts(.obj.cpt))[.cpts[i]+1]})
      .win.to <- sapply(seq(ncpt),function(i) {
        seq_along(changepoint::data.set.ts(.obj.cpt))[.cpts[i+1]]})
      ## cbind(.win.from,.win.to)
      .start <- .start[.win.from]
      .end <- .end[.win.to]
    } else{
      stop('Invalid Changepoint Type for plotting.\n Can only plot mean, variance, mean and variance')
    }
    ret <- cbind(chr,.win.from,.win.to,.start,.end,.means,.plottype)
    ##save(file='~/.obj.cpt.RData',.obj.cpt)
    ##logme(chr)
    ##logme(ret)
    colnames(ret) <- c('chr','win.from','win.to','start','end','mean','plottype')
    return(ret)
  } ## XB

  tmp <- parallel::mcmapply(
    FUN=.unit,
    obj.cpt,names(obj.cpt),
    MoreArgs=list(.data=data),
    SIMPLIFY=FALSE,
    mc.cores=pcThreads)
  ##logme(names(obj.cpt))
  ##logme(tmp)
  ret <- do.call(rbind,tmp)
  ret <- data.frame(ret,stringsAsFactors=FALSE)
  ## ret <- colclasses(ret,c("character",rep("numeric",ncol(ret))))
  for (i in seq(2,ncol(ret))){
    ret[,i] <- as.numeric(ret[,i])
  }
  return(ret)  
}



##' Unit plot segmentation
##'
##' 
##' @title Unit plot segmentation
##' @param data.segm the output of .proc_segm
##' @param ... 
##' @return NULL
##' @author Xiaobei Zhao
.unit_plot_segm <- function(chr,data.segm=NULL,...){
  if (!length(data.segm)){
    message('.unit_plot_segm | data.segm is empty')
    return()
  }
  .data.segm <- data.segm[data.segm[,'chr']==chr,]
  ncpt <- nrow(.data.segm)
  for (i in seq(ncpt)){
    if (.data.segm[i,'plottype']==1){
      graphics::segments(
        .data.segm[i,'start'],.data.segm[i,'mean'],.data.segm[i,'end'],.data.segm[i,'mean'],...
        )
    } else if (.data.segm[i,'plottype']==2){
      ## ##logme(.data.segm)
      if (i != ncpt){ # only internal breaks, i.e. excluding the end
        graphics::abline(v=.data.segm[i,'end'],...)
      }
    }
  }
}


.parse_mfrow <- function(.n){
  .mfrow <- switch(
    as.character(.n),
    "1"=c(1,1),
    "2"=c(1,2),
    c(ceiling(.n/2),2)
    )
  return(.mfrow)
}

.unit_plot_out <- function(
  chr,chrlength,
  data.cn,
  data.segm,
  xlim,ylim,
  xlab,ylab,
  ...,
  MoreArgs.plotcn,
  MoreArgs.plotseg,
  scales=c("fixed","free_y","free_x","free")
  )
{
  scales <- match.arg(scales)
  if(!length(xlab)) xlab='Positon'
  if(!length(ylab)) ylab='CNR'
  if(!length(MoreArgs.plotcn)) MoreArgs.plotcn=list(pch=20,col='#BEBEBE7F',cex=1)
  if(!length(MoreArgs.plotseg)) MoreArgs.plotseg=list(col='red',lwd=3,lty=1)
    
  .data.cn <- data.cn[data.cn$chr==chr,]

  if (! scales %in% c("free_y","free")){
    if (!length(ylim)) ylim <- range(data.cn[,'cnr'])
  }
  if (! scales %in% c("free_x","free")){
    if (!length(xlim)) xlim <- c(0,chrlength)
  }
  .n <- length(data.segm)
  if (!.n){
    .pos <- .data.cn[,'pos']
    .cn <- .data.cn[,'cnr']
    if (!length(ylim)) ylim <- range(.data.cn[,'cnr'])
    graphics::plot(xlim,ylim,xlim=xlim,ylim=ylim,type="n",main=chr,xlab=xlab,ylab=ylab,...)
    MoreArgs.plotcn$x <- .pos
    MoreArgs.plotcn$y <- .cn
    do.call(graphics::points,MoreArgs.plotcn)    
  }
  .mfrow <- .parse_mfrow(.n)
  par(mfrow=.mfrow)
  for (e in names(data.segm)){
    message('chr=',chr,', ','segtype=',e)
    .pos <- .data.cn[,'pos']
    .cn <- .data.cn[,'cnr']
    if (!length(ylim)) ylim <- range(.data.cn[,'cnr'])
    ##logme(xlim)
    ##logme(ylim)
    ##logme(chrlength)
    graphics::plot(xlim,ylim,xlim=xlim,ylim=ylim,type="n",main=sprintf('%s (%s)',chr,e),xlab=xlab,ylab=ylab,...)
    MoreArgs.plotcn$x <- .pos
    MoreArgs.plotcn$y <- .cn
    do.call(graphics::points,MoreArgs.plotcn)
    MoreArgs.plotseg$chr <- chr
    MoreArgs.plotseg$data.segm <- data.segm[[e]]
    do.call(.unit_plot_segm,MoreArgs.plotseg)
  }
}


.plot_out <- function(
  reflength,
  data.cn,
  data.segm,
  pdfFpath,
  width,height,
  xlim,ylim,
  xlab,ylab,
  ...,
  MoreArgs.plotcn,
  MoreArgs.plotseg,
  scales=c("fixed","free_y","free_x","free")
  ){
  if (missing(data.segm)) data.segm <- NULL
  if (missing(width)) width <- NULL
  if (missing(height)) height <- NULL
  if (missing(xlim)) xlim <- NULL
  if (missing(ylim)) ylim <- NULL
  if (missing(xlab)) xlab <- NULL
  if (missing(ylab)) ylab <- NULL
  if (missing(MoreArgs.plotcn)) MoreArgs.plotcn <- NULL
  if (missing(MoreArgs.plotseg)) MoreArgs.plotseg <- NULL
  
  .n <- length(data.segm)
  .mfrow <- .parse_mfrow(.n)
  if (!length(width)) width <- 7*.mfrow[2]
  if (!length(height)) height <- 7*.mfrow[1]

  ## proc_cn()
  chr_v <- unique(data.cn[,'chr'])
  pdf(file=pdfFpath,width=width,height=height)
  for (.chr in chr_v){
    .refid <- data.cn[data.cn[,'chr']==.chr,'refid'][1]
    .chrlength <- as.numeric(reflength[as.character(.refid)])
    ##logme(reflength)
    ##logme(.refid)
    ##logme(.chrlength)
    .unit_plot_out(
      chr=.chr,chrlength=.chrlength,
      data.cn=data.cn,
      data.segm=data.segm,
      xlim=xlim,ylim=ylim,
      xlab=xlab,ylab=ylab,
      ...,
      MoreArgs.plotcn=MoreArgs.plotcn,
      MoreArgs.plotseg=MoreArgs.plotseg,
      scales=scales
      )
  }
  dev.off()
  logsave(pdfFpath)
}




## ## ************************************************************************
## ## Below are
## ## external/internal utils
## ## some are adapted from Xmisc 
## ## ************************************************************************

vconcat <- function(x,sep=", ",capsule=FALSE,quote=FALSE){
  if (!length(x)){
    if (!capsule){
      return("")
    } else {
      return("c()")
    }
  }
  if (!is.vector(x)){
    stop('x must be a vector')
  }
  .flag <- is.character(x)
  if (quote & !.flag){
    warning('Surround non-character elements by double quotes. Try quote=FALSE.')
  }
  .x <- as.character(x)
  if (quote){
    .collapse <- sprintf('"%s"',sep)
  } else {
    .collapse <- sep
  }
  ret <- paste(.x,sep='',collapse=.collapse)
  if (quote){
    ret <- paste('"',ret,'"',sep='')    
  }
  if (capsule){
    ret <- sprintf('c(%s)',ret)
  }
  return(ret)
}


## ************************************************************************
## 
## ************************************************************************


## .norm_cn <- function(){ ##Caution: whether to adjust by median??
  
## }

## ##' @title RD/RC statistic 
## ##' @author Xiaobei Zhao
## .stat_rc <- function(){
  
## }



## ## ------------------------------------------------------------------------
## ## Below is relevant for single-sample input
## ## TBD
## ## ------------------------------------------------------------------------

## ##' @title correction by GC content
## ##' @author Xiaobei Zhao
## .corr_gc <- function(){
  
## }


## ##' @title correction by region mappability
## ##' @author Xiaobei Zhao
## .corr_rm <- function(){
  
## }

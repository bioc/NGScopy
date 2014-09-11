
## ************************************************************************
##
##
##
## (c) Xiaobei Zhao
##
## Tue May 20 16:01:39 EDT 2014 -0400 (Week 20)
##
##
## Reference:
## require(parallel)
##
##
## ************************************************************************



##' Detection of copy number variations in next generation sequencing (NGS)
##'
##' 
##' @title Detection of copy number variations in next generation sequencing (NGS)
##' @description
##' NGScopy is a reference class to detect copy number variations by ``restriction-imposed windowing'' in next generation sequencing (NGS).
##' 
##' @author Xiaobei Zhao
##' @field inFpathN character, 
##' the file path to the control (normal) sample
##' @field inFpathT character, 
##' the file path to the case (tumor) sample
##' @field outFpre character, 
##' the file path of the directory for output. 
##' @field libsizeN numeric, the library size of the control (normal) sample.
##' @field libsizeT numeric, the library size of the case (tumor) sample.
##' @field mindepth numeric, the minimal depth of reads per window.
##' @field minsize numeric, the minimal size of a window.
##' @field regions data.frame of three columns (chr/start/end),
##' the regions to study.
##' It follows the BED format: zero-based, half-open; (start,end].
##' @field segtype character, the type of change to capture during segmentation,
##' mean and/or variance, normal or nonparametric distributions.
##' A character vector with a single or multiple values from
##' c("mean.norm","meanvar.norm","mean.cusum","var.css"). 
##' see \code{changepoint}.
##' @field dsN integer, the downsampling factor of the control (normal) sample.
##' @field dsT integer, the downsampling factor of the case (tumor) sample.
##' @field MoreArgs.cn list, a list of arguments for method `calc_cn'.
##' See \code{set_MoreArgs.cn}.
##' @field MoreArgs.seg list, a list of arguments for method `calc_seg'.
##' See \code{set_MoreArgs.seg}.
##' @field pcThreads integer,
##' the number of processors performing the parallel computing.
##' @field auto.save logical, whether to save (any completed results)
##' automatically.
##' @field auto.load logical, whether to load (any previously completed results)
##' automatically.
##' @field force.rerun character, the names of methods to rerun regardless of
##' any previous runs, default to c().
##' @field out list, the output.
##' 
##' @examples
##' require(NGScopy)
##' require(NGScopyData)
##' 
##' ## Create an instance of `NGScopy' class
##' obj <- NGScopy$new(
##'   outFpre="ngscopy-case1",         # specified directory for output
##'   inFpathN=tps_N8.chr6()$bamFpath, # normal sample: tps_90.chr6.sort.bam
##'   inFpathT=tps_90.chr6()$bamFpath, # tumor sample:  tps_N8.chr6.sort.bam
##'   libsizeN=5777087,                # the library size of the normal sample
##'   libsizeT=4624267,                # the library size of the tumor sample
##'   mindepth=20,                     # the minimal depth of reads per window
##'   minsize=20000,                   # the minimal size of a window
##'   pcThreads=1                      # the number of processors for computing
##'   )
##' 
##' obj$show()                         # print the instance
##' 
##' \dontrun{
##'
##' ## Compute the copy number and save it
##' ## A data.frame will be saved to file `ngscopy_cn.txt' in the output directory
##' obj$write_cn()
##' 
##' ## Compute the segmentation and save it
##' ## A data.frame will be saved to file `ngscopy_seg.txt' in the output directory
##' obj$write_seg()
##' 
##' ## Visualization
##' ## A figure will be saved to file `ngscopy_out.pdf' in the output directory
##' obj$plot_cn()
##' }
##' 
##' @seealso \code{NGScopyData}
##' ##@exportClass NGScopy
NGScopy <- 
  setRefClass(
    'NGScopy',
    list(
      inFpathN='character',
      inFpathT='character',
      outFpre='character',
      libsizeN='numeric',
      libsizeT='numeric',
      mindepth='numeric',
      minsize='numeric',
      regions='data.frame',
      segtype='character',
      dsN='integer',
      dsT='integer',
      MoreArgs.cn='list',
      MoreArgs.seg='list',
      pcThreads='integer',
      auto.save='logical',
      auto.load='logical',
      force.rerun='character',
      out='list'
      ),
    contains='xRefClass',
    methods=list(
      initialize=function(...){
        .idx <- c(outFpre=1)
        callSuper(...,.index=.idx)
        setme()
      },
      setme=function(...){
        ##
        if (!length(auto.save)){
          auto.save <<- FALSE
          stampmsg(sprintf(
            'auto.save not specified, using a default of %s.',auto.save))
        }
        message('auto.save: ',.self$auto.save)
        
        ##
        if (!length(auto.load)){
          auto.load <<- FALSE
          stampmsg(sprintf(
            'auto.load not specified, using a default of %s.',auto.load))
        }
        message('auto.load: ',.self$auto.load)
        
        ## 
        if (length(inFpathN)){
          .set_inFpathN()
        }
        if (length(inFpathT)){
          .set_inFpathT()
        }
        if (length(outFpre)){
          .set_outFpre()
        }
        if (length(libsizeN)){
          .set_libsizeN()
        }
        if (length(libsizeT)){
          .set_libsizeT()
        }
        if (length(mindepth)){
          .set_mindepth()
        }
        if (length(minsize)){
          .set_minsize()
        }
        if (length(regions)){
          .set_regions()
        }
        if (length(segtype)){
          .set_segtype()
        }
        if (length(dsN)){
          .set_dsN()
        }
        if (length(dsT)){
          .set_dsT()
        }
        if (length(MoreArgs.cn)){
          .set_MoreArgs.cn()
        }
        if (length(MoreArgs.seg)){
          .set_MoreArgs.seg()
        }
        if (length(force.rerun)){
          .set_force.rerun()
        }
        if (length(pcThreads)){
          .set_pcThreads()
        }
        
        ##
        invisible()        
      },
      .set_inFpathN=function(){
        'For internal use only.'
        stampmsg("set_normal | inFpathN")
        if (!length(.self$inFpathN)){
          stop('inFpathN must be specified.')
        }
        if (!is.file(.self$inFpathN) | !is.file(paste(.self$inFpathN,'.bai',sep=''))){
          stop('inFpathN and its index file (.bai) must be valid.')          
        }
        out$inFpathN <<- .self$inFpathN
        ##
        message('inFpathN: ',.self$inFpathN)
        message('')
        invisible()
      },
      .set_mindepth=function(){
        'For internal use only.'
        stampmsg("set_normal | mindepth")
        if (!length(.self$mindepth)){
          message('mindepth not specified, using a default of 20.')
          .self$mindepth <- 20
        }
        out$mindepth <<- mindepth
        ##
        .scipen <- options()$scipen
        options(scipen=9)
        message('mindepth: ',.self$mindepth) 
        options(scipen=.scipen)    
        message('')   
        invisible()
      },
      .set_minsize=function(){
        'For internal use only.'
        stampmsg("set_normal | minsize")
        if (!length(.self$minsize)){
          message('minsize not specified, using a default of 20000.')
          .self$minsize <- 20000
        }
        out$minsize <<- minsize
        ##
        .scipen <- options()$scipen
        options(scipen=9)
        message('minsize: ',.self$minsize)
        options(scipen=.scipen)
        message('')   
        invisible()
      },
      set_inFpathN=function(inFpathN){
        'Set a control (normal) sample. inFpathN: The file path to the control (normal) sample.'
        if (!missing(inFpathN)){
          if(length(.self$inFpathN)){
            stop('set_inFpathN | inFpathN must be set only once.')
          }
          .self$inFpathN <- inFpathN
          .set_inFpathN()
        }
        invisible()        
      },
      set_mindepth=function(mindepth){
        'Set the minimal depth per window. mindepth: the minimal depth of reads per window in the control (normal) sample.'
        if (!missing(mindepth)) {
          if(length(.self$mindepth)){
            stop('set_mindepth | mindepth must be set only once.')
          } 
          .self$mindepth <- mindepth
          .set_mindepth()
        }
        invisible()
      },
      set_minsize=function(minsize){
        'Set the minimal size per window. minsize: the minimal size of a window in the control (normal) sample.'
        if (!missing(minsize)) {
          if(length(.self$minsize)){
            stop('set_minsize | minsize must be set only once.')
          }
          .self$minsize <- minsize
          .set_minsize()
        }
        invisible()
      },
      .set_normal=function(){
        'For internal use only.'
        .set_inFpathN(inFpathN)
        .set_mindepth(mindepth)
        .set_minsize(minsize)
        invisible()
      },
      set_normal=function(inFpathN,mindepth,minsize){
        'Set a control (normal) sample and minimal depth/size per window. See `set_inFpathN\', `set_mindepth\', `set_minsize\'.'
        set_inFpathN(inFpathN)
        set_mindepth(mindepth)
        set_minsize(minsize)
        invisible()
      },
      .set_inFpathT=function(){
        'For internal use only.'
        stampmsg("set_tumor | inFpathT")
        if (!length(.self$inFpathT)){
          stop('inFpathT must be specified.')
        }
        if (!is.file(.self$inFpathT) | !is.file(paste(.self$inFpathT,'.bai',sep=''))){
          stop('inFpathT and its index file (.bai) must be valid.')          
        }
        out$inFpathT <<- .self$inFpathT
        ##
        message('inFpathT: ',.self$inFpathT)
        message('')
        invisible()
      },
      set_inFpathT=function(inFpathT){
        'Set a case (tumor) sample. inFpathT: The file path to the case (tumor) sample.'
        if (!missing(inFpathT)){
          if(length(.self$inFpathT)){
            stop('set_inFpathT | inFpathT must be set only once.')
          } 
          .self$inFpathT <- inFpathT
          .set_inFpathT()
        }
        invisible()
      },
      set_tumor=function(inFpathT){
        'Set a case (tumor) sample. See `set_inFpathT\'.'
        set_inFpathT(inFpathT)
        invisible()
      },
      .set_outFpre=function(){
        'For internal use only.'
        stampmsg("set_outFpre")
        if (!length(.self$outFpre)){
          stop('set_outFpre | outFpre must be specified.')
        }
        ## out$outFpre <<- sprintf("%s.d%s_w%s",outFpre,mindepth,minsize)
        out$outFpre <<- .self$outFpre
        make.dir(out$outFpre)
        out$outFpath <<- sprintf("%s/ngscopy_out.RData",out$outFpre)
        ## out$objFpath <<- sprintf("%s/ngscopy_obj.RData",out$outFpre)
        message('outFpre: ',.self$outFpre)
        message('')
        invisible()        
      },
      set_outFpre=function(outFpre){
        'Set a directory for output. outFpre: the file path of the directory for output.'
        if (!missing(outFpre)) {
          if(length(.self$outFpre)){
            stop('set_outFpre | outFpre must be set only once.')
          }
          .self$outFpre <- outFpre
          .set_outFpre()
        }
        invisible()
      },
      .set_libsizeN=function(){
        'For internal use only.'
        stampmsg("set_libsize | libsizeN")
        if (!length(.self$libsizeN)){
          message('libsizeN not specified, assuming 10**6.')
          .self$libsizeN <- 10**6
        }
        out$libsizeN <<- libsizeN
        ##
        message('libsizeN: ',.self$libsizeN)
        message('')
        invisible()        
      },
      .set_libsizeT=function(){
        'For internal use only.'
        stampmsg("set_libsize | libsizeT")
        if (!length(.self$libsizeT)){
          message('libsizeT not specified, assuming 10**6.')
          .self$libsizeT <- 10**6
        }        
        out$libsizeT <<- libsizeT
        ##
        message('libsizeT: ',.self$libsizeT)
        message('')
        invisible()        
      },
      set_libsizeN=function(libsizeN){
        'Set library size of the control (normal). libsizeN: numeric, the library size of the control (normal) sample.'
        if (!missing(libsizeN)) {
          if(length(.self$libsizeN)){
            stop('set_libsizeN | libsizeN must be set only once.')
          } 
          .self$libsizeN <- libsizeN
          .set_libsizeN()
        }
        invisible()
      },
      set_libsizeT=function(libsizeT){
        'Set library size of the case (tumor). libsizeT: numeric, the library size of the case (tumor) sample.'
        if (!missing(libsizeT)) {
          if(length(.self$libsizeT)){
            stop('set_libsizeT | libsizeT must be set only once.')
          } 
          .self$libsizeT <- libsizeT
          .set_libsizeT()
        }
        invisible()
      },
      .set_libsize=function(){
        'For internal use only.'
        .set_libsizeN(libsizeN)
        .set_libsizeT(libsizeT)
        invisible()        
      },
      set_libsize=function(libsizeN,libsizeT){
        'Set library sizes. See `set_libsizeN\', `set_libsizeT\'.'
        set_libsizeN(libsizeN)
        set_libsizeT(libsizeT)
        invisible()
      },
      .set_regions=function(){
        'For internal use only.'
        stampmsg("set_regions")
        if(nrow(.self$regions)){
          if (! 'refid' %in% .self$regions){
            .refname <- unique(.self$regions[,'chr'])
            .refname0 <- unique(.ref_to_refname(.get_ref()))
            if (! all(.refname %in% .refname0) ){
              stop(sprintf(
                '.set_regions | inconsistency in chromosome names: the input (%s) should be all or a subset of the reference genome (%s).',
                vconcat(.refname,quote=TRUE),
                vconcat(.refname0,quote=TRUE)
                ))
            }
            .self$regions <- .refname_to_refid(.self$regions,.get_ref(),'chr')
          } 
        }else{
          .self$regions <- .ref_to_regions(.get_ref())
        }
        
        ## logme(.self$regions)
        .self$regions <-
          .trim_regions(.sort_regions(.self$regions,.get_ref()),.get_ref())
        out$regions <<- .self$regions
        ##
        .scipen <- options()$scipen
        options(scipen=9)
        message('regions: \n',dfconcat(.refid_to_refname(
          head(.self$regions),.get_ref(),refname='chr'
          )))
        options(scipen=.scipen)    
        if (nrow(.self$regions)>6){
          message('... ...')
        }
        message('')
      },
      set_regions=function(regions){
        'Set regions. regions: the regions in study, matrix, data.frame, character, file or connection with three columns (chr/start/end).'
        if(!missing(regions)){
          if (nrow(.self$regions)){
            stop('set_regions | regions must be set only once.')
          }
          .regions <- read_regions(regions)
          .self$regions <- .regions
          .set_regions()
        } 
      },
      .set_segtype=function(){
        'For internal use only.'
        stampmsg("set_segtype")
        .self$segtype <- .check_segtype(.self$segtype)
        if (!length(.self$segtype)){
          ##.default.segtype <- .parse_segtype()
          ##.default.segtype <- "mean.cusum"
          .default.segtype <- "mean.norm"
          message(sprintf(
            'segtype not specified, using a default of %s.',
            vconcat(.default.segtype,capsule=TRUE,quote=TRUE)))
          .self$segtype <- .default.segtype
        }
        out$segtype <<- .self$segtype
        ##
        message('segtype: ',vconcat(.self$segtype,capsule=TRUE,quote=TRUE))
        message('')
        invisible()
      },
      set_segtype=function(segtype){
        'Set the type of segmentation. segtype: a character vector with a single or multiple values from c("mean.norm","meanvar.norm","mean.cusum","var.css"). See `changepoint\'.'
        if (!missing(segtype)){
          if(length(.self$segtype)){
            message(sprintf(
              'set_segtype | segtype (%s) is overridden by (%s).',
              vconcat(.self$segtype,quote=TRUE),
              vconcat(segtype,quote=TRUE)
              ))
          }
          .self$segtype <- segtype
          .set_segtype()
        }
        invisible()
      },
      .set_dsN=function(){
        'For internal use only.'
        stampmsg("set_ds | dsN")
        if (!length(.self$dsN)){
          message('dsN not specified, assuming 1.')
          .self$dsN <- as.integer(1)
        }
        if (.self$dsN < 1) {
          stop('dsN must be an integer no less than 1.')
        }
        out$dsN <<- .self$dsN
        ##
        message('dsN: ',.self$dsN)
        message('')
        invisible()        
      },
      set_dsN=function(dsN){
        'Set downsampling factor of the control (normal). dsN: numeric, the library size of the control (normal) sample.'
        if (!missing(dsN)) {
          if(length(.self$dsN)){
            stop('set_dsN | dsN must be set only once.')
          }
          .self$dsN <- as.integer(dsN)
          .set_dsN()
        }
        invisible()
      },
      .set_dsT=function(){
        'For internal use only.'
        stampmsg("set_ds | dsT")
        if (!length(.self$dsT)){
          message('dsT not specified, assuming 1.')
          .self$dsT <- as.integer(1)
        }        
        if (.self$dsT < 1) {
          stop('dsT must be an integer no less than 1.')
        }
        out$dsT <<- .self$dsT
        ##
        message('dsT: ',.self$dsT)
        message('')
        invisible()        
      },
      set_dsT=function(dsT){
        'Set downsampling factor of the case (tumor). dsT: numeric, the library size of the case (tumor) sample.'
        if (!missing(dsT)) {
          if(length(.self$dsT)){
            stop('set_dsT | dsT must be set only once.')
          }
          .self$dsT <- as.integer(dsT)
          .set_dsT()
        }
        invisible()
      },
      set_ds=function(dsN,dsT){
        'Set downsampling factors. See `set_dsN\', `set_dsT\'.'
        set_dsN(dsN)
        set_dsT(dsT)
        invisible()
      },
      .set_MoreArgs.cn=function(){
        'For internal use only.'
        stampmsg("set_MoreArgs.cn")
        if (!length(.self$MoreArgs.cn)){
          message('MoreArgs.cn not specified, assuming list(pseudocount=1,logr=TRUE).')
          .self$MoreArgs.cn <- list(pseudocount=1,logr=TRUE)
        }
        out$MoreArgs.cn <<- MoreArgs.cn
        ##
        message(cat('MoreArgs.cn: '))
        message(str(.self$MoreArgs.cn))
        message('')
        invisible()        
      },
      set_MoreArgs.cn=function(...){
        'Set MoreArgs.cn. ..., (pseudocount: the pseudocounts added to the observed depth per window, default to 1; logr: logical, whether to make log2 transformation of the ratios, default to TRUE.)'
        .MoreArgs.cn <- list(...)
        if (length(.MoreArgs.cn)) {
          if(length(.self$MoreArgs.cn)){
            message('set_MoreArgs.cn | MoreArgs.cn is updated.')
          }
          ##.self$MoreArgs.cn <- .MoreArgs.cn
          ##.self$MoreArgs.cn <- utils::modifyList(.self$MoreArgs.cn,.MoreArgs.cn)
          .self$MoreArgs.cn <-
            utils::modifyList(list(pseudocount=1,logr=TRUE),.MoreArgs.cn)
          .set_MoreArgs.cn()
        }
        invisible()
      },
      .set_MoreArgs.seg=function(){
        'For internal use only.'
        stampmsg("set_MoreArgs.seg")
        if (!length(.self$MoreArgs.seg)){
          message('MoreArgs.seg not specified, assuming list().')
          .self$MoreArgs.seg <- list()
          ## .self$MoreArgs.seg <- list(maxnseg=10)
        }
        out$MoreArgs.seg <<- MoreArgs.seg
        ##
        message(cat('MoreArgs.seg: '))
        message(str(.self$MoreArgs.seg))
        message('')
        invisible()        
      },
      set_MoreArgs.seg=function(...){
        'Set MoreArgs.seg. ..., a list of other arguments to the funciton of segmentation given by segtype. See `help_segtype\'.'
        .MoreArgs.seg <- list(...)
        if (length(.MoreArgs.seg)) {
          if(length(.self$MoreArgs.seg)){
            message('set_MoreArgs.seg | MoreArgs.seg is updated.')
          }
          ##.self$MoreArgs.seg <- .MoreArgs.seg
          ##.self$MoreArgs.seg <- utils::modifyList(.self$MoreArgs.seg,.MoreArgs.seg)
          .self$MoreArgs.seg <- .MoreArgs.seg
          .set_MoreArgs.seg()
        }
        invisible()
      },
      .set_pcThreads=function(){
        'For internal use only.'
        stampmsg("set_pcThreads")
        if (!length(.self$pcThreads)){
          message('pcThreads not specified, assuming 1.')
          .self$pcThreads <- as.integer(1)
        }
        out$pcThreads <<- pcThreads
        ##
        message('pcThreads: ',.self$pcThreads)
        message('')
        invisible()        
      },
      set_pcThreads=function(pcThreads){
        'Set the number of processors. pcThreads: numeric, the number of processors performing the parallel computing. It should not exceed the system\'s capacity.'
        if (!missing(pcThreads)) {
          if(length(.self$pcThreads)){
            message(sprintf(
              'set_pcThreads | pcThreads (%s) overrides a previous value (%s).',
              pcThreads,
              .self$pcThreads
              ))
          }
          .self$pcThreads <- pcThreads
          .set_pcThreads()
        }
        invisible()
      },
      .set_force.rerun=function(){
        'For internal use only.'
        stampmsg("set_force.rerun")
        if (!length(.self$force.rerun)){
          message('force.rerun not specified, assuming c("").')
          .self$force.rerun <- c("")
        }
        out$force.rerun <<- force.rerun
        ##
        message('force.rerun: ',vconcat(.self$force.rerun,capsule=TRUE,quote=TRUE))
        message('')
        invisible()        
      },
      set_force.rerun=function(force.rerun){
        'Set force.rerun. Reset it with missing input.'
        if (!missing(force.rerun)) {
          .force.rerun <- force.rerun
          if (!length(.force.rerun)){
            .force.rerun <- c("")
          }
          if(length(.self$force.rerun)){
            message(sprintf(
              'set_force.rerun | force.rerun (%s) overrides a previous value (%s).',
              vconcat(.force.rerun,capsule=FALSE,quote=TRUE),
              vconcat(.self$force.rerun,capsule=FALSE,quote=TRUE)
              ))
          } 
        } else{
          .force.rerun <- c("")
        }
        .self$force.rerun <- .force.rerun
        .set_force.rerun()
        invisible()
      },
      get_inFpathN=function(){
        'Get inFpathN'
        if (!.loadme('inFpathN')){
          .set_inFpathN()
        }
        .self$inFpathN
      },
      get_inFpathT=function(){
        'Get inFpathT'
        if (!.loadme('inFpathT')){
          .set_inFpathT()
        }
        .self$inFpathT
      },
      get_outFpre=function(){
        'Get outFpre'
        if (!.loadme('outFpre')){
          .set_outFpre()
        }
        .self$outFpre
      },
      get_libsizeN=function(){
        'Get libsizeN'
        if (!.loadme('libsizeN')){
          .set_libsizeN()
        }
        .self$libsizeN
      },
      get_libsizeT=function(){
        'Get libsizeT'
        if (!.loadme('libsizeT')){
          .set_libsizeT()
        }
        .self$libsizeT
      },
      get_mindepth=function(){
        'Get mindepth'
        if (!.loadme('mindepth')){
          .set_mindepth()
        }
        .self$mindepth
      },
      get_minsize=function(){
        'Get minsize'
        if (!.loadme('minsize')){
          .set_minsize()
        }
        .self$minsize
      },
      .get_regions=function(){
        'For internal use only.'
        if (!.loadme('regions')){
          .set_regions()
        }
        .self$regions
      },
      get_regions=function(){
        'Get regions.'
        .refid_to_refname(.get_regions(),.get_ref(),'chr')
      },
      get_segtype=function(){
        'Get segtype, segmentation type(s).'
        if (!.loadme('segtype')){
          .set_segtype()
        }
        .self$segtype
      },
      get_dsN=function(){
        'Get dsN'
        if (!.loadme('dsN')){
          .set_dsN()
        }
        .self$dsN
      },
      get_dsT=function(){
        'Get dsT'
        if (!.loadme('dsT')){
          .set_dsT()
        }
        .self$dsT
      },
      get_MoreArgs.cn=function(){
        'Get MoreArgs.cn'
        if (!.loadme('MoreArgs.cn')){
          .set_MoreArgs.cn()
        }
        .self$MoreArgs.cn
      },
      get_MoreArgs.seg=function(){
        'Get MoreArgs.seg'
        if (!.loadme('MoreArgs.seg')){
          .set_MoreArgs.seg()
        }
        .self$MoreArgs.seg
      },
      get_pcThreads=function(){
        'Get pcThreads'
        if (!.loadme('pcThreads')){
          .set_pcThreads()
        }
        .self$pcThreads
      },
      .get_windows=function(){
        'For internal use only.'
        ## via proc_normal
        if (!.loadme('windows','proc_normal')){
          proc_normal()          
        }
        if (!'windows' %in% names(out)){
          stop('NGScopy$.get_windows | windows is not available.')
        }
        out[['windows']]
      },
      get_windows=function(){
        'Get the windows.'
        res <- .get_windows()
        res <- as.data.frame(res)
        res <- .refid_to_refname(res,.get_ref(),'chr')
        res
      },
      get_size=function(){
        'Get the size per window.'
        ## via proc_normal
        if (!.loadme('size','proc_normal')){
          proc_normal()
        }
        if (!'size' %in% names(out)){
          stop('NGScopy$get_size | size is not available.')
        }
        out[['size']]
      },
      get_pos=function(){
        'Get the position (midpoint) per window.'
        ## via proc_normal
        if (!.loadme('pos','proc_normal')){
          proc_normal()
        }
        if (!'pos' %in% names(out)){
          stop('NGScopy$get_pos | pos is not available.')
        }
        out[['pos']]
      },
      .get_dataN=function(){
        'For internal use only.'
        ## via proc_normal
        if (!.loadme('dataN','proc_normal')){
          proc_normal()
        }
        if (!'dataN' %in% names(out)){
          stop('NGScopy$get_dataN | dataN is not available.')
        }
        out[['dataN']]
      },
      get_dataN=function(){
        'Get the data of the normal per window.'
        .dataN <- .get_dataN()
        colnames(.dataN) <- paste(colnames(.dataN),'.N',sep='')
        colnames(.dataN)[colnames(.dataN)=='depth.N'] <- "depthN"
        return(.dataN)        
      },
      get_depthN=function(){
        'Get the depth of the normal per window.'
        get_dataN()[,'depthN']
      },
      .get_dataT=function(){
        'For internal use only.'
        ## via proc_tumor
        if (!.loadme('dataT','proc_tumor')){
          proc_tumor()
        }
        if (!'dataT' %in% names(out)){
          stop('NGScopy$get_dataT | dataT is not available.')
        }
        out[['dataT']]
      },
      get_dataT=function(){
        'Get the data of the normal per window.'
        .dataT <- .get_dataT()
        colnames(.dataT) <- paste(colnames(.dataT),'.T',sep='')
        colnames(.dataT)[colnames(.dataT)=='depth.T'] <- "depthT"
        return(.dataT)        
      },
      get_depthT=function(){
        'Get the depth of the tumor per window.'
        get_dataT()[,'depthT']
      },
      get_cn=function(){
        'Get the copy number object.'
        ## via calc_cn
        if (!.loadme('cn','calc_cn')){
          calc_cn()
        }
        if (!'cn' %in% names(out)){
          stop('NGScopy$get_cn | cn is not available.')
        }
        out[['cn']]
      },
      get_cnr=function(){
        'Get the copy number ratios.'
        get_cn()[['cnr']]
      },      
      .get_data.cn=function(){
        'For internal use only.'
        ## via proc_cn
        if (!.loadme('data.cn','proc_cn')){
          proc_cn()
        }
        out[['data.cn']]
      },
      get_data.cn=function(as.granges=FALSE){
        'Get the data.frame of copy number.'
        ret <- .get_data.cn()
        ret <- ret[,colnames(ret)!='refid']
        if (as.granges){
          .ret <- df.to.gr(ret,chrlength=get_reflength())
          if (length(.ret)){
            ret <- .ret
          } else{
            stop('get_data.cn | Please set as.granges to FALSE.')
          }
        }
        ret
      },
      get_seg=function(){
        'Get the segmentation object.'
        if (!.loadme('seg','calc_seg')){
          calc_seg()
        }
        if (!'seg' %in% names(out)){
          stop('NGScopy$get_seg | seg is not available.')
        }
        out[['seg']]
      },
      .get_data.seg=function(){
        'For internal use only.'
        ## via proc_seg
        if (!.loadme('data.seg','proc_seg')){
          proc_seg()
        }
        if (!'data.seg' %in% names(out)){
          stop('NGScopy$.get_data.seg | data.seg is not available.')
        }
        out[['data.seg']]
      },
      get_data.seg=function(as.granges=FALSE){
        'Get the data.frame of segmentation.'
        ret <- .get_data.seg()
        ret <- lapply(ret,function(e) {e[!names(e) %in% c('plottype')]})
        ret <- do.call(rbind,lapply(names(ret),function(a) {cbind(ret[[a]],segtype=a)}))
        ret$segtype <- as.character(ret$segtype)
        if (as.granges){
          .ret <- df.to.gr(ret,chrlength=get_reflength())
          if (length(.ret)){
            ret <- .ret
          } else{
            stop('get_data.seg | Please set as.granges to FALSE.')
          }
        }
        ret
      },
      .proc_ref=function(){
        'For internal use only.'
        ## Set reference genome of the normal sample.
        if (.loadme('ref')){
          return()
        }
        stampmsg(".proc_ref | process reference genome (in bam header of the normal sample)")
        if (!length(inFpathN)){
          stop('.proc_ref | inFpathN must be set.')
        }
        out$ref <<- .parse_ref(inFpathN)
        out$refname <<- .ref_to_refname(out$ref)
        out$reflength <<- .ref_to_reflength(out$ref)
        if(auto.save){
          saveme()
        }
        invisible()  
      },
      .get_ref=function(){
        'For internal use only.'
        'Get reference genome in the normal sample.'
        .proc_ref()
        out$ref
      },
      .get_refname=function(){
        'For internal use only.'
        .proc_ref()
        out$refname
      },
      get_refname=function(){
        'Get reference genome name in the normal sample.'
        as.character(.get_refname())
      },
      .get_reflength=function(){
        'For internal use only.'
        .proc_ref()
        out$reflength
      },
      get_reflength=function(){
        'Get reference genome length in the normal sample.'
        ret <- .get_reflength()
        names(ret) <- as.character(.get_refname()[names(ret)])
        return(ret)
      },
      saveme=function(){
        'Save the output for later usage.'
        ## .self$out$inFpathN <- .self$inFpathN
        ## .self$out$inFpathT <- .self$inFpathT
        ## .self$out$outFpre <- .self$outFpre
        ## .self$out$libsizeN <- .self$libsizeN
        ## .self$out$libsizeT <- .self$libsizeT
        ## .self$out$mindepth <- .self$mindepth
        ## .self$out$minsize <- .self$minsize
        ## .self$out$regions <- .self$regions
        ## .self$out$segtype <- .self$segtype
        ## .self$out$dsN <- .self$dsN
        ## .self$out$dsT <- .self$dsT
        ## .self$out$pcThreads <- .self$pcThreads
        ## .self$out$auto.save <- .self$auto.save
        ## .self$out$auto.load <- .self$auto.load
        base::save(file=out$outFpath,out)
      },
      ## saveobj=function(){
      ##   'Save the class instance for later usage.'
      ##   obj <- .self
      ##   base::save(file=out$objFpath,obj)
      ## },
      loadme=function(){
        'Load a previously saved output.'
        if (is.file(out$outFpath)) {
          base::load(file=out$outFpath)
          .names <- names(out)
          .names <- .names[-grep('^\\.',.names)]
          for (a in .names){
            .self$out[[a]] <- out[[a]]
          }
        }
      },
      .loadme=function(x,method=NULL){
        ret <- FALSE
        flag <- FALSE
        if (!length(method)){
          flag <- TRUE
        } else{
          if(!method %in% force.rerun){
            flag <- TRUE            
          }
        }
        if (flag){
          if (!x %in% names(out)){
            if (auto.load) {
              loadme()
            }          
          }
          if (x %in% names(out)){
            ret <- TRUE
          }
        }
        return(ret)
      },
      proc_normal=function(){
        'Process the normal sample: make the windows and count the reads per window.'
        
        if (.loadme('dataN','proc_normal')) return()
        ##
        stampmsg("proc_normal | this may take a while depending on the size of your library.")
        ##logme(regions)
        .regions <- .get_regions()
        tmp <- .make_windows(regions=.regions,inFpath=get_inFpathN(),mindepth=get_mindepth(),minsize=get_minsize(),ds=get_dsN(),pcThreads=get_pcThreads())
        ##logme(head(tmp))
        stampmsg("Processed all ",if (!length(nrow(.regions))) 0 else nrow(.regions)," regions.")
        if (!nrow(tmp)){
          stop('proc_normal | no enough reads to fulfil requirements.')
        }
        out$windows <<- tmp[,1:3]
        out$size <<- tmp[,4]
        out$pos <<- ceiling((tmp[,2]+tmp[,3])/2)
        out$dataN <<- tmp[,-(1:4),drop=FALSE]
        if(auto.save){
          saveme()
        }
      }, 
      proc_tumor=function(){
        'Process the tumor sample: count the reads per window.'
        
        if (.loadme('dataT','proc_tumor')) return()
        ##
        if (!'windows' %in% names(out)){
          proc_normal()
        }
        ##
        ## .set_tumor()
        stampmsg("proc_tumor | this may take a while depending on the size of your library.")
        .windows <- .get_windows()
        tmp <- .count_starts(.windows,get_inFpathT(),ds=get_dsT(),pcThreads=get_pcThreads())
        stampmsg("Processed all ",if (!length(nrow(.windows))) 0 else nrow(.windows)," windows.")
        if (!nrow(tmp)){
          stop('proc_tumor | no enough reads to fulfil requirements.')
        }
        ## XB
        out$dataT <<- tmp
        if(auto.save){
          saveme()
        }
      },
      calc_cn=function(){
        'Calculate the relative copy number ratios (CNRs) per window, using the depth in the control (normal) sample as denominator and the case (tumor) sample as numerator.'

        if (.loadme('cn','calc_cn')) return()
        ##
        stampmsg("calc_cn")
        .MoreArgs.cn <- get_MoreArgs.cn()
        .MoreArgs.cn$depthN <- get_depthN()
        .MoreArgs.cn$depthT <- get_depthT()
        .MoreArgs.cn$libsizeN <- get_libsizeN()
        .MoreArgs.cn$libsizeT <- get_libsizeT()
        .MoreArgs.cn$dsN <- get_dsN()
        .MoreArgs.cn$dsT <- get_dsT()
        out$cn <<- do.call(.calc_cn,.MoreArgs.cn)
        
        ## out$cn <<- .calc_cn(
        ##   get_depthN(),get_depthT(),
        ##   libsizeN=get_libsizeN(),libsizeT=get_libsizeT(),
        ##   dsN=get_dsN(),dsT=get_dsT(),
        ##   ...
        ##   )
        ##logme(summary(out$cn))
        
        if(auto.save){
          saveme()
        }
        stampmsg('calc_cn | done')
        message(cat('out$cn:'))
        message(str(out$cn))
        message('')
      },
      proc_cn=function(){
        'Process the output of coy number object and return as a data.frame.'
        
        if (.loadme('data.cn','proc_cn')) return()
        ##
        stampmsg("proc_cn")
        
        tmp <- cbind(.get_windows(),size=get_size(),pos=get_pos(),get_dataN(),get_dataT(),cnr=get_cnr())
        tmp <- as.data.frame.matrix(tmp)
        ## tmp <- .refid_to_refname(tmp,.get_ref(),'chr') # replace refid
        tmp[,'chr'] <- .ref_to_refname(.get_ref())[as.character(tmp[,'refid'])] # keep refid
        tmp <- tmp[,unique(c('refid','chr','start','end','size','pos','depthN','depthT','cnr',colnames(tmp)))]
        ## logme(head(tmp))
        out$data.cn <<- tmp
        if(auto.save){
          saveme()
        }
      },
      calc_seg=function(){
        'Calculate the segment and detect the change points.'

        if (.loadme('seg','calc_seg')) return()
        ##
        if (!'data.cn' %in% names(out)){
          proc_cn()
        }
        ##
        stampmsg("calc_seg")
        
        out$seg <<- list()
        for (e in get_segtype()){
          .MoreArgs.seg <- get_MoreArgs.seg()
          ## if (e %in% parse_segtype()){        
          ##   .MoreArgs.seg$Q <- .MoreArgs.seg$maxnseg
          ##   .MoreArgs.seg$maxnseg <- NULL
          ## }
          .seg <- .calc_seg(.get_data.cn(),segtype=e,.ref=.get_ref(),.dots=.MoreArgs.seg,pcThreads=get_pcThreads())
          out$seg[[e]] <<- .seg
        }
        if(auto.save){
          saveme()
        }
        stampmsg('calc_seg | done')
        message(cat('out$seg:'))
        message(str(out$seg))
        message('')
      },
      proc_seg=function(){
        'Process the output of segmentation object and return as a data.frame.'        

        if (.loadme('data.seg','proc_seg')) return()
        ##
        stampmsg("proc_seg")
        out$data.seg <<- list()
        for (e in get_segtype()){
          tmp <- .proc_seg(get_seg()[[e]],.get_data.cn(),pcThreads=get_pcThreads())
          out$data.seg[[e]] <<- tmp
        }
        if(auto.save){
          saveme()
        }
      },
      write_cn=function(cnFpath,...){
        'Write the output of copy numbers as a data.frame to a tab separated file. cnFpath: a file path (relative to outFpre) for `cn\' output; ...: see `Xmisc::write.data.table\'.'
        stampmsg("write_cn")
        out$.cnFpath <<- sprintf("%s/ngscopy_cn.txt",out$outFpre)
        if (missing(cnFpath)){
          cnFpath <- out$.cnFpath
        } else {
          cnFpath <- sprintf("%s/%s",out$outFpre,cnFpath)
        }
        tmp <- get_data.cn()
        write.data.table(cnFpath,tmp,...)
      },
      write_seg=function(segFpath,...){
        'Write the output of segments as a data.frame to a tab separated file. segFpath: a file path (relative to outFpre) for `seg\' output; ...: see `Xmisc::write.data.table\'.'
        stampmsg("write_seg")
        out$.segFpath <<- sprintf("%s/ngscopy_seg.txt",out$outFpre)
        if (missing(segFpath)){
          segFpath <- out$.segFpath
        } else {
          segFpath <- sprintf("%s/%s",out$outFpre,segFpath)
        }
        tmp <- get_data.seg()
        write.data.table(segFpath,tmp,...)
      },
      plot_cn=function(
        pdfFpath,
        width,height,
        scales,
        xlim,ylim,xlab,ylab,
        ...,
        MoreArgs.points,
        MoreArgs.segments
        ){
        'Plot the output and save to a pdf file. pdfFpath: a file path (relative to outFpre) for the pdf output; width,height: see `grDevices::pdf\'; scales: are scales shared across all chromossomes (i.e. x coordinates reflect the range of genomic coordinates per chromosome, y coordinates reflect the range of of CNRs per chromosome), given no specific `xlim\' and `ylim\' (the default, "fixed"), or do they vary across x coordinates ("free_x"), y coordinates ("free_y"), or both ("free"); xlim,ylim,xlab,ylab,...: see `graphics::plot\'; MoreArgs.points: additional arguments as in `graphics::points\'; MoreArgs.segments: additional arguments as in `graphics::segments\'.'
        stampmsg("plot_cn")
        
        if (missing(width)) width <- NULL
        if (missing(height)) height <- NULL
        if (missing(xlim)) xlim <- NULL
        if (missing(ylim)) ylim <- NULL
        if (missing(xlab)) xlab <- NULL
        if (missing(ylab)) ylab <- NULL
        if (missing(MoreArgs.points)) MoreArgs.points <- NULL
        if (missing(MoreArgs.segments)) MoreArgs.segments <- NULL
        
        scales <- .parse_scales(scales)[1]
        ## 
        reflength <- .get_reflength()
        data.cn <- .get_data.cn()
        data.seg <- .get_data.seg()
        out$.pdfFpath <<- sprintf("%s/ngscopy_out.pdf",out$outFpre)
        if (missing(pdfFpath)){
          pdfFpath <- out$.pdfFpath
        } else {
          pdfFpath <- sprintf("%s/%s",out$outFpre,pdfFpath)
        }
        ## 
        .plot_cn(
          reflength=reflength,
          data.cn=data.cn,
          data.seg=data.seg,
          pdfFpath=pdfFpath,
          width=width,
          height=height,
          xlim=xlim,
          ylim=ylim,
          xlab=xlab,
          ylab=ylab,
          ...,
          MoreArgs.points=MoreArgs.points,
          MoreArgs.segments=MoreArgs.segments,
          scales=scales
          )
      },
      save_normal=function(){
        'Get the output of the normal for later usage.'
        ## save using a reserved variable name 
        ##..ngscopy.normal.. <- out[c('windows','dataN','size')]
        out$.normalFpath <<- sprintf("%s/ngscopy_normal.RData",get_outFpre())
        ..ngscopy.normal.. <-
          list(
            inFpathN=get_inFpathN(),
            libsizeN=get_libsizeN(),
            mindepth=get_mindepth(),
            minsize=get_minsize(),
            regions=.get_regions(),
            segtype=get_segtype(),
            dsN=get_dsN(),
            windows=.get_windows(),
            size=get_size(),
            pos=get_pos(),
            depthN=get_depthN(),
            dataN=.get_dataN()
            ) ##XB
        base::save(file=out$.normalFpath,..ngscopy.normal..)
        logsave(out$.normalFpath)
      },
      load_normal=function(normalDpath){
        'Load a previously saved output of the normal. normalDpath: the path to the .RData file for the output of the normal.'
        out$.normalFpath <<- sprintf("%s/ngscopy_normal.RData",normalDpath)
        if (!is.file(out$.normalFpath)){
          stop('NGScopy$load_normal | the output of the normal is not available.')
        }
        base::load(file=out$.normalFpath)
        out$inFpathN <<- ..ngscopy.normal..$inFpathN
        out$libsizeN <<- ..ngscopy.normal..$libsizeN
        out$mindepth <<- ..ngscopy.normal..$mindepth
        out$minsize <<- ..ngscopy.normal..$minsize
        out$regions <<- ..ngscopy.normal..$regions
        out$segtype <<- ..ngscopy.normal..$segtype
        out$dsN <<- ..ngscopy.normal..$dsN
        out$windows <<- ..ngscopy.normal..$windows
        out$size <<- ..ngscopy.normal..$size
        out$pos <<- ..ngscopy.normal..$pos
        out$depthN <<- ..ngscopy.normal..$depthN
        out$dataN <<- ..ngscopy.normal..$dataN

        .self$inFpathN <- ..ngscopy.normal..$inFpathN
        .self$libsizeN <- ..ngscopy.normal..$libsizeN
        .self$mindepth <- ..ngscopy.normal..$mindepth
        .self$minsize <- ..ngscopy.normal..$minsize
        .self$regions <- ..ngscopy.normal..$regions
        .self$segtype <- ..ngscopy.normal..$segtype
        .self$dsN <- ..ngscopy.normal..$dsN        
      },
      show=function(){
        'Show the class instance friendly.'
        .scipen <- options()$scipen
        options(scipen=9)
        if (length(.self$inFpathN)) cat0('inFpathN: ',.self$inFpathN)
        if (length(.self$inFpathT)) cat0('inFpathT: ',.self$inFpathT)
        if (length(.self$outFpre)) cat0('outFpre: ',.self$outFpre)
        if (length(.self$libsizeN)) cat0('libsizeN: ',.self$libsizeN)
        if (length(.self$libsizeT)) cat0('libsizeT: ',.self$libsizeT)
        if (length(.self$mindepth)) cat0('mindepth: ',.self$mindepth) 
        if (length(.self$minsize)) cat0('minsize: ',.self$minsize)
        if (nrow(.self$regions)) {
          cat0('regions: \n',dfconcat(.refid_to_refname(
          head(.self$regions),.get_ref(),refname='chr'
          )))
          if (nrow(.self$regions)>6){
            cat0('... ...')
          }
        }
        if (length(.self$segtype)) cat0('segtype: ',vconcat(.self$segtype,capsule=TRUE,quote=TRUE))
        if (length(.self$dsN)) cat0('dsN: ',.self$dsN)
        if (length(.self$dsT)) cat0('dsT: ',.self$dsT)
        if (length(.self$MoreArgs.cn)){
          cat0('MoreArgs.cn:')
          cat0(str(.self$MoreArgs.cn))
        }
        if (length(.self$MoreArgs.seg)){
          cat0('MoreArgs.seg:')
          cat0(str(.self$MoreArgs.seg))
        }
        if (length(.self$auto.save)) cat0('auto.save: ',.self$auto.save)
        if (length(.self$auto.load)) cat0('auto.load: ',.self$auto.load)        
        options(scipen=.scipen)    
      }
      )
    )



##' Read regions from a data.frame, a file or a connection.
##'
##' 
##' @title Read regions from a data.frame, a file or a connection.
##' @param x data.frame or character, the regions in study.
##' A data.frame with three columns (chr/start/end), a file path or a
##' connection with three columns (chr/start/end) without column names.
##' Return an empty data.frame if "", "NULL" or NULL.
##' @return data.frame of three columns (chr/start/end)
##' @author Xiaobei Zhao
##' @examples
##' read_regions("
##' chr1 0 249250621
##' chr2 0 243199373
##' chr3 0 198022430
##' chr4 0 191154276
##' chr5 0 180915260
##' chr6 0 171115067
##' ")
read_regions <- function(x){
  x <- .check_regions(x)
  if (is.data.frame(x)){
    return(x)
  }
  ret <- read.table(
    x,header=FALSE,
    col.names=c('chr','start','end'),
    colClasses=c('character','numeric','numeric')
    )
  try(close(x),silent=TRUE)
  return(ret)
}


##' Parse the type of segmentation or return all available types
##' with missing input.
##' 
##' @title Parse the type of segmentation
##' @param segtype character, the type of segmentation.
##' Return all available types if missing.
##' @return the type of segmentation
##' @author Xiaobei Zhao
##' @examples
##' parse_segtype()
parse_segtype <- function(segtype){
  .parse_segtype(segtype)
}


##' Get help for segmentation functions given the type of
##' segmentation.
##'
##' 
##' @title Get help for segmentation functions
##' @param segtype character, the type of segmentation.
##' See \code{parse_segtype()}
##' @return NULL
##' @author Xiaobei Zhao
##' @examples
##' \dontrun{
##' help_segtype(parse_segtype()[1])
##' }
help_segtype <- function(segtype){
  .help_segtype(segtype)
}


##' Convert a data.frame to a GRanges object
##'
##' 
##' @title Convert a data.frame to a GRanges object
##' @param x data.frame or matrix
##' @param which.chr character, the column name of `chr`
##' @param which.start character, the column name of `start`
##' @param which.end character, the column name of `end`
##' @param which.width character, the column name of `width`
##' @param which.name character, the column name of `name`
##' @param chrlength numeric, the lengths of the chromosomes
##' @param start0 logical, wheter the `start` is 0-based.
##' @return a GRanges object
##' @author Xiaobei Zhao
##' @examples
##' \dontrun{
##' x <- data.frame(
##'   chr=c('chr1','chr2'),start=0,end=100,
##'   name=paste0('ID',1:2),score=1:2
##'   )
##' df.to.gr(x)
##' df.to.gr(x,chrlength=c(chr1=1000,chr2=1200))
##' }
df.to.gr <- function(
  x,
  which.chr='chr',which.start='start',which.end='end',
  which.width='width',which.name='name',
  chrlength=NULL,
  start0=TRUE
  ){
  if (!.check.packages("GenomicRanges")){
    message('df.to.gr | Package GenomicRanges is not installed.')
    return()
  }
  ## logme(search(),'df.to.gr | 1')
  x <- as.data.frame(x)
  n <- nrow(x)
  .x <- as.list(x)
  .x[[which.chr]] <- NULL
  .x[[which.start]] <- NULL
  .x[[which.end]] <- NULL
  .x[[which.width]] <- NULL
  .x[[which.name]] <- NULL
  if (length(x[[which.start]]) & start0){
    x[[which.start]] <- x[[which.start]] + 1
  }
  .x$seqnames <- Rle(x[[which.chr]],rep(1,n))
  .x$ranges <- IRanges(x[[which.start]],x[[which.end]],width=x[[which.width]],names=x[[which.name]])
  
  gr <- do.call(GRanges,.x)
  if(length(chrlength) & length(x[[which.chr]])){
    if(length(names(chrlength))){
      chrlength <- chrlength[names(chrlength) %in% as.character(x[[which.chr]])]
    }
    seqlengths(gr) <- chrlength
  }
  ## logme(search(),'df.to.gr | 2')
  return(gr)  
}


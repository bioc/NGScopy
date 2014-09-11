
## ************************************************************************
## 
## 
## 
## (c) Xiaobei Zhao
## 
## Wed Aug 06 14:01:26 EDT 2014 -0400 (Week 31)
## 
## 
## Reference: 
## 
## 
## ************************************************************************




##' Run an example of NGScopy at UNIX-like command line
##'
##' 
##' @title Run an example of NGScopy at UNIX-like command line
##' @param ifRun logical, whether to run
##' @param pcThreads numeric,
##' the number of processors performing the parallel computing.
##' @return character, the command line.
##' @author Xiaobei Zhao
##' @examples
##' ## To run at R prompt
##' ngscopy_cmdline_example()
##' \dontrun{
##' ngscopy_cmdline_example(ifRun=TRUE) # Note: this would take a while. And
##'                                     # R will create a folder for output
##'                                     # at current working directory
##' }
##' 
##' ## ## To run at Unix-like command line (not at R prompt)
##' ## Rscript -e "require(methods);NGScopy::ngscopy_cmdline_example(TRUE)"
ngscopy_cmdline_example <- function(ifRun=FALSE,pcThreads=1){
  if (!check.packages("NGScopyData")) {
    warning('ngscopy_cmdline_example',' | R package `NGScopyData` must be available to run the example.')
    return()
  } else {
    ngscopyCmd <- Xmisc::get_executable('NGScopy')
    inFpathN <- Xmisc::get_extdata('NGScopyData','tps_N8.chr6.sort.bam')
    inFpathT <- Xmisc::get_extdata('NGScopyData','tps_90.chr6.sort.bam')
    libsizeN <- 5777087
    libsizeT <- 4624267
    regions <- system.file('extdata','hg19_chr6_0_171115067.txt',package='NGScopy', mustWork=TRUE)
    dsN <- 1
    dsT <- 1
    outFpre <- 'ngscopy_cmdline_example'

    cmd <- lprintf('time Rscript %(ngscopyCmd)s --inFpathN=%(inFpathN)s --inFpathT=%(inFpathT)s --outFpre="%(outFpre)s" --libsizeN=%(libsizeN)s --libsizeT=%(libsizeT)s --regions=%(regions)s --dsN=%(dsN)s --dsT=%(dsT)s --pcThreads=%(pcThreads)s')

    ## ## ##
    ## ## sink('~/tmp.sink.out')
    ## ## logme(getwd())
    ## ## logme(ngscopyCmd)
    ## ## logme(cmd)
    ## ## sink()

    ##
    if(ifRun){
      message('Running a command line:\n',cmd)
      system(cmd)
    }
    invisible(cmd)
  }
}



##' A wrapper to run unit testing of NGScopy
##'
##' 
##' @title A wrapper to run unit testing of NGScopy
##' @param ... additional arguments to \code{UnitTest$new}
##' @return NULL
##' @author Xiaobei Zhao
##' @examples
##' ## To run at R prompt
##' \dontrun{
##' ngscopy_unittest() # Note: R will create a folder for output
##'                    # at current working directory
##' }
##' 
##' ## ## To run at Unix-like command line (not at R prompt)
##' ## Rscript -e "require(methods);NGScopy::ngscopy_unittest()"
ngscopy_unittest <- function(...){
  pkg <- 'NGScopy'
  test.obj <- UnitTest$new(pkg=pkg,...)
  test.obj$runme()
  invisible()
}

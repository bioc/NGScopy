#!/usr/bin/env Rscript

## ************************************************************************
## This is an executable R script for running `Bioconductor NGScopy'
## at commond line.
## 
## 
## 
## (c) Xiaobei Zhao
## 
## Sat Jun 28 21:19:32 EDT 2014 -0400 (Week 25)
## 
## 
## Reference:
## Bioconductor NGScopy
## http://www.bioconductor.org/packages/release/bioc/html/NGScopy.html
## 
## Get help
## ./ngscopy -h
## $(Rscript -e "cat(Xmisc::get_executable('NGScopy'))") -h
## ************************************************************************


require(methods)
require(Xmisc)

PARSEME <- function(){
  parser <- ArgumentParser$new()

  parser$add_usage('ngscopy [options]')
  parser$add_description(
    'ngscopy scans a pair of input BAM files to detect relative copy number variants.')

  parser$add_argument(
    '--h',type='logical',
    action='store_true',
    help='Print the help page'
  )

  parser$add_argument(
    '--help',type='logical',
    action='store_true',
    help='Print the help page'
  )

  parser$add_argument(
    '--inFpathN',type='character',
    help='The file path to the control (normal) sample'
  )
  
  parser$add_argument(
    '--inFpathT',type='character',
    help='The file path to the case (tumor) sample'
  )
  parser$add_argument(
    '--outFpre',type='character',
    help='The file path of the directory for output'
  )
  parser$add_argument(
    '--libsizeN',type='numeric',
    help='The library size of the control (normal) sample'
  )
  parser$add_argument(
    '--libsizeT',type='numeric',
    help='The library size of the case (tumor) sample'
  )
  parser$add_argument(
    '--mindepth',type='numeric',
    default=20,
    help='The minimal depth of reads per window'
  )
  parser$add_argument(
    '--minsize',type='numeric',
    default=20000,
    help='The minimal size of a window'
  )

  parser$add_argument(
    '--regions',type='character',
    default=NULL,
    help='The regions, a text string or a file path indicating the regions
with three columns (chr/start/end) and without column names.'
  )

  parser$add_argument(
    '--segtype',type='character',
    default="mean.norm",
    help='The type of segmentation. One of c(
"mean.norm","meanvar.norm","mean.cusum","var.css"
). Multiple values are allowed and separated by ",".'
  )
  
  parser$add_argument(
    '--dsN',type='integer',
    default=1,
    help='The downsampling factor of the control (normal) sample'
  )
  
  parser$add_argument(
    '--dsT',type='integer',
    default=1,
    help='The downsampling factor of the test (tumor) sample'
  )
  
  parser$add_argument(
    '--pcThreads',type='integer',
    default=1,
    help='The number of processors performing the parallel computing'
  )

  parser$add_argument(
    '--auto.save',type='logical',
    default=FALSE,
    help='Whether to save (any completed results) automatically'
  )
  
  parser$add_argument(
    '--auto.load',type='logical',
    default=FALSE,
    help='Whether to load (any previously completed results) automatically'
  )
  return(parser)
}


main <- function(){
  parser <- PARSEME()
  parser$helpme()
  
  require(NGScopy)
  
  if (length(regions)){
    regions <- read_regions(regions)
  }
  
  obj <- NGScopy$new(
    outFpre=outFpre,                 
    inFpathN=inFpathN, 
    inFpathT=inFpathT, 
    libsizeN=libsizeN,                
    libsizeT=libsizeT,                
    mindepth=mindepth,                     
    minsize=minsize,
    regions=regions,
    segtype=segtype,
    dsN=dsN,                
    dsT=dsT,     
    pcThreads=pcThreads,
    auto.save=auto.save,
    auto.load=auto.load
  )
  
  obj$write_cn()
  obj$write_seg()
  ## obj$plot_cn()
  obj$plot_cn(ylim=c(-3,3))
  
  invisible()
}


main()

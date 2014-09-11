#!/usr/bin/env Rscript

## ************************************************************************
## This is the R script for running `Case Study I' in `Bioconductor NGScopy'.
## 
## 
## 
## (c) Xiaobei Zhao
## 
## Mon Aug 11 13:44:26 EDT 2014 -0400 (Week 32)
## 
## Reference:
## Bioconductor NGScopy
## http://www.bioconductor.org/packages/release/bioc/html/NGScopy.html
## Case study I: copy number detection of a tumor sample compared to
##               a pooled normal sample
## 
## Runme:
## Rscript ngscopy-case1.R &> ngscopy-case1.out
## 
## ************************************************************************

require(methods)
require(NGScopy)
require(NGScopyData)

## ------------------------------------------------------------------------
## Create an instance of NGScopy class
## ------------------------------------------------------------------------
obj <- NGScopy$new(
  outFpre="ngscopy-case1",         # specified directory for output
  inFpathN=tps_N8.chr6()$bamFpath, # normal sample: tps_90.chr6.sort.bam
  inFpathT=tps_90.chr6()$bamFpath, # tumor sample:  tps_N8.chr6.sort.bam
  libsizeN=5777087,                # the library size of the normal sample
  libsizeT=4624267,                # the library size of the tumor sample
  mindepth=20,                     # the minimal depth of reads per window
  minsize=20000,                   # the minimal size of a window
  regions=read_regions("chr6 41000000 81000000"),
                                   # the regions, set to NULL for the entire genome
  segtype="mean.norm",             # the type of segmentation
  dsN=1,                           # the downsampling factor of the normal 
  dsT=1,                           # the downsampling factor of the tumor
  pcThreads=1                      # the number of processors for computing
  )

obj$show()                         # print the instance

## ------------------------------------------------------------------------
## Review the input
## ------------------------------------------------------------------------

## Get the input files
obj$get_inFpathN()
obj$get_inFpathT()

## Get the library sizes
obj$get_libsizeN()
obj$get_libsizeT()

## Get the windowing parameters
obj$get_mindepth()
obj$get_minsize()

## Get the regions
head(obj$get_regions())

## Get the segmentation type(s)
head(obj$get_segtype())

## Get the down sampling factors
obj$get_dsN()
obj$get_dsT()

## Get the number of processors
obj$get_pcThreads()

## Get the chromosome names of the reference genome 
obj$get_refname()

## Get the chromosome lengths of the reference genome 
obj$get_reflength()


## ------------------------------------------------------------------------
## Process reads in the control (normal) sample (Make windows as well)
## ------------------------------------------------------------------------
obj$proc_normal()                  # this may take a while

## ------------------------------------------------------------------------
## Process reads in the case (tumor) sample
## ------------------------------------------------------------------------
obj$proc_tumor()                   # this may take a while

## ------------------------------------------------------------------------
## Compute/Process the relative copy number ratios and save it
## ------------------------------------------------------------------------

## A data.frame will be saved to file `ngscopy_cn.txt' in the output directory
obj$calc_cn()
obj$proc_cn()
obj$write_cn()

## ------------------------------------------------------------------------
## Compute/Process the segmentation and save it
## ------------------------------------------------------------------------

## A data.frame will be saved to file `ngscopy_seg.txt' in the output directory
obj$calc_seg()
obj$proc_seg()
obj$write_seg() 

## ------------------------------------------------------------------------
## Save the output for later reference
## ------------------------------------------------------------------------
## The NGScopy output is saved as a ngscopy_out.RData file in the output directory
obj$saveme()

## ------------------------------------------------------------------------
## Load and review the output
## ------------------------------------------------------------------------

## Load the output
## (optional if the previous steps have completed in the same R session)
obj$loadme()

## Get the output directory
obj$get_outFpre()

## Get the windows
head(obj$get_windows())

## Get the window sizes
head(obj$get_size())

## Get the window positions (midpoints of the windows)
head(obj$get_pos())

## Get the number of reads per window in the normal
head(obj$get_depthN())

## Get the number of reads per window in the tumor
head(obj$get_depthT())

## Get the copy number
str(obj$get_cn())

## Get the segmentation
str(obj$get_seg())

## Get the data.frame of copy number calling
data.cn <- obj$get_data.cn()
head(data.cn)

## Get the data.frame of segmentation calling
data.seg <- obj$get_data.seg()
head(data.seg)


## ------------------------------------------------------------------------
## Visualize the output
## ------------------------------------------------------------------------
## A figure will be saved to file `ngscopy_cn.pdf' in the output directory
obj$plot_cn(ylim=c(-3,3))       # reset `ylim' to allow full-scale display

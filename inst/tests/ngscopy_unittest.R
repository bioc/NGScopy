#!/usr/bin/env Rscript

## ************************************************************************
## This is the R script for running `Unit Test' in `Bioconductor NGScopy'.
## 
## 
## 
## (c) Xiaobei Zhao
## 
## Wed Aug 06 11:51:20 EDT 2014 -0400 (Week 31)
## 
## Reference:
## Bioconductor NGScopy
## http://www.bioconductor.org/packages/release/bioc/html/NGScopy.html
## 
## Runme:
## Rscript -e "require(methods);NGScopy::ngscopy_unittest()" &> ngscopy_unittest.out
## 
## ************************************************************************

test.ngscopy.1 <- function(){
  ## 
  require("NGScopyData")
  obj <- NGScopy$new(
    outFpre="ngscopy_unittest",      # specified directory for output
    inFpathN=tps_N8.chr6()$bamFpath, # normal sample: tps_90.chr6.sort.bam
    inFpathT=tps_90.chr6()$bamFpath, # tumor sample:  tps_N8.chr6.sort.bam
    libsizeN=5777087,                # the library size of the normal sample
    libsizeT=4624267,                # the library size of the tumor sample
    mindepth=20,                     # the minimal depth of reads per window
    minsize=20000,                   # the minimal size of a window
    regions=read_regions("chr6 41000000 81000000"),
                                     # the regions. NULL for the entire genome
    segtype="mean.norm",             # the type of segmentation
    dsN=1,                           # the downsampling factor of the normal 
    dsT=1,                           # the downsampling factor of the tumor
    pcThreads=1                      # the number of processors for computing
    )

  obj$write_cn()
  obj$write_segm()

  obj$plot_out()                                         # fixed y coordinates
  obj$plot_out("ngscopy_out_-3_3.pdf",ylim=c(-3,3))      # limited y coordinates
  obj$plot_out("ngscopy_out_free_y.pdf",scales='free_y') # free y coordinates
  
}


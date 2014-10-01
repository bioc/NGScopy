#!/usr/bin/env Rscript

## ************************************************************************
## This is the R script for running `Case Study II' in `Bioconductor NGScopy'.
## 
## 
## 
## (c) Xiaobei Zhao
## 
## Mon Aug 11 13:45:24 EDT 2014 -0400 (Week 32)
## 
## Reference:
## Bioconductor NGScopy
## http://www.bioconductor.org/packages/release/bioc/html/NGScopy.html
## Case study II: copy number detection of multiple tumor samples compared to
##                a common pooled normal sample
## 
## Runme:
## Rscript ngscopy-case2.R &> ngscopy-case2.out
## 
## ************************************************************************

require(methods)
require(NGScopy)
require(NGScopyData)

## ------------------------------------------------------------------------
## Process the normal
## ------------------------------------------------------------------------

## Create an instance of `NGScopy' class
obj <- NGScopy$new(pcThreads=1)

## Set the normal sample
obj$set_normal(tps_N8.chr6()$bamFpath)

## Set the regions
regions <- read_regions("chr6 41000000 81000000")
obj$set_regions(regions)

## Set the library size of the normal
obj$set_libsizeN(5777087)

## Specify a directory for output
## It will be saved as "ngscopy_normal.RData" in this directory
obj$set_outFpre("ngscopy-case2/tps_N8")

## Show the input
obj$show()

## Make windows and count reads in the normal
obj$proc_normal()

## Save the output of the normal for later usage
obj$save_normal()

## ------------------------------------------------------------------------
## Process a tumor
## ------------------------------------------------------------------------

## Create an instance of `NGScopy' class
obj1 <- NGScopy$new(pcThreads=1)

## Load the previously saved output of the normal
obj1$load_normal("ngscopy-case2/tps_N8")

## Set a tumor sample (ID: tps_90) and specify a directory for output
obj1$set_tumor(tps_90.chr6()$bamFpath)
obj1$set_outFpre("ngscopy-case2/tps_90")

## Set the library size of the tumor
obj1$set_libsizeT(4624267)

## Show the input
obj1$show()

## Process the tumor
obj1$proc_tumor()

## Process the copy number
obj1$proc_cn()

## Process the segmentation
obj1$proc_segm()

## Plot
obj1$plot_out(ylim=c(-3,3))

## ------------------------------------------------------------------------
## Process a second tumor
## ------------------------------------------------------------------------

## Create another instance of `NGScopy' class
obj2 <- NGScopy$new(pcThreads=1)

## Load the previously saved output of the normal
obj2$load_normal("ngscopy-case2/tps_N8")

## Set a tumor sample (ID: tps_27) and specify a directory for output
obj2$set_tumor(tps_27.chr6()$bamFpath)
obj2$set_outFpre("ngscopy-case2/tps_27")

## Set the library size of the tumor
obj2$set_libsizeT(10220498)

## Show the input
obj2$show()

## Process the tumor
obj2$proc_tumor()

## Process the copy number
obj2$proc_cn()

## Process the segmentation
obj2$proc_segm()

## Plot
obj2$plot_out(ylim=c(-3,3))

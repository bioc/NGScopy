#!/usr/bin/env bash


## ************************************************************************
## This is a commond line example to run bin/ngscopy in `Bioconductor NGScopy'.
## 
## 
## (c) Xiaobei Zhao
## 
## v0.99.1
## Mon Aug 11 09:57:53 EDT 2014 -0400 (Week 32)
## [2014-08-11 09:59:52 EDT] add regions, dsN, dsT
## 
## v0.99.0
## Sat Jun 28 21:19:32 EDT 2014 -0400 (Week 25)
## 
## 
## Reference: 
## Bioconductor NGScopy
## http://www.bioconductor.org/packages/release/bioc/html/NGScopy.html
## 
## Runme:
## bash ngscopy-case1b-cmdline.sh &> ngscopy-case1b-cmdline.out
## 
## ************************************************************************


## ------------------------------------------------------------------------
## The path of the executable, for instance, ${R_LIBS}/NGScopy/bin/ngscopy
## A portable solution is to extract this path by `system.file' in R
## We call such funciton from the command line below.
## ------------------------------------------------------------------------

## ngscopyCmd=$(Rscript -e "cat(system.file('bin','ngscopy',
## package='NGScopy', mustWork=TRUE))")
ngscopyCmd=$(Rscript -e "cat(Xmisc::get_executable('NGScopy'))")
echo ${ngscopyCmd}


## ------------------------------------------------------------------------
## Get help
## ------------------------------------------------------------------------

## ${ngscopyCmd} -h
## ${ngscopyCmd} --help


## ------------------------------------------------------------------------
## An example
## This is a command line version of Case Study I from NGScopy User's Guide.
## See http://www.bioconductor.org/packages/release/bioc/html/NGScopy.html
## ------------------------------------------------------------------------

## a normal sample, given the NGScopyData package is installed
## inFpathN=$(Rscript -e "cat(system.file('extdata','tps_N8.chr6.sort.bam',
## package='NGScopyData', mustWork=TRUE))")
inFpathN=$(Rscript -e "cat(Xmisc::get_extdata('NGScopyData','tps_N8.chr6.sort.bam'))")

## echo ${inFpathN}
## ls -l ${inFpathN}

## a tumor sample, given the NGScopyData package is installed
## inFpathT=$(Rscript -e "cat(system.file('extdata','tps_90.chr6.sort.bam',
## package='NGScopyData', mustWork=TRUE))")
inFpathT=$(Rscript -e "cat(Xmisc::get_extdata('NGScopyData','tps_90.chr6.sort.bam'))")
## echo ${inFpathT}
## ls -l ${inFpathT}

## set pre-computed libsizes based on the original bam files of all chromosomes
libsizeN=5777087
libsizeT=4624267

## set the regions, given the NGScopy package is installed
## regions=$(Rscript -e "cat(system.file('extdata','hg19_chr6_0_171115067.txt',
## package='NGScopy', mustWork=TRUE))")
regions=$(Rscript -e "cat(Xmisc::get_extdata('NGScopy','hg19_chr6_0_171115067.txt'))")
## echo ${regions}

## set downsampling factor; we do not downsample here by setting them to 1.
dsN=1
dsT=1

## set a directory for output
outFpre="ngscopy-case1b-cmdline"

## Run NGScopy given arguments and time it 
time ${ngscopyCmd} --inFpathN=${inFpathN} --inFpathT=${inFpathT} --outFpre="${outFpre}" \
--libsizeN=${libsizeN} --libsizeT=${libsizeT} --regions=${regions} \
--dsN=${dsN} --dsT=${dsT} --pcThreads=1

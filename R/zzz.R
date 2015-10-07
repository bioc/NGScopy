
## ************************************************************************
## 
## 
## 
## (c) Xiaobei Zhao
## 
## Wed Oct 07 12:31:21 EDT 2015 -0400 (Week 40)
## 
## 
## Reference: 
## 
## 
## ************************************************************************

.onAttach <- function(libname, pkgname)
{
  f <- read.dcf(file.path(libname, pkgname, "DESCRIPTION"),
                c("Version", "Date"))
  url <- sprintf('https://www.bioconductor.org/packages/devel/bioc/html/%s.html',pkgname)
  packageStartupMessage(
    'Successfully loaded ', pkgname ,' package version ',f[1,1],'\n',
    'Created on ', f[1,2],'\n',
    'Please refer to the package NEWS for details of changes.','\n',
    'You may find the latest version of ', pkgname ,' at: ','\n',
    url,'\n',
    'To suppress this message: \n',
    'suppressPackageStartupMessages(require(',
    pkgname,'))\n'
    )
}

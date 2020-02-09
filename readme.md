## 0) Prerequisites: Rtools, LaTeX
##
## 1) Install 2 required packages: devtools, roxygen2
## 2) To create PACKAGE (from ICS_smoothing.R file in directory /R)
## 3) Install package "ICSsmoothing"
## 4) Knit vignette cICS_vignette.Rmd
##
## 5) CERN, Check+Test package, Document, GitHub
## 6) Forecasting - Csaba dokoncit!

## 0) Prerequisites: Rtools, LaTeX
##    Package Development Prerequisites
##    https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites
##    http://www.mjdenny.com/R_Package_Pictorial.html
##
## Version of R
version # Ctrl+Enter
#> version.string R version 3.6.1 (2019-07-05)
##
## Install:
## - Rtools from: https://cran.rstudio.com/bin/windows/Rtools/
## - the MikTeX LaTeX distributionfrom: http://miktex.org/download.
## P.S. Rtool is needed only on Windows, on Linux just install GCC
#### install.packages("Rtools", force = TRUE)
#### library("Rtools")
#### Path
#### sysdm.cpl - Advanced - (dole)Environment variables - (dole)Path-Edit

## Notice: devtools installation is broken (On Ubuntu: sudo apt install libxml2-dev libssl-dev libcurl4-openssl-dev libopenblas-dev r-base r-base-dev)
## 1) Install 2 required packages: devtools, roxygen2 (+ digest)
## 1a) devtools for development packages
install.packages("devtools")
library("devtools")
## 1b) roxygen for documentation
##devtools::install_github("klutometis/roxygen",)
devtools::install_github("klutometis/roxygen", force = TRUE)
library(roxygen2)

install.packages("digest")
library(digest)

## 2) To create PACKAGE (from ICS_smoothing.R file in directory /R)
## 2a)To create ICSsmoothing_1.0.tar.gz
## Build - Build Source Package
#> Source package written to C:/Csaba/Hudak_Torok
##
## 2b)To create ICSsmoothing_1.0.zip
## Build - Build Binary Package
#> Source package written to C:/Csaba/Hudak_Torok

## Toto na Serveri nebolo potrebne, len na PC na ustave
## 3) Install package "ICSsmoothing"
##    then load the namespace of the package with name package and attach it on the search list
##install.packages("ICSsmoothing")
#install.packages("C:/Csaba/Hudak_Torok/ICSsmoothing.tar.gz",
#install.packages("C:/Csaba/Hudak_Torok/ICSsmoothing.zip",
#install.packages("C:/Users/CsabaTorok/Desktop/_Cikkek/_____2019/Hudak_Torok/_Pre_Vila/x_Hudak_Kacala_Torok/ICSsmoothing",
getwd()
# On Ubuntu (not suer if needed on Windows)
install.packages("polynom")
# In RStudio: Tools->Install Packages and load the ICSsmoothing.tar.gz
#install.packages("ICSsmoothing",
#                 repos = NULL)
#                 type = "source")
#install.packages("C:/Csaba/Hudak_Torok/ICSsmoothing",
#                 repos = NULL,
#                 type = "source")
###remove.packages("ICSsmoothing", "C:/Csaba/Hudak_Torok/ICSsmoothing")
##remove.packages("ICSsmoothing", lib="~/R/win-library/3.6")
library("ICSsmoothing")
detach("package:ICSsmoothing", unload = TRUE)

## 4) Knit vignette cICS_vignette.Rmd
## Knit - OK

## 5) CERN, Check+Test packege, Document, GitHub
## AirPassengers !?
## CERN !?
##
## Check Package - OK
## Test Package 2 OK, 1 failed - ?
## Document - ?
## GitHub - ?

## 6) Forecasting - Csaba dokoncit!?

# title: "Interpolating cubic splines and data smoothing"
# author: "Juraj Hudák, Viliam Kačala, Csaba Török"




 ICSsmoothing
============

Library for data smoothing using interpolating splines

Building
----------
###  Prerequisites
* [Rtools](https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites)
* LaTeX

### Depends on R packages
* devtools
* roxygen2

### Installation
If using Windows run in R (or RStudio)
```R
install.packages("Rtools", force = TRUE)
library("Rtools")
```
Otherwise if using Linux (Ubuntu) install these packages
```bash
sudo apt install libxml2-dev libssl-dev libcurl4-openssl-dev libopenblas-dev r-base r-base-dev
```

In both cases then run in R 
```R
devtools::install_github("klutometis/roxygen", force = TRUE)
library(roxygen2)
install.packages("digest")
library(digest)
install.packages("polynom")
```
In RStudio install this library using top menu button *Build*->*Install and Restart*.
And use this library as any other
```R
library("ICSsmoothing")
```

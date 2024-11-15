# SEAwise_ecoMSE

*Integrating environmental impacts on stock productivity in Management Strategy Evaluation models* ([ICES webpage](https://www.ices.dk/events/Training/Pages/MSEmodels24.aspx))

Repository for training course material 

## R packages that need to be installed beforehand:

In order to save some time during the course, make sure you have the following R-packages installed. We grouped those by each practical/exercise in case you have some troubles installing a few packages, so you could at least follow the other tutorials. 

### day 1 - 01. practical on env. data pre-processing

- from CRAN:
```
install.packages( c("raster","terra","maps","sf","kohonen","corrplot","PCDimension","devtools"))
# additionally some packages are already archived on CRAN, so you need to install the last available version
devtools::install_version("maptools", version = "1.1.8", repos = "http://cran.us.r-project.org")
devtools::install_version("rgeos", version = "0.6.4", repos = "http://cran.us.r-project.org")
```
- from github: 
```
devtools::install_github("BernhardKuehn/marmalaid")
devtools::install_github("coatless-rd-rcpp/rcpp-and-doparallel")
devtools::install_github("BernhardKuehn/fast.EOT")
```
### day 1 - 02. practical on bias-correction

### day 1 - 03. practical on integrating env. uncertainty via BVARs

### day 2 - 04. practical on ...

### day 2 - 05. practical on ...

...




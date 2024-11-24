# SEAwise_ecoMSE
**Warning: this information will not be complete until Friday Nov 22nd, so students will need to check back at that time.**

*Integrating environmental impacts on stock productivity in Management Strategy Evaluation models* ([ICES webpage](https://www.ices.dk/events/Training/Pages/MSEmodels24.aspx))

Repository for training course material

## R packages that need to be installed beforehand:

In order to save some time during the course, make sure you have the following R-packages installed. We grouped those by each practical/exercise in case you have some troubles installing a few packages, so you could at least follow the other tutorials. 

## Day 1
-----

### 01. practical on env. data pre-processing

- from CRAN:
```R
install.packages( c("raster","terra","maps","sf","kohonen","corrplot","PCDimension","devtools"))
# additionally some packages are already archived on CRAN, so you need to install the last available version
devtools::install_version("maptools", version = "1.1.8", repos = "http://cran.us.r-project.org")
devtools::install_version("rgeos", version = "0.6.4", repos = "http://cran.us.r-project.org")
```
- from github: 
```R
devtools::install_github("BernhardKuehn/marmalaid")
devtools::install_github("coatless-rd-rcpp/rcpp-and-doparallel")
devtools::install_github("BernhardKuehn/fast.EOT")
```
### 02. practical on bias-correction

- from CRAN:
```R
install.packages(c("raster","terra","maps",
                   "reshape2","qmap","ggplot2",
                   "devtools","parallel","doParallel","foreach"))
# additionally some packages are already archived on CRAN, so you need to install the last available version
devtools::install_version("maptools", version = "1.1.8", repos = "http://cran.us.r-project.org")
devtools::install_version("rgeos", version = "0.6.4", repos = "http://cran.us.r-project.org")
```
- from github:
```R
devtools::install_github("BernhardKuehn/marmalaid")
```
### 03. practical on integrating env. uncertainty via BVARs

- from CRAN:
```R
install.packages(c("raster","terra","maps","plyr",
                   "BVAR","scam","mgcv","reshape2",
                   "devtools","parallel","doParallel","doSnow",
                   "foreach","funtimes","coda","tsDyn",
                   "gsignal","FNN","Rcpp","RcppEigen",
                   "RcppArmadillo","ggplot2","colorspace","patchwork"))
```
- from github:
```
devtools::install_github("BernhardKuehn/marmalaid")
```

## Day 2
-----

### 01. practical on model fitting to growth data
```R
install.packages(c(
  "glmmTMB", "RTMB", 
  "MuMIn", "bbmle",
  "broom.mixed",
  "doParallel", "foreach",
  "ggplot2", "patchwork"
  ))

```

### 02. practical on model fitting to recruitment
```R
install.packages(c("modelr","nlstools","MuMIn","rlist","formula.tools","corrplot"))
```

### 03. practical on model selection

```R
install.packages(c("corrplot", "doRNG", "dplyr", "forecast", "GA", "ggplot2", 
  "glmmTMB", "icesAdvice", "knitr", "memoise", "MuMIn", "parallel", 
  "patchwork", "tidyr"))
```

## Day 3
-----

### 01. Practical on including environmental forcing into stock recruitment in MSE

To prepare for the section where we will run some MSEs with environmentally mediated data, you will need to install and experimental version of smsR in R, and read the WKECOMSE report available on the sharepoint.
To install smsR directly from Github you need to have the library remotes installed.

```R
if('remotes' %in% installed.packages()[,"Package"] ) == FALSE) install.packages('remotes')
remotes::install_github("https://github.com/nissandjac/smsR", dependencies = TRUE) # Change this link to the development version
```
This will install smsR and the required packages needed to run the model (e.g., TMB).
This installs the package smsR, which is a seasonal stock assessment model, that can be used to regular stock assessments, but it also has the potential to 1) include seasonality, 2) include environmental covariates in the recruitment estimation, 3) various options for natural mortality estimation as random walk, or with predator input.

To test if the model is properly working you can run the sandeel stock assessment included in the package as

```R
library(smsR)

# Set up sandeel for area 1r
df.tmb <- get_TMB_parameters(
  mtrx = sandeel_1r$lhs, # List that contains M, mat, west, weca
  Surveyobs = sandeel_1r$survey, # Survey observations
  Catchobs = sandeel_1r$Catch, # Catch observations
  years = 1983:2021, # Years to run
  nseason = 2, # Number of seasons
  useEffort = TRUE, # Use effort to calculate F
  ages = 0:4, # Ages of the species
  recseason = 2, # Season where recruitment occurs
  CminageSeason = c(1, 1), # Minimum catch age per season
  Fmaxage = 3, # Fully selected fishing mortality age
  Qminage = c(0, 1), # Minimum age in surveys
  Qmaxage = c(1, 3), # Max age in surveys
  Fbarage = c(1, 2), # Age use to calculate Fbar
  effort = sandeel_1r$effort, # Effort input
  blocks = c(1983, 1999), # Blocks with unique selectivity
  nocatch = sandeel_1r$nocatch, # Seasons where F is not calculated
  surveyStart = c(0.75, 0), #
  surveyEnd = c(1, 0), #
  surveySeason = c(2, 1), #
  surveyCV = list(c(0, 1), c(1, 2)),
  catchCV = list(c(1, 3), c(1, 3)), # Catch CV groupings
  estCV = c(0, 2, 0), # Estimate CVs for 1) survey, 2) catch, 3) Stock recruitment relationship
  beta = 105809, # Hockey stick break point
  nllfactor = c(1, 1, 0.05) # Factor for relative strength of log-likelihood
)

parms <- getParms(df.tmb)
sas <- runAssessment(df.tmb, parms)

```

And if the run was successful you can plot the output results easily as

```R
plot(sas) # Plots general stock assessment output

```

This is a fairly simple model that runs without any climatic or other external data, but works as a baseline assessment model. 

### 03. practical on MSEs with FLBEIA

- from CRAN:
```R
install.packages( c("plyr","stringr","reshape2","glmmTMB",
                    "earth","ggplot2","patchwork"))
```
- from FLR:
```R
install.packages( c("FLCore", "FLFleet", "FLBEIA",
                    "FLash", "FLAssess"), 
                  repos="http://flr-project.org/R")
```

# SEAwise_ecoMSE

*Integrating environmental impacts on stock productivity in Management Strategy Evaluation models* ([ICES webpage](https://www.ices.dk/events/Training/Pages/MSEmodels24.aspx))

Repository for training course material

## day 3 - 01. Practical on including environmental forcing into stock recruitment in MSE

To prepare for the section where we will run some MSEs with environmentally mediated data, you will need to install and experimental version of smsR in R, and read the WKECOMSE report available on the sharepoint.
To install smsR directly from Github you need to have the library remotes installed.

```
if('remotes' %in% installed.packages()[,"Package"] ) == FALSE) install.packages('remotes')
remotes::install_github("https://github.com/nissandjac/smsR", dependencies = TRUE) # Change this link to the development version
```
This will install smsR and the required packages needed to run the model (e.g., TMB).
This installs the package smsR, which is a seasonal stock assessment model, that can be used to regular stock assessments, but it also has the potential to 1) include seasonality, 2) include environmental covariates in the recruitment estimation, 3) various options for natural mortality estimation as random walk, or with predator input.

To test if the model is properly working you can run the sandeel stock assessment included in the package as

```
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

```
plot(sas) # Plots general stock assessment output

```

This is a fairly simple model that runs without any climatic or other external data, but works as a baseline assessment model. 

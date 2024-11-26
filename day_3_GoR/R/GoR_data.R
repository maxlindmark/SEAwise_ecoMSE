#### Load the GoR herring data


GoR_data <- function(wd, maxage, years){

  require(tidyverse)
  ages <- 0:maxage
  nage <- length(ages)
  nyear <- length(years)


  canum <- read.table(file.path(wd,'canum.dat'), skip = 5)
  # Assign ages and years
  names(canum) <- 1:maxage
  canum$years <- years # One year is missing?


  # Ignore catch tonnes for now. Is it double used?

  mat.tmp <- read.table(file.path(wd, 'maturity.dat'), skip = 5)
  names(mat.tmp) <- 1:maxage
  mat.tmp$years <- years
  mat.tmp <- addYear(mat.tmp, nyear)
  # Add one forecast year\

  # Natural mortaliyt

  M2 <- read.table(file.path(wd, 'natural_mortality.dat'), skip = 5)
  names(M2) <- 1:maxage
  M2$years <- years
  M2 <- addYear(M2, nyear)

  # Proportion of F before spawning
  propFin  <- read.table(file.path(wd, 'propF.dat'), skip = 5)
  names(propFin) <- 1:maxage
  propFin$years <- c(years)
  propFin <- addYear(propFin, nyear)



  # Proportion of F before spawning
  propMin  <- read.table(file.path(wd, 'propM.dat'), skip = 5)
  names(propMin) <- 1:maxage
  propMin$years <- c(years)
  propMin <- addYear(propMin, nyear)

  # Survey
  # Just do the newly organized ones for now

  survey.trap <- read.csv(file.path(wd,'trap_net_fixed.csv'), sep = ';')
  survey.trap$survey <- 'trap'

  # Multiply the effort
  survey.trap[,2:8] <- survey.trap[,2:8]*survey.trap[,1]

  survey.acou <- read.csv(file.path(wd, 'acoustic_survey.csv'), sep = ';')
  survey.acou$survey <- 'acoustic'

  # Combine the surveys
  # survey <- bind_rows(survey.trap, survey.acou)
  # survey <- survey[,c(11,2:8,1,9:10)]
  # Remove 2021 for now
  #survey <- survey[-which(survey$years == 2021),]
  survey <- survey.acou[,-1] # Remove the ones


  # weight in catch
  weca.tmp <- read.table(file.path(wd, 'weca.dat'), skip = 5)
  names(weca.tmp) <- 1:maxage
  weca.tmp$years <- years
  weca.tmp <- addYear(weca.tmp, nyear, method = 'average', avg_years = 3)

  # Weight in survey
  west.tmp <- read.table(file.path(wd, 'west.dat'), skip = 5)
  names(west.tmp) <- 1:maxage
  west.tmp$years <- years
  west.tmp <- addYear(west.tmp, nyear, method = 'average', avg_years = 3)


  # Catch observations  (dimensions age, year, quarter)
  Catchobs <- array(NA, dim = c(length(ages), length(years), 1))
  Catchobs[2:length(ages),,] <- as.matrix(t(canum[,1:8]))

  # Survey observations

  # Survey observations (dimensions age, year, quarter, number of surveys)
  surveys <- unique(survey$survey)
  Surveyobs <- array(NA, dim = c(length(ages), length(years), length(surveys)))
  #Effort <- array(NA, dim = c(length(years), length(surveys)))

  for(i in 1:length(surveys)){


    sstmp <- survey[survey$survey == surveys[i],]
    Surveyobs[2:length(ages), which(years %in% sstmp$years), i] <- as.matrix(t(sstmp[,1:(length(ages)-1)]))



  }

  Surveyobs[is.na(Surveyobs)] <- -1
  dimnames(Surveyobs) <- list(ages,
                              years,
                              'acoustic')

  # Reformulate to matrices for TMB
  weca <- west<- mat <- M <- propF <- propM  <- array(0, dim = c(length(ages), length(years)+1, 1))
  weca[2:length(ages),,] <- as.matrix(t(weca.tmp[,1:maxage]))
  west[2:length(ages),,] <- as.matrix(t(west.tmp[,1:maxage]))
  M[2:length(ages),,] <- as.matrix(t(M2[,1:maxage]))
  M[1,,] <- M[2,,]
  mat[2:length(ages),,] <- as.matrix(t(mat.tmp[,1:maxage]))
  propF[2:length(ages),,] <- as.matrix(t(propFin[,1:maxage]))
  propM[2:length(ages),,] <- as.matrix(t(propMin[,1:maxage]))

  mtrx <- list(weca = weca,
               west = west,
               M = M,
               mat = mat,
               Catchobs = Catchobs,
               Surveyobs = Surveyobs)

  return(list(mtrx = mtrx,
              propF = propF,
              propM = propM)
  )
}

mtrx_forward <- function(mtrx, nyear,nage,nseason,
                         vnew){

  mtrx.new <- mtrx


  mtrx.new$weca <- addYear(mtrx$weca, nyear, method = 'average')
  mtrx.new$west <- addYear(mtrx$west, nyear, method = 'average')
  mtrx.new$M <-  addYear(mtrx$M, nyear, method = 'last')
  mtrx.new$mat <- addYear(mtrx$mat, nyear, method = 'last')
  mtrx.new$Fsel <- addYear(mtrx$Fsel, nyear, method = 'last')
  mtrx.new$Fin <- addYear(mtrx$Fin, nyear, method = 'last')
  mtrx.new$propF <- addYear(mtrx$propF, nyear, method = 'last')
  mtrx.new$propM <- addYear(mtrx$propM, nyear, method = 'last')




  return(mtrx.new)
}


estimateRecruitment <- function(df.mse, SSB, b = 1, method = 'hockey', env_new = NULL, beta = NULL,
                                alpha = df.mse$R0){

  err <- rnorm(1, mean = 0, sd = df.mse$SDR)

  if(length(b) != length(SSB)){
    warning('bias adjustment has wrong length. Assuming its 1 for all years')
    b <- rep(b, length(SSB))
  }


  if(method == 'hockey'){

    if(is.null(alpha)){
      stop('provide numerical value plateau parameter for hockey stick')
    }

    R  <- log(alpha) + log(SSB)

    if(SSB > df.mse$betaSR){
      R <- log(alpha)+log(df.mse$betaSR)
    }


    R <- exp(R)* exp(-0.5*b*df.mse$SDR^2+err)

  }

  if(method == 'median'){

    R <- median(df.mse$Rin)* exp(-0.5*b*df.mse$SDR^2+err)

  }

  if(method =='hockey_env'){

    R <- alpha * exp(err+sum(beta*env_new)-0.5*b*df.mse$SDR^2)
    # if(SSB < 60000){ # Fix this later
    #
    # }


  }

  if(method == 'BH_env'){

    env_tot <- sum(df.mse$beta_env*env_new)

    R <- (4*df.mse$h*df.mse$R0*SSB/(as.numeric(df.mse$SSB0)*(1-df.mse$h)+ SSB*(5*df.mse$h-1)))*exp(-0.5*df.mse$SDR+err+env_tot)

  }



  return(R)
}


OM.to.tmb <- function(OM,
                      df.mse,
                      df.assess){

  #
  nyear <- length(df.mse$years)
  Surveyobs <- array(-1 , dim = c(df.assess$nage, nyear, df.assess$nsurvey))

  for(k in 1:df.assess$nsurvey){
    Surveyobs[,1:(nyear-1),k] <- df.assess$Surveyobs[,,k]
    Surveyobs[,nyear,k] <- OM$survey[,nyear,k]
  }

  Surveyobs[is.na(Surveyobs)] <- -1


  Catchobs <- array(0 , dim = c(df.assess$nage, nyear, df.assess$nseason))
  Catchobs[,1:(nyear-1),] <- df.assess$Catchobs
  Catchobs[,nyear,] <- OM$CatchN.save.age[,nyear,1,]
  Catchobs[is.na(Catchobs)] <- -1
  # Group index and effort #
  weca <- df.mse$weca
  west <- df.mse$west

  # propM and propF
  df.assess$propF <- abind::abind(df.assess$propF, array(df.assess$propF[,df.assess$nyears-1,], dim = c(df.assess$nage, 1, df.assess$nsurvey)), along = 2)
  df.assess$propM <- abind::abind(df.assess$propM, array(df.assess$propM[,df.assess$nyears-1,], dim = c(df.assess$nage, 1, df.assess$nsurvey)), along = 2)


  # Insert new values into old data fra


  df.assess$weca <- weca
  df.assess$west <- west
  df.assess$Surveyobs <- Surveyobs
  df.assess$Catchobs <- Catchobs
  df.assess$years <- c(df.assess$years,max(df.assess$years)+1)
  df.assess$nyears <- length(df.assess$years)
  df.assess$scv <- abind::abind(df.assess$scv, array(0, dim = c(df.assess$nage, 1, df.assess$nsurvey)), along = 2)
  df.assess$M_matrix <- abind::abind(df.assess$M_matrix,
                                     array(0, dim = c(df.assess$M_nparms, 1)), along = 2)

  # Catches
  df.assess$effort <- rbind(df.assess$effort, c(1)) # Fix these fishing mortalities later
  df.assess$nocatch <- rbind(df.assess$nocatch, c(1))
  df.assess$bidx <- c(df.assess$bidx, max(df.assess$bidx))
  df.assess$M <- df.mse$M
  df.assess$Mat <- df.mse$mat
  df.assess$env_matrix <- df.mse$env


  return(df.assess)

}


#' Helper function for adding new years to a stock assessment object. Useful for simulation.
#'
#' @param dat
#' @param nyears
#' @param method
#' @param avg_years
#'
#' @return
#'
#'
#' @examples
#'
addYear <- function(dat, nyears, method = 'last', avg_years = 3){

  if(is.data.frame(dat)){

    if(method == 'last'){
      tmp <- dat[nyears, ]
      tmp$years <- tmp$years+1
      dat.out <- rbind(dat, tmp)
    }

    if(method == 'average'){

      tmp <- (colMeans(dat[(nyears-avg_years):nyears,]))
      tmp <- as.data.frame(t(tmp))
      tmp$years <- max(dat$years)+1
      dat.out <- rbind(dat, tmp)

    }

  }

  if(is.array(dat)){

    if(method == 'last'){
      nage <- dim(dat)[1]

      tmp <- array(dat[,nyears,1], dim = c(nage, 1, 1))
      dat.out <- abind::abind(dat, tmp, along = 2)

    }

    if(method == 'average'){

      if(dim(dat)[3] > 1){
        stop('only 1 season supported currently')
      }

      nage <- dim(dat)[1]
      tmp <- array(((rowMeans(dat[,(nyears-avg_years):nyears,1]))), dim = c(nage, 1, 1))
      dat.out <- abind::abind(dat, tmp, along = 2)

    }



  }



  return(dat.out)
}


#' Helper function to turn list into MSE dataframe
#'
#' @param mse
#' @param nruns
#'
#' @return
#' @export
#'
#' @examples
unlist_MSE <- function(mse, nruns){



  for(i in 1:nruns){

    if(i == 1){
      df.out <- mse[[i]][[1]]
      df.N <- mse[[i]][[2]]
      df.parms <- mse[[i]][[3]] %>% mutate(year_asses = df.out$year_assess[1])
    }else{
      df.out <- rbind(df.out, mse[[i]][[1]])
      df.N <- rbind(df.N, mse[[i]][[2]])
      df.parms <- rbind(df.parms, mse[[i]][[3]] %>% dplyr::mutate(year_asses = df.out$year_assess[nrow(df.out)]))
    }

  }


  return(list(df.out,
              df.N,
              df.parms))
}



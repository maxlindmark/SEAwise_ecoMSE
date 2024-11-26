run.sms.MSE <- function(df.tmb,
                        sas,
                        mse_years = 30,
                        nruns = 100,
                        HCR = 'Fmsy',
                        recruitment = 'hockey',
                        shortcut = FALSE,
                        n.cores = 4,
                        seeds = NULL,
                        wd){

  # Initialize the model

  ptm <- proc.time()

  my.cluster <- parallel::makeCluster(
    n.cores,
    type = "PSOCK"
  )

  plist <- c('TMB',
             'tidyverse',
             'smsR','abind','dplyr')

  doParallel::registerDoParallel(cl = my.cluster)


  if(is.null(seeds)){
    seeds <-  round(runif(nruns, min = 1, max = 1e6))
  }

  # Survey catchability
  Q <- getCatchability(df.tmb, sas)
  Q <- array(Q$Q, dim = c(df.tmb$nage, df.tmb$nsurvey), )
  Q[is.na(Q)] <- 0

  ages <- df.tmb$age
  nage <- df.tmb$nage

  # Overwrite later
  h <- NULL
  # Estimated parameters to use in OM
  parms.true <- getEstimatedParms(sas)

  if(df.tmb$recmodel == 1){
    alpha <- exp(parms.true$value[parms.true$parameter == 'logalpha'])
    Ninit <- exp(parms.true$value[parms.true$parameter == 'logNinit'])
    SDR <- exp(parms.true$value[parms.true$parameter == 'logSDrec'])
  }



  if(df.tmb$recmodel == 2){
    alpha <- exp(parms.true$value[parms.true$parameter == 'logalpha'])
    beta <- parms.true$value[parms.true$parameter == 'env']
    Ninit <- exp(parms.true$value[parms.true$parameter == 'logNinit'])
    SDR <- exp(parms.true$value[parms.true$parameter == 'logSDrec'])
  }

  if(df.tmb$recmodel == 3){
    alpha <- exp(parms.true$value[parms.true$parameter == 'logR0'])
    SDR <- exp(parms.true$value[parms.true$parameter == 'logSDrec'])
    Ninit <- alpha * exp(-(df.tmb$M[2:df.tmb$nage,1,1]*ages[2:df.tmb$nage])-0.5*SDR) *
      exp(sas$reps$par.random[names(sas$reps$par.random) == 'logNinit'])
    # Plus groups
    Ninit[df.tmb$nage-1] <- alpha*exp(-(df.tmb$M[nage,1,1]*(df.tmb$nage-1))-0.5*SDR)/(1-exp(-df.tmb$M[nage,1,1]))*
      exp(sas$reps$par.random[names(sas$reps$par.random) == 'logNinit'])[df.tmb$nage-1]



    # assumme steepness is mapped for now
    h <- exp(sas$obj$env$parList()$logh)



  }


  if(is.null(h)){
    h <- exp(parms$logh)
  }

  F0 <- getF(df.tmb, sas)

  Fsel <- getSel(df.tmb,sas)
  Fsel <- array(Fsel$Fsel, dim = c(df.tmb$nage,df.tmb$nyears,1))

  # Turn F0 into array
  F0 <- array(F0$F0, dim = c(df.tmb$nage, df.tmb$nyears, 1))


  mtrx <- list(weca = df.tmb$weca,
               west = df.tmb$west,
               M = df.tmb$M,
               mat = df.tmb$Mat,
               Fsel = Fsel,
               Fin = F0,
               propF = df.tmb$propF,
               propM = df.tmb$propM)

  mtrx$weca[is.na(mtrx$weca)] <- 0
  mtrx$west[is.na(mtrx$west)] <- 0

  move <- FALSE

  R_est <- getR(df.tmb, sas)

  df.OM <- list(
    years = df.tmb$years,
    last_year = max(df.tmb$years),
    nseason = df.tmb$nseason,
    age = df.tmb$age,
    nage = length(df.tmb$age),
    F0 = mtrx$Fin,
    sel = mtrx$Fsel,
    M = mtrx$M,
    mat = mtrx$mat,
    weca = mtrx$weca,
    west = mtrx$west,
    propF = df.tmb$propF,
    propM = df.tmb$propM,
    Fbarage = df.tmb$Fbarage,
    betaSR = df.tmb$betaSR,
    nsurvey = df.tmb$nsurvey,
    surveyStart = df.tmb$surveyStart,
    surveyEnd = df.tmb$surveyEnd,
    surveySD = 0.4,
    surveySeason = df.tmb$surveySeason,
    Q = Q,
    recruitment = 'estimated',
    rseason = df.tmb$recseason,
    rec.space = 1,
    h = h,
    Fmodel = 'est',
    Ninit = c(0,
              Ninit),
    Rin = R_est$R,
    move = move,
    R0 =  alpha,
    SDR = SDR,
    b = rep(0, df.tmb$nyears)
  )

  tmp <-run.agebased.sms.op(df.OM)

  Fmsy <- df.tmb$Fmsy


  if(df.tmb$recmodel == 3){
    df.OM$beta_env <- parms.true$value[parms.true$parameter == 'env']
    df.OM$SSB0 <- sas$reps$value[names(sas$reps$value) == 'SSB0']

  }


  ls <- foreach(
    runs = 1:nruns,
    # .combine = 'rbind',
    .packages = plist) %dopar% {

      set.seed(seeds[runs])


      source(file.path(wd,'R/mse_tools.R'))

      for(yr in 1:mse_years){

        if(yr == 1){
          mtrx.new <- mtrx
          df.mse <- df.OM
          df.assess <- df.tmb
          parms.mse <- getParms(df.tmb)
          tmp <- tmp
          Fnew <- df.OM$F0
          rec <- df.OM$Rin
          env <- df.tmb$env_matrix
          Fmsy <- df.tmb$Fmsy
        }else{
          df.assess <- OM.to.tmb(tmp, df.mse, df.assess)

        }


        #if(shortCut == FALSE){
        parms.mse <- getParms(df.assess) # Get new parameter sizes

        assess<- runAssessment(df.assess, parms.mse)

        # Throw a warning if it doesn't converge
        if(assess$reps$pdHess == FALSE | max(assess$reps$gradient.fixed) > 1){
          warning('not converged')
          badmodel <- 1
          # stop()
        }else{
          badmodel <- 0
        }

        # #
        # }
        nyear <- length(df.mse$years)
        # run the operating model one more year
        years.new <- c(df.mse$years, df.mse$years[nyear]+1)

        # Predict recruitment in the new year in the operating model

        # Add new environmental variables

        # Figure out what recruitment model is used

        mps <- getMPS(df.assess, parms.mse )
        env_new <- apply(df.assess$env_matrix,MAR = 1,FUN = sample, size = 1)

        if(df.assess$recmodel == 3){
        Rnew <- estimateRecruitment(df.mse, as.numeric(tmp$SSB[nyear]), method = 'BH_env', env_new = env_new)
        }

        if(df.assess$recmodel == 1){
          Rnew <- estimateRecruitment(df.mse, as.numeric(tmp$SSB[nyear]), method = 'hockey')
        }


        rec <- c(rec, Rnew)
        # Add another year to all matrices
        mtrx.new <- mtrx_forward(mtrx = mtrx.new,
                                 nyear = nyear)

        # Insert the HCR here. Just use Fmsy for easies
        # Calculate Fmsy

        if(HCR == 'Fmsy'){

          if(is.null(Fmsy)){
            stop('Provide Fmsy for this option')
          }


          fcast <- getTAC(df.assess,
                          assess,
                          recruitment = 'mean',
                          avg_R = 5,
                          HCR = 'Fmsy')
          Fsel <- fcast$Fnext
          TAC <- fcast$TAC
        }


        if(HCR == 'Fmean'){
          fcast <- getTAC(df.assess,
                          assess,
                          recruitment = 'mean',
                          avg_R = 5,
                          HCR = 'Fmean')

          TAC <- fcast$TAC
          Fsel <- fcast$Fnext


        }

        if(HCR == 'Fmean_last'){
          if(yr == 1){
            fcast <- getTAC(df.assess,
                            assess,
                            recruitment = 'mean',
                            avg_R = 5,
                            HCR = 'Fmean')

            TAC <- fcast$TAC
            Fsel <- fcast$Fnext
          }

        }

        if(HCR == 'Fmsy_ext'){
            fcast <- getTAC(df.assess,
                            assess,
                            recruitment = 'mean',
                            avg_R = 5,
                            HCR = 'Fmean')

            TAC <- fcast$TAC
            Fsel <- fcast$Fnext
            fmul <- df.tmb$Fmsy/mean(Fsel[(df.mse$Fbarage[1]+1):(df.mse$Fbarage[2]+1)])
            Fsel <- Fsel*fmul

        }

        Fnew <- abind::abind(Fnew, Fsel, along = 2)
        # Run the rest of the year in the OM
        env[,nyear-1] <- env_new
        env <- cbind(env, rep(1, length(env_new)))

        df.mse <-  list(
          years = years.new,
          nseason = df.tmb$nseason,
          age = ages,
          nage = length(ages),
          F0 = Fnew,
          sel = mtrx.new$Fsel,
          M = mtrx.new$M,
          mat = mtrx.new$mat,
          weca = mtrx.new$weca,
          west = mtrx.new$west,
          sel = mtrx.new$Fsel,
          propF = mtrx.new$propF,
          propM = mtrx.new$propM,
          betaSR = df.OM$betaSR,
          nsurvey = df.tmb$nsurvey,
          Fbarage = df.tmb$Fbarage,
          surveySeason = df.tmb$surveySeason,
          surveyStart = df.tmb$surveyStart,
          surveyEnd = df.tmb$surveyEnd,
          env = env,
          surveySD = 0.4,
          Q = Q,
          recruitment = 'estimated', # Recruitment is provided above
          rseason = df.tmb$recseason,
          Fmodel = 'est',
          beta_env = df.OM$beta_env,
          SSB0 = df.OM$SSB0,
          Ninit = df.OM$Ninit,
          h = df.OM$h,
          last_year = max(years.new),
          rec.space = 1,
          Rin = rec,
          move = FALSE,
          R0 =  df.OM$R0,
          SDR = df.OM$SDR,
          b = rep(0, length(years.new))

        )
        # Find OM F

        if(HCR == 'Fmean_last'){
          df.mse$F0[,length(years.new),] <- Fsel
        }

        if(HCR == 'Fmean' | HCR == 'Fmsy'){
          Fgo <- getOM_FTAC(TAC, df.mse)
          df.mse$F0[,length(years.new),] <- Fgo
        }

        if(HCR =='Fmsy_ext'){

          df.mse$F0[,length(years.new),] <- Fsel
          #TAC <- NA
        }
        # Run the OM for the new

        # Run the OM for the new year
        tmp <- run.agebased.sms.op(df.mse)
        #print(tmp$Fbar[length(df.mse$years)])



        SSB.tmp <- getSSB(df.assess,assess)
        R.tmp <- getR(df.assess,assess)
        Fbar.tmp <- getFbar(df.assess, assess)
        Catch.tmp <- getCatch(df.assess, assess)
        # Dont include the predicted year from the model
        N <- getN(df.assess, assess)
        CatchN <- getCatchN(df.assess, assess)


        if(yr == 1){
          df.tac <- data.frame(TAC = TAC,
                               Ctrue = tmp$Catch[length(df.mse$years)],
                               Cest = Catch.tmp$Catch[df.assess$nyears],
                               yr = yr)
        }else{
          df.tac <- rbind(data.frame(TAC = TAC,
                                     Ctrue = tmp$Catch[length(df.mse$years)],
                                     Cest = Catch.tmp$Catch[df.assess$nyears],
                                     yr = yr), df.tac)
        }

        # Save results for R and SSB
        if(yr == mse_years){

          # Calc Fmsy for the last year
          # fcast <- getTAC(df.assess,
          #                 assess,
          #                 recruitment = 'mean',
          #                 avg_R = 5,
          #                 HCR = 'Fmsy_ext'
          # )


          assess.out <- data.frame(years = SSB.tmp$years[1:df.assess$nyears],
                                   SSB = SSB.tmp$SSB[1:df.assess$nyears],
                                   R = R.tmp$R[1:df.assess$nyears],
                                   Fbar = Fbar.tmp$Fbar,
                                   Fmsy = Fmsy,
                                   Catch = Catch.tmp$Catch,
                                   model = 'EM',
                                   year_assess = paste('run',runs,sep = '-'))

          N.out <- data.frame(years = N$years,
                              N = N$N,
                              Ntrue = as.numeric(tmp$N.save.age[,1:(length(df.mse$years)-1),1,1]),
                              weight = as.numeric(df.mse$west[,1:length(df.mse$years)-1,])
          )
          #
          #



          df.out <- rbind(
            data.frame(
              years = df.assess$years,
              SSB = tmp$SSB[1:df.assess$nyears,1],
              R = tmp$R.save[1:df.assess$nyears],
              Fbar = tmp$Fbar[1:df.assess$nyears],
              Fmsy= NA, # Add the OM Fmsy later
              Catch = tmp$Catch[1:df.assess$nyears],
              model = 'OM',
              year_assess = paste('run', runs, sep = '-')
            ),
            assess.out)

          # Estimated parameterr

          parms.mse.out <- getEstimatedParms(assess)

          #
          #           assessN.out <- data.frame(years = N$years,
          #                                    N = N$N,
          #                                    model = 'EM',
          #                                    year_assess = paste('run',runs,sep = '-'))
          #
          #           assessN.out.OM <- rbind(data.frame(years = N$years,
          #                                     N = as.numeric(tmp$N.save[,1:df.assess$nyears]),
          #                                     model = 'EM',
          #                                     year_assess = paste('run',runs,sep = '-')),
          #                                   assessN.out)






        }#else{
        #   df.out <- rbind(df.out,
        #                   data.frame(
        #                     years = df.assess$years,
        #                     SSB = tmp$SSB[1:df.assess$nyears,1],
        #                     R = tmp$R.save[1:df.assess$nyears],
        #                     Fbar = tmp$Fbar[1:df.assess$nyears],
        #                     model = 'OM',
        #                     year_assess = paste('year', yr, sep = '-')
        #                   ),
        #                   assess.out
        #   )
        # }



      }


      return(list(df = df.out,
                  N.out = N.out,
                  parms = parms.mse.out))

    }

  parallel::stopCluster(cl = my.cluster)

  tt <- proc.time() - ptm
  print(paste('#### run time = ',round(tt[3],1), 's ####', sep = ''))


  return(ls)
}





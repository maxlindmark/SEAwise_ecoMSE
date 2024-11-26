library(RTMB)
TapeConfig(atomic = "disable")# to make %*% maintain sparseness of matricies (expands operations to notice zeros)
#run with and without this TapeConfig to check consistency

# function definitions --------------------------
nll = function(params, data, formLinf, mod) {
  getAll(data, params, warn=TRUE)
  l2 = OBS(l2)
  
  #initialize joint negative log likelihood inside the nll funciton
  nllin = 0
  #combine predictors for Linf
  X = model.matrix(formLinf, data=data)
  Zcohort = model.matrix(~as.factor(cohort)+0, data=data)
  Zyear = model.matrix(~as.factor(year)+0, data=data)
  
  XbetaZb = X %*% beta + Zcohort %*% bcohort + Zyear %*% byear #fixed and random effects
  Xbeta = X %*% beta #only fixed effects
  
  #transform parameters
  K = exp(logK)
  D = exp(logD)
  
  #models
  predl2 = switch(mod,
                  #special VBGF: l2 ~ l+((Linf + C1*covar1)-l) * (1-exp(-K*dt))
                  svbgf = l+ (exp(XbetaZb)-l) * (1-exp(-K*dt)),
                  
                  #generalized VBGF: l2 ~ ((Linf+covar1*C1)^(1/D)*(1-exp(-K*dt))+l^(1/D)*exp(-K*dt))^D
                  gvbgf = ((exp(XbetaZb))^(1/D)*(1-exp(-K*dt))+l^(1/D)*exp(-K*dt))^D,
                  
                  #Gompertz: l2 ~ exp(log(Linf+C1*covar1)*(1-exp(-K*dt))+log(l)*exp(-K*dt))
                  gomp = exp(log(exp(XbetaZb))*(1-exp(-K*dt))+log(l)*exp(-K*dt))
  )
  nllin = nllin - sum(dnorm(bcohort, mean=0, sd = exp(logsdcohort), log=TRUE))
  nllin = nllin - sum(dnorm(byear, mean=0, sd = exp(logsdyear), log=TRUE))
  nllin = nllin - sum(dnorm(l2, mean = predl2, sd = exp(logsdresid), log=TRUE), na.rm=TRUE)
  
  #calculate population level predictions for newdata
  totvar = exp(logsdcohort)^2 + exp(logsdyear)^2 
  for(i in 1:length(l2)){
    if(is.na(l2[i])) {
      predl2[i]=switch(mod,
                       #special VBGF
                       svbgf = l[i]+((exp(Xbeta[i])+0.5*totvar)-l[i]) * (1-exp(-K*dt[i])),
                       
                       #generalized VBGF
                       gvbgf = ((exp(Xbeta[i])+0.5*totvar)^(1/D)*(1-exp(-K*dt[i]))+l[i]^(1/D)*exp(-K*dt[i]))^D,
                       
                       #Gompertz
                       gomp = exp(log(exp(Xbeta[i])+0.5*totvar)*(1-exp(-K*dt[i]))+log(l[i])*exp(-K*dt[i]))
      )
    }
  }
  
  # get predicted lengths
  ADREPORT(predl2)
  ## return
  nllin
}

cmb = function(f, d, fl, mo) function(p) f(p, d, fl, mo)

fitlaa = function(mod = "gomp", formLinf = ~0, dat.train=dat.train) {
  X = model.matrix(formLinf, data=dat.train)
  parameters = list(
    beta = rep(0.1, ncol(X)), 
    logK = log(.1),
    logD = log(.5),
    logsdresid = 0,
    logsdcohort = 0,
    logsdyear = 0,
    bcohort = rep(0, length(unique(dat.train$cohort))), # random Linf by cohort
    byear = rep(0, length(unique(dat.train$year))) # random Linf by year
  )
 
  if(mod != "gvbgf"){
  obj = RTMB::MakeADFun(cmb(nll, dat.train, formLinf, mod), 
                        parameters, random=c("byear", "bcohort"), map = list(logD=factor(NA)))
  }else{
    obj = RTMB::MakeADFun(cmb(nll, dat.train, formLinf, mod), 
                          parameters, random=c("byear", "bcohort"))
  }
    
  opt = nlminb(obj$par, obj$fn, obj$gr)
  
  sdr = sdreport(obj)
  l=opt[["objective"]]
  k=length(opt[["par"]])
  AICc = 2*k + 2*l + 2*k*(k+1)/(nrow(dat.train)-k-1)
  if(!sdr$pdHess) AICc =NA
    
  return(list("obj"=obj, "opt"=opt, "parameters"=parameters, 
              "mod"=mod, "formLinf" = formLinf, "AICc"=AICc, 
  						"dat.train"=dat.train)) #object
  
}

predlaa = function(object, newdata, std.err=FALSE) {
  obj = object$obj
  opt = object$opt
  parameters = object$parameters
  mod = object$mod
  formLinf = object$formLinf
  dat.train = object$dat.train
  
  if(!all.equal(names(dat.train), names(newdata))) stop("newdata must have same columns as data used in fitting")
  #FIXME select columns of dat.train to make augdata (below) more flexibly
  
  oldPar = opt$par
  
  H = optimHess(oldPar, obj$fn, obj$gr)#just to save time
  
  sdrep = sdreport(obj, hessian.fixed=H) # extract estimates and Hessian-based and delta approximated sd
  pl = as.list(sdrep, what="Est")     # format parameter estimates as list of initials
  
  pred_cohort = length(parameters$bcohort)+1
  pred_year = length(parameters$byear)+1
  
  newdata$year = pred_year #overwrite year to get generic year predictions
  newdata$cohort = pred_cohort #overwrite year to get generic cohort predictions
  newdata$l2 = NA # to remove it from the likelihood
  newdata$pred = 1 #indicator 
  dat.train$pred = 0 #indicator 
  augdata = rbind(dat.train, newdata)
  
  newParameters = parameters
  newParameters$bcohort = pl$bcohort
  newParameters$byear = pl$byear 
  
  newParameters$bcohort[pred_cohort] = 0 #population level predictions
  newParameters$byear[pred_year] = 0 #population  level predictions
  
  mapArg = list()
  mapArg$bcohort = factor(rep(NA,length(newParameters$bcohort)))
  mapArg$byear = factor(rep(NA,length(newParameters$byear)))
  
  if(mod != "gvbgf"){
    newObj = RTMB::MakeADFun(cmb(nll, augdata, formLinf, mod), 
                             newParameters, random=c("bcohort", "byear"), map = list(logD=factor(NA)))
  }else{
    newObj = RTMB::MakeADFun(cmb(nll, augdata, formLinf, mod), 
                             newParameters, random=c("bcohort", "byear"))
  }
  
  newObj$fn(oldPar)  ## call once to update internal structures
  
  sdr = sdreport(newObj, oldPar, hessian.fixed=H, getReportCovariance=TRUE)
  augdata$pred.l2 = as.list(sdr, "Est", report=TRUE)$predl2
  if(std.err){
    augdata$std.err.l2 = as.list(sdr, "St", report=TRUE)$predl2
    return(subset(augdata, pred==1)[,c("pred.l2", "std.err.l2")])
  }#else
  return(subset(augdata, pred==1)[,c("pred.l2")])
}

selectlaa = function(dat.train, forms, mods) {
  grid = expand.grid(f = 1:length(forms), m=1:length(mods))
  tmp = foreach(i=1:nrow(grid)) %dopar% fitlaa(formLinf = as.formula(forms[[grid[i, "f"]]]), mod = mods[grid[i, "m"]], dat.train = dat.train)$AICc[1]
  tmp2 = do.call(c, tmp)
  grid$AICc = tmp2
  i=which.min(tmp2)
  formLinf = as.formula(forms[[mod=grid[i, "f"]]])
  mod = mods[grid[i, "m"]]
  return(fitlaa(formLinf = formLinf, mod = mod, dat.train = dat.train))
}

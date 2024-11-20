## ----setup, include=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  cache = T,
  echo = T,
  fig.width = 6, fig.height = 5,
  fig.pos = "H"
)


## ----required_packages, message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------------------------------------
# install.packages(c("corrplot", "doRNG", "dplyr", "forecast", "GA", "ggplot2", "glmmTMB", "knitr", "memoise", "MuMIn", "parallel", "patchwork", "tidyr"))

library(corrplot)
library(doRNG)
library(dplyr)
library(forecast)
library(GA)
library(ggplot2)
library(glmmTMB)
library(knitr)
library(memoise)
library(MuMIn)
library(parallel)
library(patchwork)
library(tidyr)


## ----genCovs, fig.height = 6----------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(1235)
n <- 70
ncov <- 12
AR <- runif(ncov, min = 0, max = 0.9) # auto-regression coefficients

envvars <- paste0("x", formatC(x = seq(ncov), digits = 1, flag = 0))
tmp <- lapply(seq(ncov), FUN = function(i){
  # Generate the time series
  ts_data <- arima.sim(model = list(ar = AR[i]), n = n)
  return(ts_data)
})

tmp <- as.data.frame(matrix(unlist(tmp), nrow = n, ncol = ncov))
names(tmp) <- envvars
dat <- cbind(data.frame(t = seq(n)), tmp)

# add a trend to x01 and x02
dat$x01 <- dat$x01 + diff(range(dat$x01))*0.01*(dat$t-mean(dat$t)) # 1% of sd per year
dat$x02 <- dat$x02 + diff(range(dat$x02))*0.01*(dat$t-mean(dat$t)) # 1% of sd per year

# plot time series
df <- tidyr::pivot_longer(dat, cols = seq(ncol(dat))[-1])
ggplot(df) + aes(x = t, y = value, fill = name, color = name) + 
  geom_area(show.legend = F) + 
  facet_wrap(~ name, ncol = 2, scales = "free_y") + 
  theme_bw()


## ----genBioVars, fig.height=3---------------------------------------------------------------------------------------------------------------------------------------------------
# generate ssb (random walk)
set.seed(236)
dat$ssb <- 30000 - cumsum(rnorm(n, 0, 3000))

# generate recruitment (Ricker model + env mediation from x01 & x03)
dat$rec <- exp(log(5) + log(dat$ssb) - 2.5e-5*dat$ssb + 0.05*dat$x01 + 0.05*dat$x03)
dat$rec <- dat$rec * rlnorm(n, 0, 0.2) # additional random noise

p1 <- ggplot(dat) + aes(x = t, y = ssb) + 
  geom_line()
p2 <- ggplot(dat) + aes(x = t, y = rec) + 
  geom_line()
p3 <- ggplot(dat) + aes(x = ssb, y = rec) + 
  geom_point() +
  lims(x = c(0, NA), y = c(0, NA))

layout <- "
AC
BC
"
p <- p1 + p2 + p3 + plot_layout(design = layout, axis_titles = "collect") & 
  theme_bw()
p


## ----corrPlot-------------------------------------------------------------------------------------------------------------------------------------------------------------------
corrplot.mixed(cor(dat[,-1]))


## ----fullModel------------------------------------------------------------------------------------------------------------------------------------------------------------------
fmla <- formula(paste(c("rec ~ offset(log(ssb)) + ssb", envvars), collapse = " + "))
fmla
fit0 <- glm(fmla, data = dat, family = gaussian(link = "log"))
summary(fit0)


## ----stepwiseRemoval------------------------------------------------------------------------------------------------------------------------------------------------------------
# model with all original variables in fit0
all_variables <- "rec ~ ."
# Variables that should always be in the model (Ricker without env.)
fixed_variables <- "rec ~ offset(log(ssb)) + ssb"  

fit <- step(fit0, 
  scope = list(
    lower = as.formula(fixed_variables), 
    upper = as.formula(all_variables)
  ), 
  direction = "both",
  trace = FALSE # set to TRUE if you would like to see the model selections steps
)

summary(fit)
c(fit0=AIC(fit0), fit=AIC(fit)) # comparison of AIC


## ----dredge---------------------------------------------------------------------------------------------------------------------------------------------------------------------
fmla <- formula(paste(c("rec ~ offset(log(ssb)) + ssb", envvars), collapse = " + "))
# argument na.action = na.fail required for fair comparisons
fittmp <- glm(fmla, data = dat, family = gaussian(link = "log"), na.action = na.fail)

t1 <- Sys.time()
dd <- dredge(fittmp, fixed = c("offset(log(ssb))", "ssb"), trace = FALSE)
t2 <- Sys.time()
# t2 - t1

# subset(dd[1:100,], delta < 4)
plot(dd[1:100,], labAsExpr = TRUE)


## ----plotRicker-----------------------------------------------------------------------------------------------------------------------------------------------------------------
envvarsKeep <- intersect(names(coef(fit)), envvars) # maintained env. covariates
# envvarsKeep <- c("x01", "x03")
focusVar <- envvarsKeep[1] # choose covariate of focus

## generate data for prediction
# focus covariate generated at 5 quantile levels, other env covariates use median 
L <- lapply(envvarsKeep, FUN = function(x){
  if(x == focusVar){probs = c(0.05,0.25,0.5,0.75,0.95)}else{probs <- 0.5}
  quantile(dat[,x], probs = probs)
})
names(L) <- envvarsKeep
L[["ssb"]] <-  seq(0, max(dat$ssb), len = 100)

nd <- expand.grid(L)
nd$rec <- predict(fit, newdata = nd, type = "response")

ggplot(dat) + aes(x = ssb, y = rec, col = get(focusVar)) + 
  geom_point() + 
  scale_color_continuous(name = focusVar, type = "viridis") +
  geom_line(data = nd, mapping = aes(group = get(focusVar))) + 
  theme_bw()


## ----metricFuns-----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Mean Squared Error (MSE)
mse <- function(y, yhat){mean((y - yhat)^2, na.rm = TRUE)}

# Root Mean Squared Error (RMSE)
rmse <- function(y, yhat){sqrt(mean((y - yhat)^2, na.rm = TRUE))}

# Mean Absolute Percentage Error (MAPE)
mape <- function(y, yhat){mean(abs(y - yhat) / y, na.rm = TRUE)}

# Median Absolute Percentage Error (MAPE)
mdape <- function(y, yhat){median(abs(y - yhat) / y, na.rm = TRUE)}

# Mean Absolute Error (MAE)
mae <- function(y, yhat){mean(abs(y - yhat), na.rm = TRUE)}

# Mean Absolute Scaled Error (MASE)
mase <- function(y, yhat, yhat_naive){mean(abs(y - yhat), na.rm = TRUE) /
    mean(abs(y - yhat_naive), na.rm = TRUE)}


## ----cvFuncions-----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Determines k-fold partitions for a given number of samples
# n is the number of samples; k is the number of partitions
kfold <- function(n, k = NULL){
  if(is.null(k)){k <- n} # if undefined, assume leave-one-out (LOO) CV
  fold <- vector(mode="list", k)
  n.remain <- seq(n)
  for(i in seq(k)){
    samp <- sample(seq(length(n.remain)), ceiling(length(n.remain)/(k-i+1)))
    fold[[i]] <- n.remain[samp]
    n.remain <- n.remain[-samp]
  }
  return(fold)
}

# perform repeated k-fold cross validation
# n and k are passed to kfold; seed is a integer value passed to set.seed
cvFun <- function(model = NULL, nperm = 10, k = 5, seed = 1){
  modelData <- model.frame(model)
  y <- model.response(model.frame(model))
  fmla <- formula(model)
  # when k = number of sample, a single leave-one-out cross validation is performed (LOOCV)
  if(k == length(y)){nperm <- 1} 
  res <- vector("list", nperm)
  set.seed(seed)
  for(j in seq(nperm)){ # for each permutation j
    fold <- kfold(n = nrow(modelData), k = k)
    res[[j]] <- data.frame(perm = j, fold = NaN, y = y, yhat = NaN)
    for(i in seq(fold)){ # for each fold
      train <- -fold[[i]]
      valid <- fold[[i]]
      # update model with training data
      fit.fold <- update(object = model, data = modelData[train,], 
        formula = as.formula(fmla, env = environment())) 
      # predict validation data
      res[[j]]$yhat[valid] <- predict(object = fit.fold, 
        newdata = modelData[valid,], type = "response")
      res[[j]]$fold[valid] <- i
    }
  }
  res <- do.call("rbind", res)
  return(res)
}


## ----kfoldExample, echo=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------
fold <- kfold(n = nrow(dat), k = 5)
df <- as.data.frame(do.call("cbind", fold))
names(df) <- seq(ncol(df))

df <- pivot_longer(df, cols = 1:5)
df <- df[order(df$value),]
names(df) <- c("fold", "t")

df = merge(dat, df, all = TRUE)

newdat <- data.frame(ssb = seq(0, max(df$ssb), len = 100))
pred <- vector("list", length(fold))
for (i in seq_along(fold)) {
  train <- -fold[[i]]
  valid <- fold[[i]]
  fit.fold <- update(fit0, data = df[train,], formula = as.formula(fixed_variables, env = environment()))
  newdat.i <- newdat
  newdat.i$fold <- as.character(i)
  newdat.i$rec <- predict(fit.fold, newdata = newdat.i, type = "response")
  pred[[i]] <- newdat.i
}
pred <- do.call("rbind", pred)

p1 <- ggplot(df) + aes(x = t, y = 0, col = fold) + 
  geom_point() +
  scale_color_brewer(name = "Fold", palette = "Set1") +
  labs(y = "") +
  scale_y_continuous(limits = c(0, 0), breaks = NULL) + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), axis.line.y = element_blank())

p2 <- ggplot(df) + aes(x = ssb, y = rec, col = fold) + 
  geom_point() +
  geom_line(data = pred, show.legend = F) +
  lims(x = c(0,NA), y = c(0,NA)) +
  scale_color_brewer(name = "Fold", palette = "Set1")

layout = "
AA
BB
BB
"

p <- (p1 / p2) + plot_layout(design = layout, guides = "collect") &
  theme_bw()
print(p)


## ----kfoldErrorEst, echo=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------
n <- 20
set.seed(123)
x <- runif(n, min = 0, max = 40)
err <- rnorm(n, sd = 15)
y <- 10 + 2*x + err
df <- data.frame(x, y)
fitlm <- lm(y ~ x, df)

res <- data.frame(k = c(round(exp(seq(log(2),log(n/2), len = 10))),n), rmse = NaN, se = NaN)
for(i in seq(nrow(res))){
  tmp <- cvFun(model = fitlm, nperm = 50, k = res$k[i], seed = 1)
  res$rmse[i] <- rmse(y = tmp$y, yhat = tmp$yhat) # sd of residuals
  byPerm <- tmp %>% group_by(perm) %>% summarise(rmse = rmse(y = y, yhat = yhat))
  res$se[i] <- sd(byPerm$rmse)/sqrt(nrow(byPerm))
  if(res$k[i]==n){ res$se[i] <- 0}
}

p1 <- ggplot(df) + aes(x = x, y = y) + 
  geom_point() + 
  geom_abline(intercept = coef(fitlm)[1], slope = coef(fitlm)[2])

p2 <- ggplot(res) + aes(x = k, y = rmse) + 
  geom_line() + 
  geom_point() + 
  geom_segment(mapping = aes(yend = rmse - se)) +
  geom_segment(mapping = aes(yend = rmse + se)) +
  geom_hline(yintercept = sqrt(mean(fitlm$residuals^2)), linetype = 2, col = 2) +
  geom_hline(yintercept = sqrt(mean(err^2)), linetype = 2, col = 4) + 
  geom_text(x = n/2, y = sqrt(mean(fitlm$residuals^2)), 
    label = "Model error", col = 2, vjust = -0.25) + 
  geom_text(x = n/2, y = sqrt(mean(err^2)), 
    label = "True error", col = 4, vjust = -0.25)

p <- (p1 | p2) & theme_bw()
print(p)
  


## ----cvModelCompare, warning=FALSE----------------------------------------------------------------------------------------------------------------------------------------------
res0 <- cvFun(model = fit0, nperm = 5, k = 5, seed = 1111)
res <- cvFun(model = fit, nperm = 5, k = 5, seed = 1111)

# mixed model with autoregressive term on time (t)
fmlatmp <- formula(fit)
fmlaMM <- as.formula(paste(c(deparse(fmlatmp), "+ ar1(t + 0 | 1)")))
fitMM <- glmmTMB(fmlaMM, data = dat, family = gaussian(link = "log"))
resMM <- cvFun(model = fitMM, nperm = 5, k = 5, seed = 1111)

cvStats0 <- with(res0, c(
  mse = mse(y = y, yhat = yhat),
  rmse = rmse(y = y, yhat = yhat),
  mape = mape(y = y, yhat = yhat),
  mdape = mdape(y = y, yhat = yhat),
  mae = mae(y = y, yhat = yhat)
))

cvStats <- with(res, c(
  mse = mse(y = y, yhat = yhat),
  rmse = rmse(y = y, yhat = yhat),
  mape = mape(y = y, yhat = yhat),
  mdape = mdape(y = y, yhat = yhat),
  mae = mae(y = y, yhat = yhat)
))

cvStatsMM <- with(resMM, c(
  mse = mse(y = y, yhat = yhat),
  rmse = rmse(y = y, yhat = yhat),
  mape = mape(y = y, yhat = yhat),
  mdape = mdape(y = y, yhat = yhat),
  mae = mae(y = y, yhat = yhat)
))

tmp <- cbind(cvStats0, cvStats, cvStatsMM)
tmp <- as.data.frame(t(apply(tmp, 1, FUN = function(x){
  paste0(format(x, scientific = TRUE, digits = 3), 
    ifelse(x == min(x, na.rm = T), "*", " "))
})))

names(tmp) <- c("res0", "res", "resMM")
tmp


## ----fitnessFun-----------------------------------------------------------------------------------------------------------------------------------------------------------------
fitnessFun <- function(genes = NULL, data = NULL, nperm = 10, k = 5, 
  ffSeed = 1, fitnessOnly = TRUE){

  fmla <- as.formula(
    paste(c("rec ~ offset(log(ssb)) + ssb", envvars[as.logical(genes)]), collapse = " + "), 
    env = environment())
  fit <- glm(formula = fmla, data = data, family = gaussian(link = "log"))
  
  res <- cvFun(model = fit, nperm = nperm, k = k, seed = ffSeed)
  
  fitness <- -1 * mape(y = res$y, yhat = res$yhat) # to maximize, so reverse sign 
  if(!fitnessOnly){
    byGroup <- res %>% group_by(perm) %>% summarise(mape = mape(y = y, yhat = yhat))
    se <- sd(byGroup$mape)/sqrt(nrow(byGroup))
  }
  if(fitnessOnly){
    out <- fitness
  }else{
    out <- list(fitness = fitness, byGroup = byGroup, se = se)
  }
  return(out)
}


## ----startPop-------------------------------------------------------------------------------------------------------------------------------------------------------------------
popSize <- max(30, length(envvars)) # number of individuals

PARS <- seq(envvars)*0
startPop <- matrix(0, ncol = length(PARS), nrow = popSize)
diag(startPop[seq(PARS),]) <- 1
if(popSize > length(PARS)){
  hit <- seq(nrow(startPop))[-seq(PARS)]
  startPop[hit, ] <- t(apply(startPop[hit, ], 1, 
    FUN = function(x){x[sample(seq(PARS), size = 2)] <- 1; x}))
}

tmp <- startPop
tmp[] <- c("off", "on")[c(tmp)+1]
df <- as.data.frame(tmp)
names(df) <- envvars
df$ind <- seq(nrow(df))

df_long <- df %>%
  pivot_longer(cols = -ind, names_to = "variable", values_to = "value")
  
ggplot(df_long, aes(x = variable, y = ind, fill = factor(value))) +
  geom_tile(colour = 1) +
  scale_fill_manual(name = "Gene state", 
    values = c("black", "white"), breaks = c("on", "off")) +
  scale_y_reverse(expand = c(0, 0)) + 
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = "Genes (covariates)", y = "Individuals", 
    title = "Starting population") +
  theme_bw()



## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

tmp1 <- apply(startPop, 1, FUN = function(x){fitnessFun(x, data = dat)})
tmp2 <- apply(startPop, 1, FUN = function(x){paste(envvars[as.logical(x)], collapse = "+")})
tmp <- data.frame(ind = seq(tmp1), fitness = tmp1, covariates = tmp2)

ggplot(tmp) + aes(x = ind, y = fitness, label = covariates) + 
  geom_point() + 
  geom_text(angle = 90, hjust = -0.1, vjust = 0.5) + 
  scale_y_continuous(expand = expansion(mult = c(0.05,0.25))) + 
  theme_bw()



## ----message=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------
# extra diagnostic data object to save all generations
postfit <- function(object, ...){
  pop <- as.data.frame(object@population)
  names(pop) <- paste0("par", seq(ncol(pop)))
  pop$fitness <- object@fitness
  pop$iter <- object@iter
  # update info
  if(!exists(".pop", envir = globalenv()))
    assign(".pop", NULL, envir = globalenv())
  .pop <- get(".pop", envir = globalenv())
  assign(".pop", rbind(.pop, pop), envir = globalenv()) 
  # output the input ga object (this is needed!)
  object 
}

mfitnessFun <- memoise(fitnessFun) # a memoised version of the fitness function
.pop <- NULL # make sure diagnostic output is empty
t1 <- Sys.time()
ga.fit <- ga(
  type = "binary", # for covariate selection
  fitness = mfitnessFun, # fitness function
  nBits = length(envvars), # number of genes (covariates)
  suggestions = startPop, # starting population suggestions
  popSize = popSize, # population size
  pcrossover = 0.8, # probability of gene crossover
  pmutation = 0.2,  # probability of mutation
  elitism = popSize*0.1, # number of best fitness individuals to survive at each generation
  parallel = FALSE, # set to TRUE for parallel computing (memoisation may be less effective)
  maxiter = 100, # max number of generations
  run = 40, # the max number of generations without any improvement
  monitor = FALSE, # set to TRUE to monitor fitness evolution
  data = dat, # input data
  seed = 1, # for reproducibility
  postFitness = postfit # extra diagnostics function
)
t2 <- Sys.time()
# t2 - t1 # elapsed time
tmp <- forget(mfitnessFun) # to stop memoisation

plot(ga.fit) # plot fitness development


## ----paretoPlot-----------------------------------------------------------------------------------------------------------------------------------------------------------------
.pop$nVar <- rowSums(.pop[,seq(envvars)])
.pop$fracVar <- .pop$nVar / length(envvars)
names(.pop)[seq(envvars)] <- envvars

# pareto front
agg <- aggregate(fitness ~ fracVar, data = .pop, FUN = max)

# determine standard error of pareto front models
parFront <- merge(x = agg, y = .pop, all.x = T)
parFront <- unique(parFront[,c(envvars, "nVar", "fracVar", "fitness")])
parFront <- parFront[order(parFront$nVar),]
parFront

parFront$se <- NaN
for(i in seq(nrow(parFront))){
  genes.i <- parFront[i,envvars]
  tmp <- fitnessFun(genes = genes.i, data = dat, nperm = 10, k = 5, ffSeed = 1, fitnessOnly = F)
  parFront$se[i] <- tmp$se
}


upr <- parFront$fitness + parFront$se
lwr <- parFront$fitness - parFront$se
mid <- parFront$fitness
# check if model fitness is greater than the upper bound of the simpler model
better <- mid > c(-Inf, upr[-length(lwr)]) 
parFront$best <- FALSE
hit <- max(which(better))
parFront$best[hit] <- TRUE

ggplot(.pop) + aes(x = -fitness, y = fracVar, color = iter) + 
  geom_point(size = 3) + 
  # scale_color_viridis_c(name = "Iteration", option = "magma", direction = -1) +
  scale_color_continuous(name = "Iteration", type = "viridis", direction = -1) +
  geom_path(data = parFront[better,], color = 4, linetype = 2) +
  geom_label(data = parFront[1,], color = 1, label = "Ricker", 
    hjust = 0.5, vjust = -0.5) +
  geom_segment(data = parFront, aes(x = -fitness-se, xend = -fitness+se, y = fracVar, yend = fracVar), 
    color = 4, linewidth = 1) +
  geom_text(data = subset(parFront, best), color = "red", label = "*", size = 10, vjust = 0.25) +
  labs(x = "MAPE", y = "Fraction of environmental covariates", title = "Pareto front") + 
  theme_bw()



## ----paretoFront----------------------------------------------------------------------------------------------------------------------------------------------------------------
parFront


## ----bestModel, include=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------
best <- subset(parFront, best)[envvars]
best <- paste(c(fixed_variables, names(best)[best == 1]), collapse = " + ")


## ----acfRec---------------------------------------------------------------------------------------------------------------------------------------------------------------------
acf(dat$rec)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
fitRicker <- glm(rec ~ offset(log(ssb)) + ssb, data = dat, gaussian(link = "log"))
acf(resid(fitRicker))


## ----randomBlocks, echo=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------
minSize <- 6

blockfold <- function(n, minSize){
  nBlock <- n %/% minSize # number of blocks
  leftover <- n %% minSize # remainder
  blockSize <- rep(minSize, nBlock) # define size of each block
  if(leftover > 0){ # add any leftovers
    blockSize <- blockSize + c(rep(1, leftover), rep(0, length(blockSize)-leftover))
  }
  blockSize <- sample(blockSize, length(blockSize)) # random shuffle
  idx <- inverse.rle(list(values = seq(blockSize), lengths = blockSize)) # fold index
  fold <- split(x = seq(n), idx) # create folds
  return(fold)
}

fold <- blockfold(n = nrow(dat), minSize = minSize)
df <- data.frame(t = dat$t, 
  fold = inverse.rle(list(values = seq(fold), lengths = sapply(fold, length))))

df = merge(dat, df, all = TRUE)

newdat <- data.frame(ssb = seq(0, max(df$ssb), len = 100))
pred <- vector("list", length(fold))
for (i in seq_along(fold)) {
  train <- -fold[[i]]
  valid <- fold[[i]]
  fit.fold <- update(fit0, data = df[train,], formula = as.formula(fixed_variables, env = environment()))
  newdat.i <- newdat
  newdat.i$fold <- as.character(i)
  newdat.i$rec <- predict(fit.fold, newdata = newdat.i, type = "response")
  pred[[i]] <- newdat.i
}
pred <- do.call("rbind", pred)

p1 <- ggplot(df) + aes(x = t, y = 0, col = factor(fold)) + 
  geom_point() +
  scale_color_brewer(name = "Fold", palette = "Paired") +
  labs(y = "") +
  scale_y_continuous(limits = c(0, 0), breaks = NULL) + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), axis.line.y = element_blank())

p2 <- ggplot(df) + aes(x = ssb, y = rec, col = factor(fold)) + 
  geom_point() +
  geom_line(data = pred, show.legend = F) +
  lims(x = c(0,NA), y = c(0,NA)) +
  scale_color_brewer(name = "Fold", palette = "Paired")

layout = "
AA
BB
BB
"

p <- (p1 / p2) + plot_layout(design = layout, guides = "collect") &
  theme_bw()
print(p)




## ----cvStatsBlock, eval=FALSE, include=FALSE------------------------------------------------------------------------------------------------------------------------------------
## 
## cvFunBlock <- function(model = NULL, minSize = 5, seed = 1){
##   modelData <- model.frame(model)
##   y <- model.response(model.frame(model))
##   fmla <- formula(model)
##   set.seed(seed)
##   fold <- blockfold(n = nrow(modelData), minSize = minSize)
##   res <- data.frame(fold = NaN, y = y, yhat = NaN)
##   for(i in seq(fold)){ # for each fold
##     train <- -fold[[i]]
##     valid <- fold[[i]]
##     # update model with training data
##     fit.fold <- update(object = model, data = modelData[train,],
##       formula = as.formula(fmla, env = environment()))
##     # predict validation data
##     res$yhat[valid] <- predict(object = fit.fold,
##       newdata = modelData[valid,], type = "response")
##     res$fold[valid] <- i
##   }
##   return(res)
## }
## 
## res0 <- cvFunBlock(model = fit0, minSize = 6, seed = 1)
## res <- cvFunBlock(model = fit, minSize = 6, seed = 1)
## 
## 
## cvStats0 <- with(res0, c(
##   mse = mse(y = y, yhat = yhat),
##   rmse = rmse(y = y, yhat = yhat),
##   mape = mape(y = y, yhat = yhat),
##   mdape = mdape(y = y, yhat = yhat),
##   mae = mae(y = y, yhat = yhat)
## ))
## 
## cvStats <- with(res, c(
##   mse = mse(y = y, yhat = yhat),
##   rmse = rmse(y = y, yhat = yhat),
##   mape = mape(y = y, yhat = yhat),
##   mdape = mdape(y = y, yhat = yhat),
##   mae = mae(y = y, yhat = yhat)
## ))
## 
## tmp <- cbind(cvStats0, cvStats)
## tmp <- as.data.frame(t(apply(tmp, 1, FUN = function(x){
##   paste0(format(x, scientific = TRUE, digits = 3),
##     ifelse(x == min(x, na.rm = T), "*", " "))
## })))
## 
## names(tmp) <- c("res0", "res")
## tmp
## 


## ----testHoldout, echo=FALSE, fig.height=2--------------------------------------------------------------------------------------------------------------------------------------
n <- nrow(dat) # length of time series
nTest <- 5 # number of test holdout values
nValid <- 15 # number of test holdout values
nTrain <- n - nValid - nTest # number of training values
df <- dat
df$split <- c(rep("Training", nTrain), rep("Validation", nValid), rep("Testing", nTest))
df$split <- factor(df$split, levels = c("Training", "Validation", "Testing"))
ggplot(df) + aes(x = t, y = 0, col = split) + 
  geom_point() +
  scale_color_brewer(name = "Split", palette = "Set2") +
  labs(y = "") +
  scale_y_continuous(limits = c(0, 0), breaks = NULL) + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
  theme_bw()


## ----rollingTestHoldout, echo=FALSE, fig.height=2-------------------------------------------------------------------------------------------------------------------------------
n <- nrow(dat) # length of time series
nTest <- 5 # number of test holdout values
nValid <- 10 # number of test holdout values
nTrain <- 45 # number of training values
res <- vector("list", n)
for(i in seq(n)){
  startIdx <- i
  if((startIdx + nTrain + nValid + nTest - 1) < n){
    trainIdx <- seq(startIdx, len = nTrain)
    validIdx <- seq(startIdx + nTrain, len = nValid)
    testIdx <- seq(startIdx + nTrain + nValid, len = nTest)
    
    res[[i]] <- data.frame(split = c(rep("Training", nTrain), rep("Validation", nValid), rep("Testing", nTest)),
      t = c(trainIdx, validIdx, testIdx), 
      fold = i)
  }
}
df <- do.call("rbind", res)
df$split <- factor(df$split, levels = c("Training", "Validation", "Testing"))

ggplot(df) + aes(x = t, y = fold, col = split) + 
  geom_point() +
  scale_color_brewer(name = "Split", palette = "Set2") +
  scale_y_reverse() + 
  theme_bw()


## ----kfoldHindcast, echo=FALSE, fig.height=2------------------------------------------------------------------------------------------------------------------------------------
n <- nrow(dat) # length of time series
k = 5
nTest <- 5 # number of test holdout values
nTrainValid <- n - nTest
res <- vector("list", 5)
ks <- kfold(n = n-nTest, k = 5)
tmp <- vector("list", k)
for(i in seq(ks)){
  df1 <- dat[seq(nTrainValid), c("t","rec")]
  df2 <- tail(dat[, c("t","rec")], nTest)
  df1$split <- NA
  df1$split[-ks[[i]]] <- "Training" 
  df1$split[ks[[i]]] <- "Validation"
  df2$split <- "Testing"
  df <- rbind(df1, df2)
  df$fold <- i
  tmp[[i]] <- df
}
df <- do.call("rbind", tmp)

df$split <- factor(df$split, levels = c("Training", "Validation", "Testing"))
# df$fold <- factor(df$fold)
ggplot(df) + aes(x = t, y = fold, color = split) + 
  geom_point() +
  scale_color_brewer(name = "Split", palette = "Set2") +
  scale_y_reverse() + 
  theme_bw()


## ----hindcast-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# split data into training/validation & testing
nTest <- 5
datTrainValid <- subset(dat, t <= max(t)-nTest)
datTest <- subset(dat, t > max(t)-nTest)

mfitnessFun <- memoise(fitnessFun) # a memoised version of the fitness function
.pop <- NULL # make sure diagnostic output is empty
t1 <- Sys.time()
ga.fit <- ga(
  type = "binary", # for covariate selection
  fitness = mfitnessFun, # fitness function
  nBits = length(envvars), # number of genes (covariates)
  suggestions = startPop, # starting population suggestions
  popSize = popSize, # population size
  pcrossover = 0.8, # probability of gene crossover
  pmutation = 0.2,  # probability of mutation
  elitism = popSize*0.1, # number of best fitness individuals to survive at each generation
  parallel = FALSE, # set to TRUE for parallel computing (memoisation may be less effective)
  maxiter = 100, # max number of generations
  run = 40, # the max number of generations without any improvement
  monitor = FALSE, # set to TRUE to monitor fitness evolution
  data = datTrainValid, # input data *** ONLY TRAIN/VALID ***
  seed = 1, # for reproducibility
  postFitness = postfit # extra diagnostics function
)
t2 <- Sys.time()
# t2 - t1 # elapsed time
tmp <- forget(mfitnessFun) # to stop memoisation

.pop$nVar <- rowSums(.pop[,seq(envvars)])
.pop$fracVar <- .pop$nVar / length(envvars)
names(.pop)[seq(envvars)] <- envvars

# pareto front
agg <- aggregate(fitness ~ fracVar, data = .pop, FUN = max)

# determine standard error of pareto front models
parFront <- merge(x = agg, y = .pop, all.x = T)
parFront <- unique(parFront[,c(envvars, "nVar", "fracVar", "fitness")])
parFront <- parFront[order(parFront$nVar),]
parFront

parFront$se <- NaN
for(i in seq(nrow(parFront))){
  genes.i <- parFront[i,envvars]
  tmp <- fitnessFun(genes = genes.i, data = dat, nperm = 10, k = 5, ffSeed = 1, fitnessOnly = F)
  parFront$se[i] <- tmp$se
}


upr <- parFront$fitness + parFront$se
lwr <- parFront$fitness - parFront$se
mid <- parFront$fitness
# check if model fitness is greater than the upper bound of the simpler model
better <- mid > c(-Inf, upr[-length(lwr)]) 
parFront$best <- FALSE
hit <- max(which(better))
parFront$best[hit] <- TRUE
best <- subset(parFront, best)[envvars]
best <- paste(c(fixed_variables, names(best)[best == 1]), collapse = " + ")


## ----naiveCompare---------------------------------------------------------------------------------------------------------------------------------------------------------------
# refit best model
fmla <- best
fitBest <- update(fit, data = datTrainValid, 
  formula = as.formula(fmla, env = environment()))
yhat <- predict(fitBest, newdata = datTest, type = "response")

# compare to geometric mean
yhat_naive_gm <- rep(exp(mean(log(tail(datTrainValid$rec, 10)))), nrow(datTest))
MASE_gm <- mase(y = datTest$rec, yhat = yhat, yhat_naive = yhat_naive_gm)

# compare to non-env Ricker
fmla <- "rec ~ ssb + offset(log(ssb))"
fitRicker <- update(fit, data = datTrainValid, 
  formula = as.formula(fmla, env = environment()))
yhat_naive_ri <- predict(fitRicker, newdata = datTest, type = "response")
MASE_ri <- mase(y = datTest$rec, yhat = yhat, yhat_naive = yhat_naive_ri)


## ----naiveCompare2, echo=FALSE, fig.height=4------------------------------------------------------------------------------------------------------------------------------------
df1 <- dat[,c("t", "rec")]
df1$type = ifelse(df1$t %in% datTrainValid$t, "train/valid", "test")
df2 <- datTest[,c("t", "rec")]
df2$rec <- yhat
df2$type = "model"
df3 <- datTest[,c("t", "rec")]
df3$rec <- yhat_naive_gm
df3$type = "Geom. mean (naive)"
df4 <- datTest[,c("t", "rec")]
df4$rec <- yhat_naive_ri
df4$type = "Ricker (naive)"

df <- rbind(df1, df2, df3, df4)
df$type <- factor(df$type, levels = c("train/valid", "test", "model", "Ricker (naive)", "Geom. mean (naive)"))

ggplot(df) + aes(x = t, y = rec, colour = type) + 
  geom_point() + 
  geom_line() +
  annotate("text", x = min(df$t), y = max(df$rec), 
    label = paste("MASE (Ricker) =", round(MASE_ri, 3)), 
    hjust = 0, color = "orange3") + 
  annotate("text", x = min(df$t), y = max(df$rec), 
    label = paste("MASE (Geom. mean) =", round(MASE_gm, 3)), 
    hjust = 0, vjust = 2, color = "green4") + 
  scale_color_manual(name = "Data type", 
    values = c("grey30", "grey", "blue", "orange3", "green4")) +
  theme_bw()



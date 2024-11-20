## ----setup, include=FALSE------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  cache = T,
  echo = T,
  fig.width = 6, fig.height = 5,
  fig.pos = "H"
)


## ----required_packages, message=FALSE, warning=FALSE---------------------------------------------------------------------------------------------
library(parallel)
library(forecast)
library(tidyr)
library(dplyr)
library(ggplot2)
library(corrplot)
library(GA)
library(patchwork)
library(doRNG)
library(knitr)
library(glmmTMB)


## ----genCovs-------------------------------------------------------------------------------------------------------------------------------------
set.seed(1235)
n <- 70
ncov <- 12
AR <- runif(ncov, min = 0, max = 0.9) # auto-regression coefficients
MA <- runif(ncov, min = 0.1, max = 0.5) # moving average coefficients
MA[sample(ncov, size = ncov*0.5)] <- 0 # set some MA to zero for AR1 process only 

envvars <- paste0("x", formatC(x = seq(ncov), digits = 1, flag = 0))
tmp <- lapply(seq(ncov), FUN = function(i){
  # Generate the time series
  ts_data <- arima.sim(model = list(ar = AR[1], ma = MA[i]), n = n)
  return(ts_data)
})

tmp <- as.data.frame(matrix(unlist(tmp), nrow = n, ncol = ncov))
names(tmp) <- envvars
dat <- cbind(data.frame(t = seq(n)), tmp)

# add a trend to x01 and x02
dat$x01 <- dat$x01 + sd(dat$x01)*0.02*(dat$t-mean(dat$t)) # 2% of sd per year
dat$x02 <- dat$x02 + sd(dat$x02)*0.03*(dat$t-mean(dat$t)) # 3% of sd per year

# plot time series
df <- tidyr::pivot_longer(dat, cols = seq(ncol(dat))[-1])
ggplot(df) + aes(x = t, y = value, fill = name, color = name) + 
  geom_area(show.legend = F) + 
  facet_wrap(~ name, ncol = 2)




## ----genBio, fig.height=3------------------------------------------------------------------------------------------------------------------------
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

p <- p1 + p2 + p3 + plot_layout(design = layout, axis_titles = "collect") 
p


## ----corrPlot------------------------------------------------------------------------------------------------------------------------------------
corrplot.mixed(cor(dat[,-1]))


## ----fullModel-----------------------------------------------------------------------------------------------------------------------------------
fmla <- formula(paste(c("rec ~ offset(log(ssb)) + ssb", envvars), collapse = " + "))
fmla
fit0 <- glm(fmla, data = dat, family = gaussian(link = "log"))
summary(fit0)


## ----stepwiseRemoval-----------------------------------------------------------------------------------------------------------------------------
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

# compare the AIC of 'true' model
# AIC(update(fit, formula("~ . + x03")))


## ----plotRicker----------------------------------------------------------------------------------------------------------------------------------
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
  geom_line(data = nd, mapping = aes(group = get(focusVar)))


## ----metricFuns----------------------------------------------------------------------------------------------------------------------------------
# Mean Squared Error (MSE)
mse <- function(y, yhat){
  mean((y - yhat)^2, na.rm = TRUE)
}

# Root Mean Squared Error (RMSE)
rmse <- function(y, yhat){
  sqrt(mean((y - yhat)^2, na.rm = TRUE))
}

# Mean Absolute Percentage Error (MAPE)
mape <- function(y, yhat){
  mean(abs(y - yhat) / y, na.rm = TRUE)
}

# Median Absolute Percentage Error (MAPE)
mdape <- function(y, yhat){
  median(abs(y - yhat) / y, na.rm = TRUE)
}

# Mean Absolute Error (MAE)
mae <- function(y, yhat){
  mean(abs(y - yhat), na.rm = TRUE)
}

# Mean Absolute Scaled Error (MASE)
mase <- function(y, yhat, yhat_naive){
  mean(abs(y - yhat), na.rm = TRUE) /
    mean(abs(y - yhat_naive), na.rm = TRUE)
}


## ----cvFuncions----------------------------------------------------------------------------------------------------------------------------------
# Determines k-fold partitions for a given number of samples
# n is the number of samples; k is the number of partitions
kfold <- function(n, k = NULL){
  if(is.null(k)){k <- n} # if undefined, assume leave-one-out (LOO) CV
  res <- vector(mode="list", k)
  n.remain <- seq(n)
  for(i in seq(k)){
    samp <- sample(seq(length(n.remain)), ceiling(length(n.remain)/(k-i+1)))
    res[[i]] <- n.remain[samp]
    n.remain <- n.remain[-samp]
  }
  return(res)
}

# perform repeated k-fold cross validation
# n and k are passed to kfold; seed is a integer value passed to set.seed
cvFun <- function(model = NULL, nperm = 10, k = 5, seed = 1){
  
  if(any(class(model) %in% "glmmTMB")){
    modelData <- model$frame
    y <- model.response(model$frame)
  }
  if(any(class(model) %in% "lm")){
    modelData <- model$data
    y <- model$y
  }
  res <- vector("list", nperm)
  set.seed(seed)
  for(j in seq(nperm)){ # for each permutation j
    ks <- kfold(n = nrow(modelData), k = k)
    res[[j]] <- data.frame(perm = j, k = NaN, y = y, yhat = NaN)
    for(i in seq(ks)){ # for each partition k
      train <- -ks[[i]]
      valid <- ks[[i]]
      fit.k <- update(object = model, data = modelData[train,]) # update model with training data
      res[[j]]$k[valid] <- i
      res[[j]]$yhat[valid] <- predict(object = fit.k, newdata = modelData[valid,], type = "response")
    }
  }
  res <- do.call("rbind", res)
  return(res)
}


## ----kfold---------------------------------------------------------------------------------------------------------------------------------------
ks <- kfold(n = nrow(dat), k = 5)
df <- as.data.frame(do.call("cbind", ks))
names(df) <- seq(ncol(df))

df <- pivot_longer(df, cols = 1:5)
df <- df[order(df$value),]

p1 <- ggplot(dat) + aes(x = t, y = 1, col = df$name) + 
  geom_point() +
  scale_color_brewer(palette = "Set1", name = "Fold") +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), axis.line.y = element_blank())

p2 <- ggplot(dat) + aes(x = ssb, y = rec, col = df$name) + 
  geom_point() +
  lims(x = c(0,NA), y = c(0,NA)) +
  scale_color_brewer(palette = "Set1", name = "Fold")

layout = "
AA
BB
BB
"

p <- p1 + p2 + plot_layout(design = layout, guides = "collect")
p



## ----cvModelCompare------------------------------------------------------------------------------------------------------------------------------
res0 <- cvFun(model = fit0, nperm = 5, k = 5, seed = 1111)
res <- cvFun(model = fit, nperm = 5, k = 5, seed = 1111)

# mixed model with autoregressive term on time (t)
fmlaMM <- formula(rec ~ ssb + offset(log(ssb)) + x01 + ar1(t + 0 | 1))
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


## ----fitnessFun----------------------------------------------------------------------------------------------------------------------------------
fitnessFun <- function(genes = NULL, data = NULL, covars = envvars,
  nperm = 10, k = 5, ffSeed = 1, fitnessOnly = TRUE){
  if(length(genes) != length(covars)){stop("length of 'genes' must equal length of 'covars'")}

  fmlaGA <<- paste(c("rec ~ offset(log(ssb)) + ssb", covars[as.logical(genes)]), collapse = " + ")

  fitGA <- glm(formula = fmlaGA, data = data, family = gaussian(link = "log"))
  
  res <- cvFun(model = fitGA, nperm = nperm, k = k, seed = ffSeed)
  
  fitness <- -1 * mape(y = res$y, yhat = res$yhat) # to maximize, so reverse sign 
  
  if(!fitnessOnly){
    byPerm <- res %>% group_by(perm) %>% summarise(mape = mape(y = y, yhat = yhat))
  }

  if(fitnessOnly){
    out <- fitness
  }else{
    out <- list(fitness = fitness, byPerm = byPerm)
  }
  
  return(out)
}


## ----startPop------------------------------------------------------------------------------------------------------------------------------------
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
  scale_fill_manual(name = "Gene state", values = c("black", "white"), breaks = c("on", "off")) +   # Use color scale (e.g., viridis)
  scale_y_reverse(expand = c(0, 0)) + 
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = "Genes (covariates)", y = "Individuals", title = "Starting population") +
  theme_bw()



## ------------------------------------------------------------------------------------------------------------------------------------------------

tmp1 <- apply(startPop, 1, FUN = function(x){fitnessFun(x, data = dat)})
tmp2 <- apply(startPop, 1, FUN = function(x){paste(envvars[as.logical(x)], collapse = "+")})
tmp <- data.frame(ind = seq(tmp1), fitness = tmp1, covariates = tmp2)

ggplot(tmp) + aes(x = ind, y = fitness, label = covariates) + 
  geom_point() + 
  geom_text(angle = 90, hjust = -0.1, vjust = 0.5) + 
  scale_y_continuous(expand = expansion(mult = c(0.05,0.25)))



## ----geneticAlg----------------------------------------------------------------------------------------------------------------------------------
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

.pop <- NULL # make sure diagnostic output is empty
ga.fit <- ga(
  type = "binary", # for covariate selection
  fitness = fitnessFun, # fitness function
  nBits = length(envvars), # number of genes (covariates)
  suggestions = startPop, # starting population suggestions
  popSize = popSize, # population size
  pcrossover = 0.8, pmutation = 0.1, elitism = popSize*0.05, # default settings
  parallel = round(parallel::detectCores()/2), # parallel computing using half of available cores
  maxiter = 100, # max number of generations
  run = 30, # the max number of generations without any improvement
  monitor = FALSE, # set to TRUE to monitor fitness evolution
  data = dat, # input data
  seed = 1, # for reproducibility
  postFitness = postfit # extra diagnostics function
)

plot(ga.fit) # plot development


## ----paretoPlot----------------------------------------------------------------------------------------------------------------------------------
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
  tmp <- fitnessFun(genes = genes.i, data = dat, nperm = 20, k = 5, ffSeed = 1, fitnessOnly = F)
  parFront$se[i] <- sd(tmp$byPerm$mape)/sqrt(nrow(tmp$byPerm))
}
parFront$fitnessUpper <- parFront$fitness + parFront$se

which(parFront$fitness > c(-Inf, (parFront$fitnessUpper)[-nrow(parFront)]))

parFront$best <- FALSE
hit <- max(which(parFront$fitness > c(-Inf, (parFront$fitness + parFront$se)[-nrow(parFront)])))
parFront$best[hit] <- TRUE

ggplot(.pop) + aes(x = -fitness, y = fracVar, color = iter) + 
  geom_point(size = 3) + 
  # scale_color_viridis_c(name = "Iteration", option = "magma", direction = -1) +
  scale_color_continuous(name = "Iteration", type = "viridis", direction = -1) +
  geom_path(data = parFront, color = 4, linetype = 2) +
  geom_label(data = parFront[1,], color = 1, label = "Naive model", 
    hjust = -0.25, vjust = 0.5) +
  geom_segment(data = parFront, aes(x = -fitness-se, xend = -fitness+se, y = fracVar, yend = fracVar), 
    color = 4, linewidth = 1) +
  labs(x = "MAPE", y = "Fraction of covariates", title = "Pareto front")



## ----paretoFront---------------------------------------------------------------------------------------------------------------------------------
parFront


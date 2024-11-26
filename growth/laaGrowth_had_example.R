library(RTMB)
library(doParallel)
library(ggplot2); theme_set(theme_bw())
source("growth/laaGrowthRTMB.R")

# had.27.46a20-----------------
load("growth/had.27.46a20.Rdata") #dat and LWab
#for new stocks, get LWab from https://www.fishbase.se/search.php

dat = transform(subset(dat, age!=max(dat$age)),
								l = 10^(log10(w*1000/LWab$a)/LWab$b), 
								l2 = 10^(log10(w2*1000/LWab$a)/LWab$b),
								ssb.s=scale(ssb), sal.s=scale(sal), temp.s=scale(temp))


# Try a gvbgf with formLinf = ~1 + ssb.s+ temp.s
gFit = fitlaa(dat.train=dat, formLinf = ~1 + ssb.s+ temp.s, mod = "gvbgf")

#make predictions back on the observed data
dat$predl2=predlaa(gFit, newdata = dat)

## transform back to w ---------
dat$predw2 = (LWab$a * dat$predl2 ^ LWab$b)/1000

ggplot(dat, aes(x=age, colour=as.factor(cohort)))+
	geom_point(aes(y=w2))+
	geom_line(aes(y=predw2), alpha=0.25)+
	guides(colour=guide_legend(title="Cohort"))+
	ylab("weight")

# Model selection---------------------

forms = list( #formulas for modifying Linf
  ~1 + ssb.s + sal.s + temp.s,
  ~1 + ssb.s + sal.s,
  ~1 + ssb.s + temp.s,
  ~1 + sal.s + temp.s,
  ~1 + ssb.s,
  ~1 + temp.s,
  ~1 #null
)
mods = c( "svbgf", "gvbgf", "gomp")


## select forms and mods -------------------

#set up cluster 
ncores <- detectCores() # max available cores
clusterType <- ifelse(Sys.info()["sysname"] == "Windows", "PSOCK", "FORK")
mycl <- parallel::makeCluster(ncores, type=clusterType)  
doParallel::registerDoParallel(mycl)

invisible(clusterEvalQ(mycl, {
	library(RTMB , logical = TRUE)
	source("growth/laaGrowthRTMB.R")
}))

grid = expand.grid(f = 1:length(forms), m=1:length(mods))
tmp = foreach(i=1:nrow(grid)) %dopar% fitlaa(formLinf = as.formula(forms[[grid[i, "f"]]]), mod = mods[grid[i, "m"]], dat.train = dat)$AICc[1]
tmp2 = do.call(c, tmp)
grid$AICc = tmp2
i=which.min(tmp2)
(formLinf = as.formula(forms[[mod=grid[i, "f"]]]))
(mod = mods[grid[i, "m"]])

## AICc table--------------
grid = transform(grid, form=as.character(forms[f]), 
               mod=mods[m],
               deltaAIC = AICc-grid$AICc[i])
grid[order(grid$AICc),]

## Predict from the best model---------
gFit = fitlaa(dat.train=dat, formLinf = formLinf, mod = mod)

#make predictions back on the observed data
dat$predl2=predlaa(gFit, newdata = dat)

## transform back to w ---------
dat$predw2 = (LWab$a * dat$predl2 ^ LWab$b)/1000

ggplot(dat, aes(x=age, colour=as.factor(cohort)))+
  geom_point(aes(y=w2))+
  geom_line(aes(y=predw2), alpha=0.25)+
  guides(colour=guide_legend(title="Cohort"))+
  ylab("weight")

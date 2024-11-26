library(glmmTMB)
library(MuMIn)
library(broom.mixed)
library(doParallel)

load("growth/had.27.46a20.Rdata") #dat

#transform the data 
dat = transform(subset(dat, age<max(dat$age)), # remove the plus group 
							logw2=log(w2), logw=log(w), # log weight variables
							ssb.s=scale(ssb), sal.s=scale(sal), temp.s=scale(temp), #scale environment
							yearf=factor(year), cohortf=factor(cohort)) #not necessary here

# Fit a single model-----

gFit1 = glmmTMB(logw2 ~ logw + age + 
								 	sal.s + ssb.s + temp.s+
								 	(1|cohortf) + (1|yearf), data=dat)

# Plot the model estimates with CI-----
t1 = broom.mixed::tidy(gFit1, effects = "fixed", conf.int = TRUE)
ggplot(subset(t1, term !="(Intercept)"), 
			 aes(estimate, term, xmin=conf.low, xmax=conf.high))+
	geom_errorbarh(height=0)+
	geom_vline(xintercept=0,lty=2)+
	geom_point()+ylab(NULL)


# Model selection---- 
#global model of weight-at-age with random effects
fullformRE = logw2 ~ logw + age+ logw:age +  I(age^2) + I(logw^2) +  logw:I(age^2)+
	sal.s + ssb.s + temp.s +I(temp.s^2) + 
	logw:sal.s + logw:ssb.s + logw:temp.s +(1|cohortf) + (1|yearf)

#fit the global model
global.model = glmmTMB(fullformRE,
												data=dat, na.action = na.fail)

#fit all subsets of fixed effects and make AICc table
#set up cluster 
ncores = detectCores()-1 # max available cores
clusterType = ifelse(Sys.info()["sysname"] == "Windows", "PSOCK", "FORK")
mycl = parallel::makeCluster(ncores, type=clusterType)  
doParallel::registerDoParallel(mycl)
invisible(clusterEvalQ(mycl, {
	library(glmmTMB , logical = TRUE)
}))

clusterExport(mycl, "dat")
dtab = dredge(global.model, rank = "AICc", cluster = mycl)
head(dtab)

gFit = get.models(dtab, subset= 1)[[1]]

# Plot best model with CI-----
t1 = broom.mixed::tidy(gFit, effects = "fixed", conf.int = TRUE)
ggplot(subset(t1, term !="(Intercept)"), 
			 aes(estimate, term, xmin=conf.low, xmax=conf.high))+
	geom_errorbarh(height=0)+
	geom_vline(xintercept=0,lty=2)+
	geom_point()+ylab(NULL)
#Be aware that these CI do not have 95% coverage because of model selection	

# Predicting back on the real scale---------------
newdata = dat

newdata$cohortf = NA #omit RE to get mean predictions
newdata$yearf = NA #omit RE to get mean predictions

newdata$predlogw2 = predict(gFit, newdata = newdata, allow.new.levels=TRUE) # predict size2
tot.var = as.numeric(VarCorr(gFit)$cond$cohortf+VarCorr(gFit)$cond$yearf+attr(VarCorr(gFit)$cond,"sc")^2)
#predict on real scale with all 3 variance components included
newdata$predw2 = exp(newdata$predlogw2 + tot.var/2) #lognormal mean

ggplot(newdata, aes(year+1, w2))+geom_point(alpha=0.5)+
			 	geom_line(data=newdata, aes(y=predw2))+
	facet_wrap(~age, labeller = label_both)+
	ylab("weight-at-age") 

ggplot(newdata, aes(w2, predw2))+geom_point(alpha=0.25)+
	geom_abline(intercept = 0, slope = 1)+
	facet_wrap(~age, labeller = label_both, scale="free")

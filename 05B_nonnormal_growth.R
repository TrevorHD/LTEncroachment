## Purpose: trying to model non-normal growth for creosote SIPM
## Author: Tom Miller
## Last modified: 4.22.2021

library(mgcv)
library(sgt)
library(maxLik)
library(zoo)
source("C:/Users/tm9/Dropbox/github/IPM_size_transitions/Diagnostics.R")
  # Create unique transect as interaction of transect and site
  # First-time only; this bit of code in 06_BootRes does it all other times
  # No transplants in this data set
LATR_full <- CData %>% 
  mutate(unique.transect = interaction(transect, site))

## set the gamma argument of gam() -- gamma>1 generates smoother fits, less likely to be overfit
gamma = 1.2

##### Growth model ----------------------------------------------------------------------------------------

# Prepare a data subset for growth that drops rows missing either t or t1 size data
# Also create log_volume as a new variable because GAM doesn't like functions of variables as variables
LATR_grow <- LATR_full  %>% drop_na(volume_t, volume_t1) %>%
  mutate(log_volume_t = log(volume_t),
         log_volume_t1 = log(volume_t1))

# Create empty list to populate with model results
LATR_gam_models <- list()

# Three candidate models for the mean: size only, size + density, or size, density, and size:density
# Three candidates for variance: size only, size + density, fitted value (all the covariates plus rfx)

# constant sigma
LATR_gam_models[[1]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(unique.transect, bs = "re"), ~1),
                            data = LATR_grow, gamma = gamma, family = gaulss())
LATR_gam_models[[2]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(weighted.dens) + s(unique.transect, bs = "re"), ~1), 
                            data = LATR_grow, gamma = gamma, family = gaulss())                
LATR_gam_models[[3]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(weighted.dens) + weighted.dens:log_volume_t + s(unique.transect, bs = "re"), ~1), 
                            data = LATR_grow, gamma = gamma, family = gaulss()) 
# sigma as function of initial size
LATR_gam_models[[4]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(unique.transect,bs = "re"), ~s(log_volume_t)), 
                            data = LATR_grow, gamma = gamma, family = gaulss())
LATR_gam_models[[5]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(weighted.dens) + s(unique.transect, bs = "re"), ~s(log_volume_t)), 
                            data = LATR_grow, gamma = gamma, family = gaulss())                
LATR_gam_models[[6]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(weighted.dens) + weighted.dens:log_volume_t + s(unique.transect,bs = "re"), ~s(log_volume_t)), 
                            data = LATR_grow, gamma = gamma, family = gaulss()) 
# sigma depends on both initial size and density
LATR_gam_models[[7]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(unique.transect, bs = "re"), ~s(log_volume_t) + s(weighted.dens)), 
                            data = LATR_grow, gamma = gamma, family = gaulss())
LATR_gam_models[[8]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(weighted.dens) + s(unique.transect, bs = "re"), ~s(log_volume_t) + s(weighted.dens)), 
                            data = LATR_grow, gamma = gamma, family = gaulss())                
LATR_gam_models[[9]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(weighted.dens) + weighted.dens:log_volume_t + s(unique.transect, bs = "re"), ~s(log_volume_t) + s(weighted.dens)), 
                            data = LATR_grow, gamma = gamma, family = gaulss()) 
# Collect model AICs into a single table
grow_aic <- AICtab(LATR_gam_models, base = TRUE, sort = FALSE)
## define best Gaussian model
LATR_gam_model <- LATR_gam_models[[which.min(grow_aic$AIC)]]
## extract fitted values
LATR_grow$fitted_vals = predict(LATR_gam_model,type="response",data=LATR_grow)
## extract the linear predictor for the mean
LATR_Xp <- predict.gam(LATR_gam_model,type="lpmatrix")
##################################################################  
# Extract values of the fitted splines to explore their properties 
##################################################################
fitted_terms = predict(LATR_gam_model,type="terms")  

##### effect of initial size on mean of final size 
plot(LATR_grow$log_volume_t, fitted_terms[,"s(log_volume_t)"]) 
##### effect of weighted.dens on mean of final size 
plot(LATR_grow$weighted.dens, fitted_terms[,"s(weighted.dens)"]); ## complicated 

##### sigma versus size - presumably log scale
plot(LATR_grow$log_volume_t, fitted_terms[,"s.1(log_volume_t)"]) 
##### sigma versus density 
plot(LATR_grow$weighted.dens, fitted_terms[,"s.1(weighted.dens)"]) 

##################################################################  
# Inspect scaled residuals to evaluate the pilot model: FAILS 
##################################################################
fitted_all = predict(LATR_gam_model,type="response")  
sigma.hat = 1/fitted_all[,2]
scaledResids = residuals(LATR_gam_model,type="response")/sigma.hat;  # note the 'type' argument is needed
par(mfrow=c(1,2))
plot(fitted_all[,1], scaledResids,xlab="Fitted values", ylab="Scaled residuals") 

qqPlot(scaledResids) # really bad in both tails
jarque.test(scaledResids) # normality test: FAILS, P < 0.001 
anscombe.test(scaledResids) # kurtosis: FAILS, P < 0.001 
agostino.test(scaledResids) # skewness: FAILS, P<0.001 

######## rolling NP moments diagnostic: skew is small but variable, tails are fat. 
px = LATR_grow$fitted_vals[,1]; py=scaledResids; 
z = rollMomentsNP(px,py,windows=8,smooth=TRUE,scaled=TRUE) 

###########################################################
## fit sgt
###########################################################
sgtLogLik=function(pars,response){
  val = dsgt(x = response, 
             mu=pars[1]+pars[2]*LATR_grow$log_volume_t+pars[3]*LATR_grow$weighted.dens,
             sigma=exp(pars[4]+pars[5]*LATR_grow$log_volume_t+pars[6]*LATR_grow$weighted.dens),
             lambda=pars[7],
             p=exp(pars[8]),
             q=exp(pars[9]),
             log=T) 
  return(val); 
}

## initial parameter values
p0=c(0,1,0,1,0,0,0,1,1) 
paranoid_iter <- 3
coefs = list(paranoid_iter); LL=numeric(paranoid_iter);  
for(j in 1:paranoid_iter) {
  out=maxLik(logLik=sgtLogLik,start=p0*exp(0.2*rnorm(length(p0))), response=LATR_grow$log_volume_t1,
             method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE); 
  
  out=maxLik(logLik=sgtLogLik,start=out$estimate,response=LATR_grow$log_volume_t1,
             method="NM",control=list(iterlim=5000,printLevel=1),finalHessian=FALSE); 
  
  out=maxLik(logLik=sgtLogLik,start=out$estimate,response=LATR_grow$log_volume_t1,
             method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=FALSE); 
  
  coefs[[j]] = out$estimate; LL[j] = out$maximum;
  cat(j, "#--------------------------------------#",out$maximum,"\n"); 
}

j = min(which(LL==max(LL))) 
out=maxLik(logLik=sgtLogLik,start=coefs[[j]],response=LATR_grow$log_volume_t1,
           method="BHHH",control=list(iterlim=5000,printLevel=2),finalHessian=TRUE) 

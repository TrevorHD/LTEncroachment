## Author: Tom Miller
## Purpose: This script builds a local (non-spatial) density-dependent IPM, which is the foundation for the SIPM
## This code requires objects create in script 04_CDataPrep.R (this could be sourced here, but for now I have just run that script manually)

library(tidyverse)
library(mgcv)
library(popbio)
library(bbmle)

# misc functions
invlogit <- function(x){exp(x)/(1+exp(x))}

## Full demography data frame
LATR_full <- CData %>% 
  mutate(
    #create unique transect as interaction of transect and site
    unique.transect = interaction(transect, site))

############################################################################
##1. Growth
############################################################################

## prep a data subset for growth that drops rows missing either t or t1 size data
LATR_grow <- LATR_full  %>% drop_na(volume_t,volume_t1) %>% 
  ## also create log_volume as a new variable because gam doesn't like functions of variables as variables
  mutate(log_volume_t = log(volume_t),
         log_volume_t1 = log(volume_t1))

## three candidate models for the mean: size only, size + density, or size, density, and size:density
## three candidates for variance: size only, size+density, fitted value (all the covariates plus rfx)
LATR_gam_models=list()
## Pilot fits, where sigma depends on initial size only
LATR_gam_models[[1]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(unique.transect,bs="re"), ~s(log_volume_t)), 
                            data=LATR_grow, gamma=1.4, family=gaulss())
LATR_gam_models[[2]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(weighted.dens) + s(unique.transect,bs="re"), ~s(log_volume_t)), 
                            data=LATR_grow, gamma=1.4, family=gaulss())                
LATR_gam_models[[3]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(weighted.dens) + weighted.dens:log_volume_t + s(unique.transect,bs="re"), ~s(log_volume_t)), 
                            data=LATR_grow, gamma=1.4, family=gaulss())  
## Fits where sigma depends on both initial size and density
LATR_gam_models[[4]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(unique.transect,bs="re"), ~s(log_volume_t) + s(weighted.dens)), 
                            data=LATR_grow, gamma=1.4, family=gaulss())
LATR_gam_models[[5]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(weighted.dens) + s(unique.transect,bs="re"), ~s(log_volume_t) + s(weighted.dens)), 
                            data=LATR_grow, gamma=1.4, family=gaulss())                
LATR_gam_models[[6]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(weighted.dens) + weighted.dens:log_volume_t + s(unique.transect,bs="re"), ~s(log_volume_t) + s(weighted.dens)), 
                            data=LATR_grow, gamma=1.4, family=gaulss()) 
## these models will be iterated to fit sigma as f(fitted value)
LATR_grow$fitted_vals = LATR_grow$log_volume_t 
LATR_gam_models[[7]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(unique.transect,bs="re"), ~s(fitted_vals)), 
                            data=LATR_grow, gamma=1.4, family=gaulss())
LATR_gam_models[[8]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(weighted.dens) + s(unique.transect,bs="re"), ~s(fitted_vals)), 
                            data=LATR_grow, gamma=1.4, family=gaulss())                
LATR_gam_models[[9]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(weighted.dens) + weighted.dens:log_volume_t + s(unique.transect,bs="re"), ~s(fitted_vals)), 
                            data=LATR_grow, gamma=1.4, family=gaulss())  
for(mod in 7:9) {
  fitGAU = LATR_gam_models[[mod]]
  fitted_all = predict(fitGAU,type="response",data=LATR);                  
  fitted_vals = new_fitted_vals = fitted_all[,1]; 
  weights = fitted_all[,2]; # what I call "weights" here are 1/sigma values; see ?gaulss for details.
  
  err=100; k=0; 
  while(err>10^(-6)) {
    LATR_grow$fitted_vals = new_fitted_vals; 
    fitGAU <- update(fitGAU); 
    fitted_all = predict(fitGAU,type="response",data=LATR_grow);   
    new_fitted_vals = fitted_all[,1]; new_weights = fitted_all[,2];
    err = weights - new_weights; err=sqrt(mean(err^2)); 
    weights = new_weights; 
    k=k+1; cat(k,err,"\n"); 
  }   
  LATR_gam_models[[mod]] =  fitGAU;
}

grow_aic <- AICtab(LATR_gam_models,base=T,sort=F) 
## Model 5 is the winner: mean ~ s(size) + s(density), sd ~ s(size) + s(density)
## define models 5 as our best Gaussian model
LATR_grow_best <- LATR_gam_models[[which.min(grow_aic$AIC)]]
LATR_grow_fitted_terms = predict(LATR_grow_best,type="terms") 
LATR_grow$pred = predict.gam(LATR_grow_best,newdata = LATR_grow, exclude = "s(unique.transect)")
##### effect of size on future size -- obviously linear
plot(LATR_grow$log_volume_t,LATR_grow_fitted_terms[,"s(log_volume_t)"]) 
#### effect of density on growth 
plot(LATR_grow$weighted.dens,LATR_grow_fitted_terms[,"s(weighted.dens)"]) 
## effect of size and density on sd(future size)
plot(LATR_grow$log_volume_t,LATR_grow_fitted_terms[,"s.1(log_volume_t)"]) 
plot(LATR_grow$weighted.dens,LATR_grow_fitted_terms[,"s.1(weighted.dens)"]) 

############################################################################
##2. Flowering probability
############################################################################
# populate year t of 2017-2018 transition year (there are no 2018 data but this way we get all four years in the reproduction models)
# I'll do this by creating the 2017-18 data as a stand-alone df then bind rows
LATR_dat_201718 <- LATR_full[LATR_full$year_t==2016 & LATR_full$survival_t1==1,]
## these are the 2017 survivors. Make their year t demography last year's data
LATR_dat_201718$year_t<-2017;LATR_dat_201718$year_t1<-2018
LATR_dat_201718$max.ht_t<-LATR_dat_201718$max.ht_t1
LATR_dat_201718$max.w_t<-LATR_dat_201718$max.w_t1
LATR_dat_201718$volume_t<-LATR_dat_201718$volume_t1
LATR_dat_201718$perp.w_t<-LATR_dat_201718$perp.w_t1
LATR_dat_201718$flowers_t<-LATR_dat_201718$flowers_t1
LATR_dat_201718$fruits_t<-LATR_dat_201718$fruits_t1
LATR_dat_201718$reproductive_fraction_t<-LATR_dat_201718$reproductive_fraction_t1
LATR_dat_201718$total.reproduction_t<-LATR_dat_201718$total.reproduction_t1
## now set all the t1 data to NA
LATR_dat_201718$max.ht_t1<-NA
LATR_dat_201718$max.w_t1<-NA
LATR_dat_201718$volume_t1<-NA
LATR_dat_201718$perp.w_t1<-NA
LATR_dat_201718$flowers_t1<-NA
LATR_dat_201718$fruits_t1<-NA
LATR_dat_201718$reproductive_fraction_t1<-NA
LATR_dat_201718$total.reproduction_t1<-NA
## bind rows
LATR_flow_dat <- bind_rows(LATR_full,LATR_dat_201718) %>% 
  select(unique.transect,volume_t,total.reproduction_t,weighted.dens) %>% drop_na()
## create log_vol as new variables (easier for gams)
LATR_flow_dat$log_volume_t <- log(LATR_flow_dat$volume_t)

LATR_flower <- list()
LATR_flower[[1]] <-  gam(total.reproduction_t>0 ~ s(log_volume_t) + s(unique.transect,bs="re"),
                         data=LATR_flow_dat, gamma=1.4, family="binomial")
LATR_flower[[2]] <-  gam(total.reproduction_t>0 ~ s(log_volume_t) + s(weighted.dens) + s(unique.transect,bs="re"),
                         data=LATR_flow_dat, gamma=1.4, family="binomial")
LATR_flower[[3]] <-  gam(total.reproduction_t>0 ~ s(log_volume_t) + s(weighted.dens) + weighted.dens:log_volume_t  + s(unique.transect,bs="re"),
                         data=LATR_flow_dat, gamma=1.4, family="binomial")
flower_aic<-AICtab(LATR_flower,base=T,sort=F)
LATR_flower_best <- LATR_flower[[which.min(flower_aic$AIC)]]
LATR_flower_fitted_terms = predict(LATR_flower_best,type="terms") 
LATR_flow_dat$pred = predict.gam(LATR_flower_best,newdata = LATR_flow_dat, exclude = "s(unique.transect)")

## Tom's practice bootstrap
#for(b in 1:n.boot){
#  LATR_flower_best <- gam(total.reproduction_t>0 ~ s(log_volume_t) + s(weighted.dens) + weighted.dens:log_volume_t  + s(unique.transect,bs="re"),
#                         data=LATR_flow_dat_boot[[b]], gamma=1.4, family="binomial")
#}

##### effect of size on pr(flower)
plot(LATR_flow_dat$log_volume_t,LATR_flower_fitted_terms[,"s(log_volume_t)"]) 
#### effect of density on pr(flower)
plot(LATR_flow_dat$weighted.dens,LATR_flower_fitted_terms[,"s(weighted.dens)"]) 

## visualize data + model
n_cuts_dens <- 3
n_cuts_size <- 5
LATR_flow_dat %>% 
  mutate(size_bin = as.integer(cut_number(log_volume_t,n_cuts_size)),
         dens_bin = as.integer(cut_number(weighted.dens,n_cuts_dens))) %>% 
  group_by(size_bin,dens_bin) %>% 
  summarise(mean_size = mean(log_volume_t),
            mean_density = mean(weighted.dens),
            mean_flower = mean(total.reproduction_t > 0),
            pred_flower = mean(pred),
            bin_n = n()) -> LATR_flow_dat_plot

## generate predictions for plotting
size_means_flow <- LATR_flow_dat_plot %>% group_by(size_bin) %>% summarise(mean_size=mean(mean_size))
LATR_flow_pred <- data.frame(
  weighted.dens = rep(seq(min(LATR_flow_dat$weighted.dens),max(LATR_flow_dat$weighted.dens),length.out = 20),times=n_cuts_size),
  log_volume_t = rep(size_means_flow$mean_size,each=20),
  unique.transect="1.FPS",
  size_bin = rep(size_means_flow$size_bin,each=20)
)
LATR_flow_pred$pred <- predict.gam(LATR_flower_best,newdata = LATR_flow_pred, exclude = "s(unique.transect)")

plot(LATR_flow_dat_plot$mean_density,LATR_flow_dat_plot$mean_flower,type="n",ylim=c(0,1),
     xlab="Weighted density",ylab="Pr(Flowering)")
for(i in 1:n_cuts_size){
  points(LATR_flow_dat_plot$mean_density[LATR_flow_dat_plot$size_bin==i],
         LATR_flow_dat_plot$mean_flower[LATR_flow_dat_plot$size_bin==i],pch=16,col=i,
         cex=(LATR_flow_dat_plot$bin_n[LATR_flow_dat_plot$size_bin==i]/max(LATR_flow_dat_plot$bin_n))*3)
  lines(LATR_flow_pred$weighted.dens[LATR_flow_pred$size_bin==i],
        invlogit(LATR_flow_pred$pred[LATR_flow_pred$size_bin==i]),col=i)
}


############################################################################
##3. Fruit production
############################################################################
LATR_fruits_dat <- subset(LATR_flow_dat,total.reproduction_t>0)
LATR_fruits <- list()
LATR_fruits[[1]] <-  gam(total.reproduction_t ~ s(log_volume_t) + s(unique.transect,bs="re"),
                         data=LATR_fruits_dat, gamma=1.4, family="nb")
LATR_fruits[[2]] <-  gam(total.reproduction_t ~ s(log_volume_t) + s(weighted.dens) + s(unique.transect,bs="re"),
                         data=LATR_fruits_dat, gamma=1.4, family="nb")
LATR_fruits[[3]] <-  gam(total.reproduction_t ~ s(log_volume_t) + s(weighted.dens) + weighted.dens:log_volume_t  + s(unique.transect,bs="re"),
                         data=LATR_fruits_dat, gamma=1.4, family="nb")
fruits_aic<-AICtab(LATR_fruits,base=T,sort=F)
LATR_fruits_best <- LATR_fruits[[which.min(fruits_aic$AIC)]]
LATR_fruits_fitted_terms = predict(LATR_fruits_best,type="terms") 
LATR_fruits_dat$pred = predict.gam(LATR_fruits_best,newdata = LATR_fruits_dat,exclude="s(unique.transect)")

##### effect of size on fruits
plot(LATR_fruits_dat$log_volume_t,LATR_fruits_fitted_terms[,"s(log_volume_t)"]) 
#### effect of density on fruits 
plot(LATR_fruits_dat$weighted.dens,LATR_fruits_fitted_terms[,"s(weighted.dens)"]) 

LATR_fruits_dat %>% 
  mutate(size_bin = as.integer(cut_number(log_volume_t,n_cuts_size)),
         dens_bin = as.integer(cut_number(weighted.dens,n_cuts_dens))) %>% 
  group_by(size_bin,dens_bin) %>% 
  summarise(mean_size = mean(log_volume_t),
            mean_density = mean(weighted.dens),
            mean_fruits = mean(total.reproduction_t),
            pred_fruits = mean(pred),
            bin_n = n()) -> LATR_fruits_dat_plot

## new data set for gam prediction
size_means_fruit <- LATR_fruits_dat_plot %>% group_by(size_bin) %>% summarise(mean_size=mean(mean_size))
LATR_fruit_pred <- data.frame(
  weighted.dens = rep(seq(min(LATR_fruits_dat$weighted.dens),max(LATR_fruits_dat$weighted.dens),length.out = 20),times=n_cuts_size),
  log_volume_t = rep(size_means_fruit$mean_size,each=20),
  unique.transect="1.FPS",
  size_bin = rep(size_means_fruit$size_bin,each=20)
)
LATR_fruit_pred$pred <- predict.gam(LATR_fruits_best,newdata = LATR_fruit_pred, exclude = "s(unique.transect)")


plot(LATR_fruits_dat_plot$mean_density,LATR_fruits_dat_plot$mean_fruits,type="n",
     xlab="Weighted density",ylab="Flowers and Fruits")
for(i in 1:n_cuts_size){
  points(LATR_fruits_dat_plot$mean_density[LATR_fruits_dat_plot$size_bin==i],
         LATR_fruits_dat_plot$mean_fruits[LATR_fruits_dat_plot$size_bin==i],pch=16,col=i,
         cex=(LATR_fruits_dat_plot$bin_n[LATR_fruits_dat_plot$size_bin==i]/max(LATR_fruits_dat_plot$bin_n))*3)
  lines(LATR_fruit_pred$weighted.dens[LATR_fruit_pred$size_bin==i],
        exp(LATR_fruit_pred$pred[LATR_fruit_pred$size_bin==i]),col=i)
}

############################################################################
##4. Survival
############################################################################
# Combine transplants with large shrubs
# Keep only location info, survival, volume, and density
CData.Transplants %>% 
  select("site", "transect", "actual.window", 
         "spring_survival_t1", "volume_t", "weighted.dens", "transplant") %>% 
  rename("survival_t1" = "spring_survival_t1") %>% 
  mutate(unique.transect = interaction(transect, site)) %>% 
  rbind(select(LATR_full, "site", "transect", "actual.window", 
               "survival_t1", "volume_t", "weighted.dens", "transplant","unique.transect")) %>% 
  mutate(log_volume_t = log(volume_t)) %>% 
  drop_na()-> LATR_surv_dat

## how much size overlap do we have between transplant experiment and observational census?
hist(log(LATR_surv_dat$volume_t[LATR_surv_dat$transplant==F]))
hist(log(LATR_surv_dat$volume_t[LATR_surv_dat$transplant==T]),add=T,col=alpha("gray",0.5))

plot(log(LATR_surv_dat$volume_t[LATR_surv_dat$transplant==F]),
     LATR_surv_dat$survival_t1[LATR_surv_dat$transplant==F])
points(log(LATR_surv_dat$volume_t[LATR_surv_dat$transplant==T]),
       LATR_surv_dat$survival_t1[LATR_surv_dat$transplant==T]-0.025,pch=2)

LATR_surv <- list()
LATR_surv[[1]] <-  gam(survival_t1 ~ s(log_volume_t) + transplant + s(unique.transect,bs="re"),
                       data=LATR_surv_dat, gamma=1.4, family="binomial")
LATR_surv[[2]] <-  gam(survival_t1 ~ s(log_volume_t) + s(weighted.dens)  + transplant + s(unique.transect,bs="re"),
                       data=LATR_surv_dat, gamma=1.4, family="binomial")
LATR_surv[[3]] <-  gam(survival_t1 ~ s(log_volume_t) + s(weighted.dens) + transplant + weighted.dens:log_volume_t + s(unique.transect,bs="re"),
                       data=LATR_surv_dat, gamma=1.4, family="binomial")
surv_aic<-AICtab(LATR_surv,base=T,sort=F)
LATR_surv_best <- LATR_surv[[which.min(surv_aic$AIC)]]
LATR_surv_fitted_terms = predict(LATR_surv_best,type="terms") 
LATR_surv_dat$pred = predict.gam(LATR_surv_best,newdata = LATR_surv_dat,exclude="s(unique.transect)")

##### effect of size on pr(survival)
plot(LATR_surv_dat$log_volume_t,LATR_surv_fitted_terms[,"s(log_volume_t)"]) 
#### effect of density on pr(survival)
plot(LATR_surv_dat$weighted.dens,LATR_surv_fitted_terms[,"s(weighted.dens)"]) 

## visualize data + model -- this is for the natural census
n_cuts_dens <- 4
n_cuts_size <- 4
LATR_surv_dat %>% 
  filter(transplant==F) %>% 
  mutate(size_bin = as.integer(cut_interval(log_volume_t,n_cuts_size)),
         dens_bin = as.integer(cut_interval(weighted.dens,n_cuts_dens))) %>% 
  group_by(size_bin,dens_bin) %>% 
  summarise(mean_size = mean(log_volume_t),
            mean_density = mean(weighted.dens),
            mean_surv = mean(survival_t1),
            pred_surv = mean(pred),
            bin_n = n()) -> LATR_surv_nat_plot

## generate predictions for plotting
size_means_surv_nat <- LATR_surv_nat_plot %>% group_by(size_bin) %>% summarise(mean_size=mean(mean_size))
LATR_surv_nat_pred <- data.frame(
  weighted.dens = rep(seq(min(LATR_surv_nat_plot$mean_density),max(LATR_surv_nat_plot$mean_density),length.out = 20),times=n_cuts_size),
  log_volume_t = rep(size_means_surv_nat$mean_size,each=20),
  unique.transect="1.FPS",
  transplant=F,
  size_bin = rep(size_means_surv_nat$size_bin,each=20)
)
LATR_surv_nat_pred$pred <- predict.gam(LATR_surv_best,newdata = LATR_surv_nat_pred, exclude = "s(unique.transect)")

plot(LATR_surv_nat_plot$mean_density,LATR_surv_nat_plot$mean_surv,type="n",ylim=c(0,1),
     xlab="Weighted density",ylab="Pr(Survival)")
for(i in 1:n_cuts_size){
  points(LATR_surv_nat_plot$mean_density[LATR_surv_nat_plot$size_bin==i],
         LATR_surv_nat_plot$mean_surv[LATR_surv_nat_plot$size_bin==i],pch=16,col=i,cex=2)
         #cex=(LATR_surv_nat_plot$bin_n[LATR_surv_nat_plot$size_bin==i]/max(LATR_surv_nat_plot$bin_n))*3)
  lines(LATR_surv_nat_pred$weighted.dens[LATR_surv_nat_pred$size_bin==i],
        invlogit(LATR_surv_nat_pred$pred[LATR_surv_nat_pred$size_bin==i]),col=i)
}

## now transplants
LATR_surv_dat %>% 
  filter(transplant==T) %>% 
  mutate(dens_bin = as.integer(cut_interval(weighted.dens,n_cuts_dens)),
         mean_size = mean(log_volume_t)) %>% 
  group_by(dens_bin) %>% 
  summarise(mean_size = unique(mean_size),
            mean_density = mean(weighted.dens),
            mean_surv = mean(survival_t1),
            bin_n = n()) -> LATR_surv_exp_plot
LATR_surv_exp_pred <- data.frame(
  weighted.dens = seq(min(LATR_surv_exp_plot$mean_density),max(LATR_surv_exp_plot$mean_density),length.out = 20),
  log_volume_t = LATR_surv_exp_plot$mean_size[1],
  unique.transect="1.FPS",
  transplant=T
)
LATR_surv_exp_pred$pred <- predict.gam(LATR_surv_best,newdata = LATR_surv_exp_pred, exclude = "s(unique.transect)")

points(LATR_surv_exp_plot$mean_density,LATR_surv_exp_plot$mean_surv,ylim=c(0,1),pch=2)
lines(LATR_surv_exp_pred$weighted.dens,invlogit(LATR_surv_exp_pred$pred),lty=2)

############################################################################
##5. Per-seed recruitment probability
############################################################################

LATR_recruits <- LATR_full %>% 
  mutate(unique.transect = interaction(transect, site)) %>% 
  group_by(year_t1,unique.transect,actual.window) %>% 
  filter(seedling_t1==1) %>% 
  summarise(recruits = n()) %>% 
  rename(window=actual.window)

## now estimate total seeds produced in each window using the known plant sizes and the fitted flowering and fruiting models
LATR_transects <- Cdata.Transects.Windows %>% 
  mutate(unique.transect = interaction(transect, site),
         log_volume_t = log(volume))
LATR_transects$seeds = ceiling(invlogit(predict.gam(LATR_flower_best,newdata = LATR_transects)) * 
                                 6*exp(predict.gam(LATR_fruits_best,newdata = LATR_transects))) ## note hard-coded # seeds per fruit
LATR_transects %>% 
  group_by(unique.transect,window) %>% 
  summarise(total_seeds=sum(seeds),
            weighted.dens = unique(weighted.dens)) -> LATR_transects

## now do something weird. take three copies of this df, assigning each one to a different year and assigning recruits to zero (for now)
LATR_recruitment <- bind_rows(LATR_transects %>% filter(unique.transect=="1.FPS"|unique.transect=="2.FPS"|unique.transect=="3.FPS") %>% 
                                mutate(year_t1=2014,recruits=0), ## only FPS for 2013-2014
                              LATR_transects %>% mutate(year_t1=2015,recruits=0),
                              LATR_transects %>% mutate(year_t1=2016,recruits=0),
                              LATR_transects %>% mutate(year_t1=2017,recruits=0)) %>% 
  left_join(.,LATR_recruits,by=c("year_t1","unique.transect","window")) %>% 
  mutate(recruits.y=replace_na(recruits.y,0),
         recruits = pmax(recruits.x,recruits.y,na.rm=T)) %>% 
  drop_na()

LATR_recruit <- list()
LATR_recruit[[1]] <-  gam(cbind(recruits,total_seeds-recruits) ~ s(unique.transect,bs="re"),
                          data=LATR_recruitment, gamma=1.4, family="binomial")
LATR_recruit[[2]] <-  gam(cbind(recruits,total_seeds-recruits) ~ s(weighted.dens) + s(unique.transect,bs="re"),
                          data=LATR_recruitment, gamma=1.4, family="binomial")
recruit_aic<-AICtab(LATR_recruit,base=T,sort=F)
LATR_recruit_best <- LATR_recruit[[which.min(recruit_aic$AIC)]]

plot(LATR_recruitment$weighted.dens,LATR_recruitment$recruits/LATR_recruitment$total_seeds)
LATR_recruitment$pred = predict.gam(LATR_recruit_best,newdata = LATR_recruitment,exclude="s(unique.transect)")
points(LATR_recruitment$weighted.dens,invlogit(LATR_recruitment$pred),col="red",pch=".")
## no evidence for density dependence in recruitment, just a really low overall recruitment rate

## just out of curiosity, the density-dependent model is a very close second...what does this look like?
LATR_recruit_fitted_terms = predict(LATR_recruit[[2]],type="terms") 
#### effect of density on pr(seedling recruitment) -- negative DD
plot(LATR_recruitment$weighted.dens,LATR_recruit_fitted_terms[,"s(weighted.dens)"]) 

############################################################################
##6. Seedling size distribution 
############################################################################
LATR_recruit_size <- LATR_full %>% 
  filter(seedling_t1==1) %>% 
  mutate(log_volume = log(volume_t1))

hist(LATR_recruit_size$log_volume)
LATR_recruit_size <- data.frame(recruit_mean = mean(LATR_recruit_size$log_volume),
                   recruit_sd = sd(LATR_recruit_size$log_volume))

############################################################################
##7. Integration limits (size bounds)
############################################################################
LATR_size_bounds <- data.frame(min_size = log(min(LATR_full$volume_t,LATR_full$volume_t1[LATR_full$transplant==F],na.rm=T)),
                   max_size = log(max(LATR_full$volume_t,LATR_full$volume_t1[LATR_full$transplant==F],na.rm=T)))

############################################################################
##8. IPM functions
############################################################################

## Growth -- Gaussian using best gam
growth_fn <- function(x,y,d){
  xb=pmin(pmax(x,LATR_size_bounds$min_size),LATR_size_bounds$max_size)
  lpmat <- predict.gam(LATR_grow_best,
                       newdata = data.frame(
                         weighted.dens = d,
                         log_volume_t = xb,
                         unique.transect="1.FPS"),
                       type="lpmatrix",
                       exclude = "s(unique.transect)")
  ## linear predictor for mean and log sigma -- need to update so these indices are not hard-coded but for now they work
  grow_mu <- lpmat[,1:19]%*%coef(LATR_grow_best)[1:19]
  grow_sigma <- exp(lpmat[,32:50]%*%coef(LATR_grow_best)[32:50])
  return(dnorm(y,mean=grow_mu,sd=grow_sigma))
}

## Survival-- prediction from naturally occuring plants (transplant=F)
survival_fn <- function(x,d){
  xb=pmin(pmax(x,LATR_size_bounds$min_size),LATR_size_bounds$max_size)
  lpmat <- predict.gam(LATR_surv_best,
                       newdata = data.frame(
                         weighted.dens = d,
                         log_volume_t = xb,
                         transplant=F,
                         unique.transect="1.FPS"),
                       type="lpmatrix",
                       exclude = "s(unique.transect)")
  pred <- lpmat[,1:21]%*%coef(LATR_surv_best)[1:21]
  return(invlogit(pred))
}

## combined growth and survival
pxy <- function(x,y,d){
  survival_fn(x,d) * growth_fn(x,y,d)
}

## Flowering
flower_fn <- function(x,d){
  xb=pmin(pmax(x,LATR_size_bounds$min_size),LATR_size_bounds$max_size)
  lpmat <- predict.gam(LATR_flower_best,
                       newdata = data.frame(
                         weighted.dens = d,
                         log_volume_t = xb,
                         unique.transect="1.FPS"),
                       type="lpmatrix",
                       exclude = "s(unique.transect)")
  pred <- lpmat[,1:20]%*%coef(LATR_flower_best)[1:20]
  return(invlogit(pred))
}

## seed production (fruits * seeds/fruit)
seeds_fn <- function(x,d,seeds.per.fruit=6){
  xb=pmin(pmax(x,LATR_size_bounds$min_size),LATR_size_bounds$max_size)
  lpmat <- predict.gam(LATR_fruits_best,
                       newdata = data.frame(
                         weighted.dens = d,
                         log_volume_t = xb,
                         unique.transect="1.FPS"),
                       type="lpmatrix",
                       exclude = "s(unique.transect)")
  pred <- lpmat[,1:19]%*%coef(LATR_fruits_best)[1:19]
  return(exp(pred)*seeds.per.fruit)
}

## Seed-to-Seedling recruitment probability
recruitment_fn <- function(d){
  lpmat <- predict.gam(LATR_recruit_best,
                       newdata = data.frame(
                         weighted.dens = d,
                         unique.transect="1.FPS"),
                       type="lpmatrix",
                       exclude = "s(unique.transect)")
  pred <- lpmat[,1]%*%coef(LATR_recruit_best)[1]
  return(invlogit(pred[[1]]))
}

## Recruit size distribution
recruit_size<-function(y){
  dnorm(x=y,mean=LATR_recruit_size$recruit_mean,sd=LATR_recruit_size$recruit_sd)
}

## combined flowering, fertility, and recruitment
fxy <- function(x,y,d){
  flower_fn(x,d) * seeds_fn(x,d) * recruitment_fn(d) * recruit_size(y)
}

#PUT IT ALL TOGETHER
## projection matrix is a function of weighted density (dens)
bigmatrix<-function(lower.extension = -8, ## needs a large lower extension because growth variance (gaussian) is greater for smaller plants 
                    upper.extension = 2,
                    min.size = LATR_size_bounds$min_size,
                    max.size = LATR_size_bounds$max_size,
                    mat.size = 200,
                    dens){
  
  n<-mat.size
  L<-min.size + lower.extension
  U<-max.size + upper.extension
  #these are the upper and lower integration limits
  h<-(U-L)/n                   #Bin size
  b<-L+c(0:n)*h;               #Lower boundaries of bins 
  y<-0.5*(b[1:n]+b[2:(n+1)]);  #Bin midpoints
  
  # Growth/Survival matrix
  Pmat<-t(outer(y,y,pxy,d=dens)) * h 
  
  # Fertility/Recruiment matrix
  Fmat<-t(outer(y,y,fxy,d=dens)) * h 
  
  # Put it all together
  IPMmat<-Pmat+Fmat
  
  #and transition matrix
  return(list(IPMmat=IPMmat,Fmat=Fmat,Pmat=Pmat,meshpts=y))
}

############################################################################
##9. IPM analysis
############################################################################
## lambda over density variation
density_dummy <- seq(min(LATR_full$weighted.dens,na.rm=T),max(LATR_full$weighted.dens,na.rm=T),length.out = 10)
lambda_density <- c()
for(d in 1:length(density_dummy)){
  print(d)
  lambda_density[d] <- lambda(bigmatrix(dens=density_dummy[d], mat.size=200)$IPMmat)
}

plot(density_dummy,lambda_density,type="l",lwd=3,xlab="Weighted density",ylab="lambda")
abline(h=1,lty=3)


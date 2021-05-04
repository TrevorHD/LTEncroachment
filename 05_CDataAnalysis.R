##### Prepare data for analysis ---------------------------------------------------------------------------

# Script authored by Tom, with some changes and additions from Trevor

# Define inverse logit and prepare data for model fitting
if(boot.switch == FALSE){
  
  # Define inverse logit function for later use; need to do this only once
  invlogit <- function(x){exp(x)/(1 + exp(x))}
  
  # Create unique transect as interaction of transect and site
  # First-time only; this bit of code in 06_BootRes does it all other times
  LATR_full <- CData %>% 
    mutate(unique.transect = interaction(transect, site))}


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

# Pilot fits, where sigma depends on initial size only
# constant sigma
LATR_gam_models[[1]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(unique.transect, bs = "re"), ~1),
                            data = LATR_grow, gamma = gamma, family = gaulss())
LATR_gam_models[[2]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(weighted.dens) + s(unique.transect, bs = "re"), ~1), 
                            data = LATR_grow, gamma = gamma, family = gaulss())                
LATR_gam_models[[3]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(weighted.dens) + ti(log_volume_t,weighted.dens) + s(unique.transect, bs = "re"), ~1), 
                            data = LATR_grow, gamma = gamma, family = gaulss()) 

LATR_gam_models[[4]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(unique.transect,bs = "re"), ~s(log_volume_t)), 
                            data = LATR_grow, gamma = gamma, family = gaulss())
LATR_gam_models[[5]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(weighted.dens) + s(unique.transect, bs = "re"), ~s(log_volume_t)), 
                            data = LATR_grow, gamma = gamma, family = gaulss())                
LATR_gam_models[[6]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(weighted.dens) + ti(log_volume_t,weighted.dens) + s(unique.transect,bs = "re"), ~s(log_volume_t)), 
                            data = LATR_grow, gamma = gamma, family = gaulss()) 

# Fits where sigma depends on both initial size and density
LATR_gam_models[[7]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(unique.transect, bs = "re"), ~s(log_volume_t) + s(weighted.dens)), 
                            data = LATR_grow, gamma = gamma, family = gaulss())
LATR_gam_models[[8]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(weighted.dens) + s(unique.transect, bs = "re"), ~s(log_volume_t) + s(weighted.dens)), 
                            data = LATR_grow, gamma = gamma, family = gaulss())                
LATR_gam_models[[9]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(weighted.dens) + ti(log_volume_t,weighted.dens) + s(unique.transect, bs = "re"), ~s(log_volume_t) + s(weighted.dens)), 
                            data = LATR_grow, gamma = gamma, family = gaulss()) 



# Collect model AICs into a single table
grow_aic <- AICtab(LATR_gam_models, base = TRUE, sort = FALSE)

# Set top model as "best"; find intercept for the sd
LATR_grow_best <- LATR_gam_models[[which.min(grow_aic$AIC)]]
grow_sd_index <- which(as.factor(names(coef(LATR_grow_best)))=="(Intercept).1")
LATR_grow_fitted_terms <- predict(LATR_grow_best, type = "terms") 
LATR_grow$pred <- predict.gam(LATR_grow_best, newdata = LATR_grow, exclude = "s(unique.transect)")

# Plot of effect of size on future size -- obviously linear
# plot(LATR_grow$log_volume_t, LATR_grow_fitted_terms[, "s(log_volume_t)"]) 

# Plot of effect of density on growth 
# plot(LATR_grow$weighted.dens, LATR_grow_fitted_terms[, "s(weighted.dens)"]) 

# Plots of effect of size and density on sd(future size)
# plot(LATR_grow$log_volume_t, LATR_grow_fitted_terms[, "s.1(log_volume_t)"]) 
# plot(LATR_grow$weighted.dens, LATR_grow_fitted_terms[, "s.1(weighted.dens)"]) 

##### Flowering probability model -------------------------------------------------------------------------

# Populate year t of 2017-2018 transition year
# There are no 2018 data but this way we get all four years in the reproduction models
# Do this by creating the 2017-18 data as a stand-alone df then bind rows
LATR_dat_201718 <- LATR_full[LATR_full$year_t == 2016 & LATR_full$survival_t1 == 1, ]

# These are the 2017 survivors; make their year t demography last year's data
LATR_dat_201718$year_t <- 2017
LATR_dat_201718$year_t1 <- 2018
LATR_dat_201718$max.ht_t <- LATR_dat_201718$max.ht_t1
LATR_dat_201718$max.w_t <- LATR_dat_201718$max.w_t1
LATR_dat_201718$volume_t <- LATR_dat_201718$volume_t1
LATR_dat_201718$perp.w_t <- LATR_dat_201718$perp.w_t1
LATR_dat_201718$flowers_t <- LATR_dat_201718$flowers_t1
LATR_dat_201718$fruits_t <- LATR_dat_201718$fruits_t1
LATR_dat_201718$reproductive_fraction_t <- LATR_dat_201718$reproductive_fraction_t1
LATR_dat_201718$total.reproduction_t <- LATR_dat_201718$total.reproduction_t1

# Now set all the t1 data to NA
LATR_dat_201718$max.ht_t1 <- NA
LATR_dat_201718$max.w_t1 <- NA
LATR_dat_201718$volume_t1 <- NA
LATR_dat_201718$perp.w_t1 <- NA
LATR_dat_201718$flowers_t1 <- NA
LATR_dat_201718$fruits_t1 <- NA
LATR_dat_201718$reproductive_fraction_t1 <- NA
LATR_dat_201718$total.reproduction_t1 <- NA

# Bind rows and create log_vol as new variables (easier for GAMs)
LATR_flow_dat <- bind_rows(LATR_full,LATR_dat_201718) %>% 
  dplyr::select(unique.transect,volume_t,total.reproduction_t,weighted.dens) %>% drop_na()
LATR_flow_dat$log_volume_t <- log(LATR_flow_dat$volume_t)

# Create empty list to populate with model results
LATR_flower <- list()

# Three candidate models for the mean: size only, size + density, or size, density, and size:density
LATR_flower[[1]] <- gam(total.reproduction_t > 0 ~ s(log_volume_t) + s(unique.transect, bs = "re"),
                        data = LATR_flow_dat, gamma = gamma, family = "binomial")
LATR_flower[[2]] <- gam(total.reproduction_t > 0 ~ s(log_volume_t) + s(weighted.dens) + s(unique.transect, bs = "re"),
                        data = LATR_flow_dat, gamma = gamma, family = "binomial")
LATR_flower[[3]] <- gam(total.reproduction_t > 0 ~ s(log_volume_t) + s(weighted.dens) + ti(log_volume_t,weighted.dens) + s(unique.transect, bs = "re"),
                        data = LATR_flow_dat, gamma = gamma, family = "binomial")

# Collect model AICs into a single table
flower_aic<-AICtab(LATR_flower, base = TRUE, sort = FALSE)

# Set top model as "best"
LATR_flower_best <- LATR_flower[[which.min(flower_aic$AIC)]]
LATR_flower_fitted_terms <- predict(LATR_flower_best, type = "terms") 
LATR_flow_dat$pred <- predict.gam(LATR_flower_best, newdata = LATR_flow_dat, exclude = "s(unique.transect)")

##### Fruit production model ------------------------------------------------------------------------------

# Create new df with plants that have produced at least one reproductive structure
LATR_fruits_dat <- subset(LATR_flow_dat, total.reproduction_t > 0)

# Create empty list to populate with model results
LATR_fruits <- list()

# Three candidate models for the mean: size only, size + density, or size, density, and size:density
LATR_fruits[[1]] <- gam(total.reproduction_t ~ s(log_volume_t) + s(unique.transect, bs = "re"),
                        data = LATR_fruits_dat, gamma = gamma, family = "nb")
LATR_fruits[[2]] <- gam(total.reproduction_t ~ s(log_volume_t) + s(weighted.dens) + s(unique.transect, bs = "re"),
                        data = LATR_fruits_dat, gamma = gamma, family = "nb")
LATR_fruits[[3]] <- gam(total.reproduction_t ~ s(log_volume_t) + s(weighted.dens) + weighted.dens:log_volume_t + s(unique.transect, bs = "re"),
                        data = LATR_fruits_dat, gamma = gamma, family = "nb")
# Collect model AICs into a single table
fruits_aic <- AICtab(LATR_fruits, base = TRUE, sort = FALSE)

# Set top model as "best"
LATR_fruits_best <- LATR_fruits[[which.min(fruits_aic$AIC)]]
LATR_fruits_fitted_terms <- predict(LATR_fruits_best, type = "terms") 
LATR_fruits_dat$pred <- predict.gam(LATR_fruits_best, newdata = LATR_fruits_dat, exclude = "s(unique.transect)")

# Plot effect of size on fruits
# plot(LATR_fruits_dat$log_volume_t, LATR_fruits_fitted_terms[, "s(log_volume_t)"]) 

# Plot effect of density on fruits 
# plot(LATR_fruits_dat$weighted.dens, LATR_fruits_fitted_terms[, "s(weighted.dens)"]) 





##### Survival model --------------------------------------------------------------------------------------

# Combine transplants with large shrubs; keep only location info, survival, volume, and density
CData.Transplants %>% 
  dplyr::select("site", "transect", "actual.window", 
         "spring_survival_t1", "volume_t", "weighted.dens", "transplant") %>% 
  rename("survival_t1" = "spring_survival_t1") %>% 
  mutate(unique.transect = interaction(transect, site)) %>% 
  rbind(dplyr::select(LATR_full, "site", "transect", "actual.window", 
               "survival_t1", "volume_t", "weighted.dens", "transplant","unique.transect")) %>% 
  mutate(log_volume_t = log(volume_t)) %>% 
  drop_na() -> LATR_surv_dat

# Investigate size overlap between transplant experiment and observational census
hist(log(LATR_surv_dat$volume_t[LATR_surv_dat$transplant == FALSE]), breaks = 25)
hist(log(LATR_surv_dat$volume_t[LATR_surv_dat$transplant == TRUE]), breaks = 10, add = TRUE, col = alpha("gray", 0.5))

# Plot survival against volume, grouped by transplant status
plot(log(LATR_surv_dat$volume_t[LATR_surv_dat$transplant == FALSE]),
     LATR_surv_dat$survival_t1[LATR_surv_dat$transplant == FALSE])
points(log(LATR_surv_dat$volume_t[LATR_surv_dat$transplant == TRUE]),
       LATR_surv_dat$survival_t1[LATR_surv_dat$transplant == TRUE] - 0.025, pch = 2)

# Create empty list to populate with model results
LATR_surv <- list()

# Three candidate models for the mean: size only, size + density, or size, density, and size:density
LATR_surv[[1]] <- gam(survival_t1 ~ s(log_volume_t) + transplant + s(unique.transect, bs = "re"),
                      data = LATR_surv_dat, gamma = gamma, family = "binomial")
LATR_surv[[2]] <- gam(survival_t1 ~ s(log_volume_t) + s(weighted.dens)  + transplant + s(unique.transect, bs = "re"),
                      data = LATR_surv_dat, gamma = gamma, family = "binomial")
LATR_surv[[3]] <- gam(survival_t1 ~ s(log_volume_t) + s(weighted.dens) + transplant + weighted.dens:log_volume_t + s(unique.transect, bs = "re"),
                      data = LATR_surv_dat, gamma = gamma, family = "binomial")

# Collect model AICs into a single table
surv_aic <- AICtab(LATR_surv, base = TRUE, sort = FALSE)

# Set top model as "best"
LATR_surv_best <- LATR_surv[[which.min(surv_aic$AIC)]]
LATR_surv_fitted_terms <- predict(LATR_surv_best, type = "terms") 
LATR_surv_dat$pred <- predict.gam(LATR_surv_best, newdata = LATR_surv_dat, exclude = "s(unique.transect)")

# Plot effect of size on pr(survival)
# plot(LATR_surv_dat$log_volume_t, LATR_surv_fitted_terms[, "s(log_volume_t)"]) 

# Plot effect of density on pr(survival)
# plot(LATR_surv_dat$weighted.dens, LATR_surv_fitted_terms[, "s(weighted.dens)"]) 





##### Per-seed recruitment probability model --------------------------------------------------------------

## number of seeds per fruit
seeds_per_fruit <- 5

# Create subset df of recruits
LATR_recruits <- LATR_full %>% 
  mutate(unique.transect = interaction(transect, site)) %>% 
  group_by(year_t1, unique.transect, actual.window) %>% 
  filter(seedling_t1 == 1)
  suppressMessages(summarise(LATR_recruits, recruits = n())) %>% 
  rename(window = actual.window) -> LATR_recruits

# Estimate total seeds produced in each window
# This is computed using the known plant sizes and the fitted flowering and fruiting models
# Note: we assume 6 seeds per fruit -- updated 18 Apr 2021 to 5
LATR_transects <- Cdata.Transects.Windows %>% 
  mutate(unique.transect = interaction(transect, site),
         log_volume_t = log(volume))
LATR_transects$seeds = ceiling(invlogit(predict.gam(LATR_flower_best, newdata = LATR_transects))* 
                                 seeds_per_fruit*exp(predict.gam(LATR_fruits_best, newdata = LATR_transects)))
LATR_transects <- group_by(LATR_transects, unique.transect, window)
suppressMessages(summarise(LATR_transects, total_seeds = sum(seeds),
                                           weighted.dens = unique(weighted.dens))) -> LATR_transects

# Take three copies of this df, assigning each one to a different year and assigning recruits to zero (for now)
LATR_recruitment <- bind_rows(LATR_transects %>% filter(unique.transect == "1.FPS" | unique.transect == "2.FPS" | unique.transect == "3.FPS") %>% 
                                mutate(year_t1 = 2014, recruits = 0), ## only FPS for 2013-2014
                              LATR_transects %>% mutate(year_t1 = 2015, recruits = 0),
                              LATR_transects %>% mutate(year_t1 = 2016, recruits = 0),
                              LATR_transects %>% mutate(year_t1 = 2017, recruits = 0)) %>% 
  left_join(., LATR_recruits, by = c("year_t1", "unique.transect", "window")) %>% 
  mutate(recruits.y = replace_na(recruits.y, 0),
         recruits = pmax(recruits.x, recruits.y, na.rm = T)) %>% 
  drop_na()

# Create empty list to populate with model results
LATR_recruit <- list()

# Two candidate models for the mean: no effect, or weighted density only
LATR_recruit[[1]] <- gam(cbind(recruits,total_seeds - recruits) ~ s(unique.transect, bs = "re"),
                         data = LATR_recruitment, gamma = gamma, family = "binomial")
LATR_recruit[[2]] <- gam(cbind(recruits,total_seeds - recruits) ~ s(weighted.dens) + s(unique.transect, bs = "re"),
                         data = LATR_recruitment, gamma = gamma, family = "binomial")

# Collect model AICs into a single table
recruit_aic <- AICtab(LATR_recruit, base = TRUE, sort = FALSE)

# Set top model as "best"
LATR_recruit_best <- LATR_recruit[[which.min(recruit_aic$AIC)]]

# Plot effect of density on pr(seedling recruitment); negative density dependence
# LATR_recruit_fitted_terms <- predict(LATR_recruit[[2]], type = "terms") 
# plot(LATR_recruitment$weighted.dens,LATR_recruit_fitted_terms[, "s(weighted.dens)"])





##### Recruit sizes and integration limits (size bounds) --------------------------------------------------

# Filter out seedlings from full data set; only do this once
if(boot.switch == FALSE){

  # Filter out seedlings and get their sizes
  LATR_recruit_size <- LATR_full %>% 
    filter(seedling_t1 == 1) %>% 
    mutate(log_volume = log(volume_t1))}

# Plot distribution of recruit sizes using hist(LATR_recruit_size$log_volume)

# Create maximum and minimum size bounds for the IPM
LATR_size_bounds <- data.frame(min_size = log(min(LATR_full$volume_t, LATR_full$volume_t1[LATR_full$transplant == FALSE], na.rm = TRUE)),
                               max_size = log(max(LATR_full$volume_t, LATR_full$volume_t1[LATR_full$transplant == FALSE], na.rm = TRUE)))

# Collect AIC table info --------------------------------------------------

grow_aic
flower_aic
fruits_aic
surv_aic
recruit_aic


# create vital rate figures ------------------------------------------------
size_breaks <- 4
density_breaks <- 5
## I get nicer colors if I double the breaks and take every other color
LATR_cols <- wes_palette("Zissou1", size_breaks, type = "continuous")
#LATR_cols <- rainbow(size_breaks)

# growth figure -----------------------------------------------------------
## starting with raw data visualization, binning sizes and densities
## this plot freaked me out because it clearly shows positive density
## dependence at the smallest size -- no matter how you slice up sizes
LATR_grow %>% 
  mutate(size_bin = as.numeric(cut(log_volume_t,breaks = size_breaks)),
         density_bin = as.numeric(cut(weighted.dens,breaks = density_breaks))) %>% 
  group_by(size_bin,density_bin) %>% 
  summarise(mean_size = mean(log_volume_t),
            mean_dens = mean(weighted.dens),
            mean_sizet1 = mean(log_volume_t1),
            sd_sizet1 = sd(log_volume_t1),
            bin_n = n()) -> LATR_grow_plot
LATR_grow_plot$size_col <- LATR_cols[LATR_grow_plot$size_bin]

plot(LATR_grow_plot$mean_dens,LATR_grow_plot$mean_sizet1,col=LATR_grow_plot$size_bin,pch=16,cex=2)
## however, I realized that the bins actually differ in initial size, and
## this can explain the apparent positive DD, here plotting the change in size to control for size differences
plot(LATR_grow_plot$mean_dens,(LATR_grow_plot$mean_sizet1-LATR_grow_plot$mean_size),col=LATR_grow_plot$size_bin,pch=16,cex=2)

## let's avoid binning the growth data for this read
## create dummy data frame for gam prediction
mean_size <- LATR_grow %>% 
  mutate(size_bin = as.numeric(cut(log_volume_t,breaks = size_breaks))) %>% 
  group_by(size_bin) %>% 
  summarise(mean_size=mean(log_volume_t))
grow_predict <- expand.grid(
  weighted.dens = seq(min(LATR_grow$weighted.dens),max(LATR_grow$weighted.dens),1),
  size_bin = mean_size$size_bin,
  unique.transect = LATR_grow$unique.transect[1]
)
grow_predict$log_volume_t<-mean_size$mean_size[grow_predict$size_bin]
grow_predict[,c("pred_mu","pred_sigma")] <- predict.gam(LATR_grow_best, newdata = grow_predict, type = "response", exclude = "s(unique.transect)")
LATR_grow$size_col <- LATR_cols[as.numeric(cut(LATR_grow$log_volume_t,breaks = size_breaks))]
grow_predict$size_col <- LATR_cols[grow_predict$size_bin]

## growth mean and SD inset
plot(LATR_grow$weighted.dens,LATR_grow$log_volume_t1,pch=16,cex=0.5,
     col=alpha(LATR_grow$size_col,0.75),xlab="Weighted density",ylab=expression(paste(Size[t+1]," (log ",cm^3,")")))
points(grow_predict$weighted.dens,grow_predict$pred_mu,col=grow_predict$size_col,pch=16,cex=0.25)
rect(100, -2.5, 220, 6, col="white")
plotInset(100, -1, 220, 6,
          expr=plot(grow_predict$weighted.dens,1/(grow_predict$pred_sigma^2),
                    col=grow_predict$size_col,pch=16,cex=0.25,xlab="",ylab=expression(paste("SD(",Size[t+1],")")),
                    cex.axis=0.5,mgp=c(3/2, 1/2, 0)),
          mar=c(0, 3, 0, 0))


# flowering figure --------------------------------------------------------
LATR_flow_dat %>% 
  mutate(size_bin = as.numeric(cut(log_volume_t,breaks = size_breaks)),
         density_bin = as.numeric(cut(weighted.dens,breaks = density_breaks))) %>% 
  group_by(size_bin,density_bin) %>% 
  summarise(mean_size = mean(log_volume_t),
            mean_dens = mean(weighted.dens),
            mean_flow = mean(total.reproduction_t > 0),
            bin_n = n()) -> LATR_flow_plot
LATR_flow_plot$size_col <- LATR_cols[LATR_flow_plot$size_bin]

mean_size <- LATR_flow_plot %>% group_by(size_bin) %>% summarise(mean_size=mean(mean_size))
flow_predict <- expand.grid(
  weighted.dens = seq(min(LATR_flow_dat$weighted.dens),max(LATR_flow_dat$weighted.dens),50),
  log_volume_t = mean_size$mean_size,
  unique.transect = LATR_flow_dat$unique.transect[1]
)
## associate sizes with bin numbers for plotting
flow_predict$pred <- predict.gam(LATR_flower_best, newdata = flow_predict, type = "response", exclude = "s(unique.transect)")
flow_predict$size_col <- LATR_cols[flow_predict$size_bin]

plot(LATR_flow_plot$mean_dens,LATR_flow_plot$mean_flow,col=LATR_flow_plot$size_col,pch=16,
     ylim=c(0,1),xlab="Weighted density",ylab="Flowering probability",
     cex=LATR_flow_plot$bin_n/max(LATR_flow_plot$bin_n)*2+1,
     xlim=c(min(LATR_flow_dat$weighted.dens),max(LATR_flow_dat$weighted.dens)))
points(flow_predict$weighted.dens,flow_predict$pred,col=flow_predict$size_col,pch=16,cex=0.25)

plot(flow_predict$weighted.dens,flow_predict$pred,col=flow_predict$size_col,pch=16,cex=0.25,type="p")



## This is Tom's file to figure out how to automate the SIPM so that we can 
## fit gams to each boostrap replicate and use the best fit, which may differ 
## across boot reps
## Uses growth as an example (it's also the most complicated bc sigma is a function of covariates)

## Load data prep
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/04_CDataPrep.R")

## Take a data subset -- 75% subsample
LATR_full <- CData %>% 
  mutate(unique.transect = interaction(transect, site))
LATR_full <- LATR_full[sample(1:nrow(LATR_full), round(nrow(LATR_full)*0.75), replace = FALSE), ]

# Prepare a data subset for growth that drops rows missing either t or t1 size data
LATR_grow <- LATR_full  %>% drop_na(volume_t, volume_t1) %>%
  mutate(log_volume_t = log(volume_t),
         log_volume_t1 = log(volume_t1))

# Create empty list to populate with model results
LATR_gam_models <- list()
  
  # Three candidate models for the mean: size only, size + density, or size, density, and size:density
  # Two candidates for variance: size only, size + density (not including fitted value in this test example)
  
  # Fits, where sigma depends on initial size only
  LATR_gam_models[[1]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(unique.transect,bs = "re"), ~s(log_volume_t)), 
                              data = LATR_grow, gamma = 1.4, family = gaulss())
  LATR_gam_models[[2]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(weighted.dens) + s(unique.transect, bs = "re"), ~s(log_volume_t)), 
                              data = LATR_grow, gamma = 1.4, family = gaulss())                
  LATR_gam_models[[3]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(weighted.dens) + weighted.dens:log_volume_t + s(unique.transect,bs = "re"), ~s(log_volume_t)), 
                              data = LATR_grow, gamma = 1.4, family = gaulss()) 
  # Fits where sigma depends on both initial size and density
  LATR_gam_models[[4]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(unique.transect, bs = "re"), ~s(log_volume_t) + s(weighted.dens)), 
                              data = LATR_grow, gamma = 1.4, family = gaulss())
  LATR_gam_models[[5]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(weighted.dens) + s(unique.transect, bs = "re"), ~s(log_volume_t) + s(weighted.dens)), 
                              data = LATR_grow, gamma = 1.4, family = gaulss())                
  LATR_gam_models[[6]] <- gam(list(log_volume_t1 ~s(log_volume_t) + s(weighted.dens) + weighted.dens:log_volume_t + s(unique.transect, bs = "re"), ~s(log_volume_t) + s(weighted.dens)), 
                              data = LATR_grow, gamma = 1.4, family = gaulss()) 

# Model selection
grow_aic <- AICtab(LATR_gam_models, base = TRUE, sort = FALSE)
## set top model as "best"
LATR_grow_best <- LATR_gam_models[[which.min(grow_aic$AIC)]]
#LATR_grow_fitted_terms <- predict(LATR_grow_best, type = "terms") 
#LATR_grow$pred <- predict.gam(LATR_grow_best, newdata = LATR_grow, exclude = "s(unique.transect)")

## set a size and density for practice
test_df <- LATR_grow[runif(1,1,nrow(LATR_grow)),c("weighted.dens","log_volume_t","unique.transect")]

## linear predictor matrix
lpmat_mean <- predict.gam(LATR_grow_best,newdata = test_df,type = "lpmatrix",exclude = "s(unique.transect)")
## names of coefficients
gam_names <- as.factor(names(coef(LATR_grow_best)))
## there will always be an intercept for the mean and I can always find its index
which(gam_names=="(Intercept)")
## there will always be an intercept for the sd and I can always find its index
which(gam_names=="(Intercept).1")


## SIPM growth function based on best gam
# Growth from size x to y at density d, using best GAM
TM.growth <- function(x, y, d){
  xb = pmin(pmax(x, LATR_size_bounds$min_size), LATR_size_bounds$max_size)
  lpmat <- predict.gam(LATR_grow_best,
                       newdata = data.frame(weighted.dens = d, log_volume_t = xb, unique.transect = "1.FPS"),
                       type = "lpmatrix",
                       exclude = "s(unique.transect)")
  # Linear predictor for mean and log sigma 
  # Need to update so these indices are not hard-coded but for now they work
  grow_mu <- lpmat[, 1:19] %*% coef(LATR_grow_best)[1:19]
  grow_sigma <- exp(lpmat[, 32:50] %*% coef(LATR_grow_best)[32:50])
  return(dnorm(y, mean = grow_mu, sd = grow_sigma))
  }

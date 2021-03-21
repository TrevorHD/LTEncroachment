## This is Tom's file to figure out how to automate the SIPM so that we can 
## fit gams to each boostrap replicate and use the best fit, which may differ 
## across boot reps
## Uses growth as an example (it's also the most complicated bc sigma is a function of covariates)

## 

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
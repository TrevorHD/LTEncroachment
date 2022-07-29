## This is Tom's analysis of the transplant data

## create a new data frame from the existing transplant data
## for analysis that will appear in an appendix
transplants.app <- CData.Transplants %>%
  # # Add additional columns to data, starting with standardised weighted density
  mutate("d.stand" = (weighted.dens - mean(weighted.dens, na.rm = TRUE)) / sd(weighted.dens, na.rm = TRUE),
         # Total number of patches that contain any type of grass
         "total.grass" = num_black_gramma_t + num_blue_gramma_t + num_other_grass_t,
         # Total number of patches that contain any type of plant (including shrubs)
         "total.plant" = num_shrub_t + total.grass + num_other_t,
         # Unique transect identifier
         "unique.transect" = interaction(transect, site),
         "volume_t1" = log(vol(h = max.ht_t1, w = max.w_t1, p = perp.w_t1)),
         "logGR" = volume_t1 - volume_t)

## Density here is measured at two scales with two different metrics
## See how tightly those scales are correlated by summarising over plots
## First, some bookkeeping. There should be four reps of each plot number on each transect
table(transplants.app$unique.transect,transplants.app$plot)

transplant.plots <- transplants.app %>% 
  group_by(unique.transect,plot) %>% 
  summarise(shrub_cover = sum(num_shrub_t)/36,## division by 36 bc this is the total number of squares in a plot
            grass_cover = sum(total.grass)/36,
            weighted.dens = unique(weighted.dens),
            fall_survivors = sum(falll_survival_t),
            spring_survivors = sum(spring_survival_t1)) 
plot(jitter(transplant.plots$weighted.dens),jitter(transplant.plots$shrub_cover))  
## This makes sense: windows with zero density always have zero plot cover,
## plots with low density often have zero plot cover,
## and plots with high density usually have non-zero plot cover

## Model selection with shrub effects at different scales
fall_surv_models <- list()
## transect effects only
fall_surv_models[[1]] <- gam(cbind(fall_survivors,4-fall_survivors) ~ 1 + s(unique.transect, bs = "re") + s(plot, bs = "re"),
                      data = transplant.plots, gamma = gamma, family = "binomial")
## window density
fall_surv_models[[2]] <- gam(cbind(fall_survivors,4-fall_survivors) ~ s(weighted.dens) + s(unique.transect, bs = "re") + s(plot, bs = "re"),
                             data = transplant.plots, gamma = gamma, family = "binomial")
## plot shrub cover
fall_surv_models[[3]] <- gam(cbind(fall_survivors,4-fall_survivors) ~ s(shrub_cover) + s(unique.transect, bs = "re") + s(plot, bs = "re"),
                             data = transplant.plots, gamma = gamma, family = "binomial")
## plot grass cover
fall_surv_models[[4]] <- gam(cbind(fall_survivors,4-fall_survivors) ~ s(grass_cover) + s(unique.transect, bs = "re") + s(plot, bs = "re"),
                             data = transplant.plots, gamma = gamma, family = "binomial")
## window density and plot shrub cover
fall_surv_models[[5]] <- gam(cbind(fall_survivors,4-fall_survivors) ~ s(weighted.dens) + s(shrub_cover) + s(unique.transect, bs = "re") + s(plot, bs = "re"),
                             data = transplant.plots, gamma = gamma, family = "binomial")
## plot shrub cover and plot grass cover
fall_surv_models[[6]] <- gam(cbind(fall_survivors,4-fall_survivors) ~ s(grass_cover) + s(shrub_cover) + s(unique.transect, bs = "re") + s(plot, bs = "re"),
                             data = transplant.plots, gamma = gamma, family = "binomial")
## window density, plot shrub and grass cover
fall_surv_models[[7]] <- gam(cbind(fall_survivors,4-fall_survivors) ~ s(weighted.dens) + s(grass_cover) + s(shrub_cover) + s(unique.transect, bs = "re") + s(plot, bs = "re"),
                             data = transplant.plots, gamma = gamma, family = "binomial")
## Model selection
AICtab(fall_surv_models)
fall_surv_fitted_terms <- predict(fall_surv_models[[2]], type = "terms") 
plot(transplant.plots$weighted.dens, fall_surv_fitted_terms[, "s(weighted.dens)"]) 

fall_surv_fitted_terms <- predict(fall_surv_models[[5]], type = "terms") 
plot(transplant.plots$weighted.dens, fall_surv_fitted_terms[, "s(weighted.dens)"]) 
plot(transplant.plots$shrub_cover, fall_surv_fitted_terms[, "s(shrub_cover)"]) 

## Same but for spring
spring_surv_models <- list()
## transect effects only
spring_surv_models[[1]] <- gam(cbind(spring_survivors,fall_survivors-spring_survivors) ~ 1 + s(unique.transect, bs = "re") + s(plot, bs = "re"),
                             data = subset(transplant.plots,fall_survivors>0), gamma = gamma, family = "binomial")
## window density
spring_surv_models[[2]] <- gam(cbind(spring_survivors,fall_survivors-spring_survivors) ~ s(weighted.dens) + s(unique.transect, bs = "re") + s(plot, bs = "re"),
                             data = subset(transplant.plots,fall_survivors>0), gamma = gamma, family = "binomial")
## plot shrub cover
spring_surv_models[[3]] <- gam(cbind(spring_survivors,fall_survivors-spring_survivors) ~ s(shrub_cover) + s(unique.transect, bs = "re") + s(plot, bs = "re"),
                             data = subset(transplant.plots,fall_survivors>0), gamma = gamma, family = "binomial")
## plot grass cover
spring_surv_models[[4]] <- gam(cbind(spring_survivors,fall_survivors-spring_survivors) ~ s(grass_cover) + s(unique.transect, bs = "re") + s(plot, bs = "re"),
                             data = subset(transplant.plots,fall_survivors>0), gamma = gamma, family = "binomial")
## window density and plot shrub cover
spring_surv_models[[5]] <- gam(cbind(spring_survivors,fall_survivors-spring_survivors) ~ s(weighted.dens) + s(shrub_cover) + s(unique.transect, bs = "re") + s(plot, bs = "re"),
                             data = subset(transplant.plots,fall_survivors>0), gamma = gamma, family = "binomial")
## plot shrub cover and plot grass cover
spring_surv_models[[6]] <- gam(cbind(spring_survivors,fall_survivors-spring_survivors) ~ s(grass_cover) + s(shrub_cover) + s(unique.transect, bs = "re") + s(plot, bs = "re"),
                             data = subset(transplant.plots,fall_survivors>0), gamma = gamma, family = "binomial")
## window density, plot shrub and grass cover
spring_surv_models[[7]] <- gam(cbind(spring_survivors,fall_survivors-spring_survivors) ~ s(weighted.dens) + s(grass_cover) + s(shrub_cover) + s(unique.transect, bs = "re") + s(plot, bs = "re"),
                             data = subset(transplant.plots,fall_survivors>0), gamma = gamma, family = "binomial")
## Model selection
AICtab(spring_surv_models)
spring_surv_fitted_terms <- predict(spring_surv_models[[6]], type = "terms") 
plot(transplant.plots$grass_cover[transplant.plots$fall_survivors>0], spring_surv_fitted_terms[, "s(grass_cover)"]) 
plot(transplant.plots$shrub_cover[transplant.plots$fall_survivors>0], spring_surv_fitted_terms[, "s(shrub_cover)"]) 

## show fitted
transplant.plots$springpred4 <- predict.gam(spring_surv_models[[4]], newdata = transplant.plots, type = "response", exclude = "s(unique.transect)")
plot(jitter(transplant.plots$grass_cover),jitter(transplant.plots$spring_survivors/4))
points(transplant.plots$grass_cover,transplant.plots$springpred4,pch=16,col="red")


## what was total fall and spring survival?
sum(transplant.plots$fall_survivors);length(transplant.plots$fall_survivors)
sum(transplant.plots$spring_survivors);sum(transplant.plots$fall_survivors)


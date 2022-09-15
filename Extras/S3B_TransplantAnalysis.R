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
fall_surv_models[[1]] <- gam(cbind(fall_survivors,4-fall_survivors) ~ 1 + s(unique.transect, bs = "re"),
                      data = transplant.plots, gamma = gamma, family = "binomial")
## window density
fall_surv_models[[2]] <- gam(cbind(fall_survivors,4-fall_survivors) ~ s(weighted.dens) + s(unique.transect, bs = "re"),
                             data = transplant.plots, gamma = gamma, family = "binomial")
## plot shrub cover
fall_surv_models[[3]] <- gam(cbind(fall_survivors,4-fall_survivors) ~ s(shrub_cover) + s(unique.transect, bs = "re"),
                             data = transplant.plots, gamma = gamma, family = "binomial")
## plot grass cover
fall_surv_models[[4]] <- gam(cbind(fall_survivors,4-fall_survivors) ~ s(grass_cover) + s(unique.transect, bs = "re"),
                             data = transplant.plots, gamma = gamma, family = "binomial")
## window density and plot shrub cover
fall_surv_models[[5]] <- gam(cbind(fall_survivors,4-fall_survivors) ~ s(weighted.dens) + s(shrub_cover) + s(unique.transect, bs = "re"),
                             data = transplant.plots, gamma = gamma, family = "binomial")
## plot shrub cover and plot grass cover
fall_surv_models[[6]] <- gam(cbind(fall_survivors,4-fall_survivors) ~ s(grass_cover) + s(shrub_cover) + s(unique.transect, bs = "re"),
                             data = transplant.plots, gamma = gamma, family = "binomial")
## window density and grass cover
fall_surv_models[[7]] <- gam(cbind(fall_survivors,4-fall_survivors) ~ s(weighted.dens) + s(grass_cover) + s(unique.transect, bs = "re"),
                             data = transplant.plots, gamma = gamma, family = "binomial")
## window density, plot shrub and grass cover
fall_surv_models[[8]] <- gam(cbind(fall_survivors,4-fall_survivors) ~ s(weighted.dens) + s(grass_cover) + s(shrub_cover) + s(unique.transect, bs = "re"),
                             data = transplant.plots, gamma = gamma, family = "binomial")

## Model selection
AICtab(fall_surv_models,weights=T)

fall_surv_fitted_terms <- predict(fall_surv_models[[7]], type = "terms") 
plot(transplant.plots$weighted.dens, fall_surv_fitted_terms[, "s(weighted.dens)"]) 
plot(transplant.plots$grass_cover, fall_surv_fitted_terms[, "s(grass_cover)"]) 

fall_surv_fitted_terms <- predict(fall_surv_models[[2]], type = "terms") 
plot(transplant.plots$weighted.dens, fall_surv_fitted_terms[, "s(weighted.dens)"]) 

fall_surv_fitted_terms <- predict(fall_surv_models[[5]], type = "terms") 
plot(transplant.plots$weighted.dens, fall_surv_fitted_terms[, "s(weighted.dens)"]) 
plot(transplant.plots$shrub_cover, fall_surv_fitted_terms[, "s(shrub_cover)"]) 

## plot fitted model
transplant.plots$fallpred2 <- predict.gam(fall_surv_models[[2]], newdata = transplant.plots, type = "response", exclude = "s(unique.transect)")
transplant.plots$fallpred5 <- predict.gam(fall_surv_models[[5]], newdata = transplant.plots, type = "response", exclude = "s(unique.transect)")
transplant.plots$fall_surv <- transplant.plots$fall_survivors/4

plot(transplant.plots$weighted.dens,transplant.plots$fall_surv)
points(transplant.plots$weighted.dens,transplant.plots$fallpred2,col="red")

plot(transplant.plots$weighted.dens,transplant.plots$fall_surv)
points(transplant.plots$weighted.dens,transplant.plots$fallpred5,col="red")

plot(transplant.plots$shrub_cover,transplant.plots$fall_surv)
points(transplant.plots$shrub_cover,transplant.plots$fallpred5,col="red")


## Same but for spring
spring_surv_models <- list()
## transect effects only
spring_surv_models[[1]] <- gam(cbind(spring_survivors,fall_survivors-spring_survivors) ~ 1 + s(unique.transect, bs = "re"),
                             data = subset(transplant.plots,fall_survivors>0), gamma = gamma, family = "binomial")
## window density
spring_surv_models[[2]] <- gam(cbind(spring_survivors,fall_survivors-spring_survivors) ~ s(weighted.dens) + s(unique.transect, bs = "re"),
                             data = subset(transplant.plots,fall_survivors>0), gamma = gamma, family = "binomial")
## plot shrub cover
spring_surv_models[[3]] <- gam(cbind(spring_survivors,fall_survivors-spring_survivors) ~ s(shrub_cover) + s(unique.transect, bs = "re"),
                             data = subset(transplant.plots,fall_survivors>0), gamma = gamma, family = "binomial")
## plot grass cover
spring_surv_models[[4]] <- gam(cbind(spring_survivors,fall_survivors-spring_survivors) ~ s(grass_cover) + s(unique.transect, bs = "re"),
                             data = subset(transplant.plots,fall_survivors>0), gamma = gamma, family = "binomial")
## window density and plot shrub cover
spring_surv_models[[5]] <- gam(cbind(spring_survivors,fall_survivors-spring_survivors) ~ s(weighted.dens) + s(shrub_cover) + s(unique.transect, bs = "re"),
                             data = subset(transplant.plots,fall_survivors>0), gamma = gamma, family = "binomial")
## plot shrub cover and plot grass cover
spring_surv_models[[6]] <- gam(cbind(spring_survivors,fall_survivors-spring_survivors) ~ s(grass_cover) + s(shrub_cover) + s(unique.transect, bs = "re"),
                             data = subset(transplant.plots,fall_survivors>0), gamma = gamma, family = "binomial")
## window density, plot shrub and grass cover
spring_surv_models[[7]] <- gam(cbind(spring_survivors,fall_survivors-spring_survivors) ~ s(weighted.dens) + s(grass_cover) + s(unique.transect, bs = "re"),
                             data = subset(transplant.plots,fall_survivors>0), gamma = gamma, family = "binomial")
## window density, plot shrub and grass cover
spring_surv_models[[8]] <- gam(cbind(spring_survivors,fall_survivors-spring_survivors) ~ s(weighted.dens) + s(grass_cover) + s(shrub_cover) + s(unique.transect, bs = "re"),
                               data = subset(transplant.plots,fall_survivors>0), gamma = gamma, family = "binomial")

## Model selection
AICtab(spring_surv_models,weights=T)
spring_surv_fitted_terms <- predict(spring_surv_models[[6]], type = "terms") 
plot(transplant.plots$grass_cover[transplant.plots$fall_survivors>0], spring_surv_fitted_terms[, "s(grass_cover)"]) 
plot(transplant.plots$shrub_cover[transplant.plots$fall_survivors>0], spring_surv_fitted_terms[, "s(shrub_cover)"]) 

## show fitted
transplant.plots$springpred4 <- predict.gam(spring_surv_models[[4]], newdata = transplant.plots, type = "response", exclude = "s(unique.transect)")
plot(jitter(transplant.plots$grass_cover),jitter(transplant.plots$spring_survivors/4))
points(transplant.plots$grass_cover,transplant.plots$springpred4,pch=16,col="red")


## what was total fall and spring survival?
sum(transplant.plots$fall_survivors);length(transplant.plots$fall_survivors)*4
sum(transplant.plots$spring_survivors);sum(transplant.plots$fall_survivors)

## figure for appendix
## to visualize model 7 I need a clean data set
model7pred_weightdens <- data.frame(weighted.dens=seq(min(transplant.plots$weighted.dens),max(transplant.plots$weighted.dens),5),
                                    grass_cover=mean(transplant.plots$grass_cover),
                                    unique.transect=transplant.plots$unique.transect[1])
model7pred_weightdens$pred <- predict.gam(fall_surv_models[[7]], newdata = model7pred_weightdens, type = "response", exclude = "s(unique.transect)")

plot(transplant.plots$weighted.dens,transplant.plots$fall_surv)
lines(model7pred_weightdens$weighted.dens,model7pred_weightdens$pred,col="red",lwd=3)

model7pred_grass <- data.frame(grass_cover=seq(min(transplant.plots$grass_cover),max(transplant.plots$grass_cover),0.01),
                               weighted.dens=mean(transplant.plots$weighted.dens),
                                    unique.transect=transplant.plots$unique.transect[1])
model7pred_grass$pred <- predict.gam(fall_surv_models[[7]], newdata = model7pred_grass, type = "response", exclude = "s(unique.transect)")

plot(transplant.plots$grass_cover,transplant.plots$fall_surv)
lines(model7pred_grass$grass_cover,model7pred_grass$pred,col="red",lwd=3)


pdf("Manuscript/Figures/survival_appendix.pdf",useDingbats = F,height=4,width=7)
par(mfrow=c(1,2))
plot(jitter(transplant.plots$weighted.dens),
     jitter(transplant.plots$fall_surv),
     ylab="July-October survival",xlab="Window weighted density")
lines(model7pred_weightdens$weighted.dens,model7pred_weightdens$pred,col="red",lwd=3)
title("A",adj=0)

plot(jitter(transplant.plots$grass_cover),
     jitter(transplant.plots$fall_surv),
     ylab="October-June survival",xlab="Local grass cover (%)")
lines(model7pred_grass$grass_cover,model7pred_grass$pred,col="red",lwd=3)
title("B",adj=0)
dev.off()

# Collect AIC table info --------------------------------------------------
fall<-as.data.frame(AICtab(fall_surv_models,weights=T,sort=F))
spring<-as.data.frame(AICtab(spring_surv_models,weights=T,sort=F))

fall$`Pr(Survival)`<- c("(1|transect)",
                                "window density + (1|transect)",
                                "shrub cover + (1|transect)",
                                "grass cover + (1|transect)",
                                "window density + shrub cover + (1|transect)",
                                "grass cover + shrub cover + (1|transect)",
                                "window density + grass cover + (1|transect)",
                                "window density + shrub cover + grass cover + (1|transect)")
spring$`Pr(Survival)`<- c("(1|transect)",
                        "window density + (1|transect)",
                        "shrub cover + (1|transect)",
                        "grass cover + (1|transect)",
                        "window density + shrub cover + (1|transect)",
                        "grass cover + shrub cover + (1|transect)",
                        "window density + grass cover + (1|transect)")

appendix_aic_tables = list(fall_aic=fall,spring_aic=spring)
write_rds(x=appendix_aic_tables,file="Manuscript/appendix_aic_tables.rds")

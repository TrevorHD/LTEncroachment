##### Initialise Data -------------------------------------------------------------------------------------

CData.Transplants.s <- CData.Transplants %>%
  
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





##### Visualise survival as function of different cover types ---------------------------------------------

# List out things that will be changed between plots
plots.xvar <- c("num_shrub_t", "num_bare_t", "total.grass", "total.plant")
plots.xlab <- c("Shrub-Covered", "Bare", "Grass-Covered", "Plant-Covered")
plots.yvar <- c("falll_survival_t", "spring_survival_t1")
plots.ylab <- c("Fall", "Spring")
plots.name <- c("Shrub", "Bare", "Grass", "AllPlants")



# Plot survival (0 or 1) and survival rate for each cover type
par(mfrow = c(2, 2))
for(i in 1:4){
  for(j in 1:2){
    plot(jitter(get(plots.xvar[i], CData.Transplants.s), amount = 0.15),
         jitter(get(plots.yvar[j], CData.Transplants.s), amount = 0.03),
         xlab = paste0(plots.xlab[i], " Sections"), ylab = paste0("Survival (", plots.ylab[j], ")"))
    mean.probs <- c()
    for(k in 0:max(get(plots.xvar[i], CData.Transplants.s))){
      mean.probs <- append(mean.probs, mean(get(plots.yvar[j], subset(CData.Transplants.s, 
                                                                      get(plots.xvar[i]) == k))))}
    plot(mean.probs, pch = 19, xlab = paste0(plots.xlab[i], " Sections"), 
         ylab = paste0("Survival Rate (", plots.ylab[j], ")"))}
  assign(paste0("SPlot.", plots.name[i]), recordPlot())}

# To see each plot, cycle through plots in the "plots" window

# Check correlation matrix between response (survival) and predictors (cover type)
cor(CData.Transplants.s[, c("spring_survival_t1", "num_bare_t", "num_shrub_t", 
                            "total.grass", "total.plant")])

# Same as above, but more visually appealing
par(mfrow = c(1, 1))
corrplot(cor(CData.Transplants.s[, c("spring_survival_t1", "num_bare_t", "num_shrub_t", 
                                     "total.grass", "total.plant")]),
         method = "square", type = "upper")

# Multicollinearity between predictors (e.g. bare patches, shrub patches, etc.) may be a problem





##### Model survival as function of different cover types -------------------------------------------------

# Univariate logistic regression
# Survival is response variable, shrub cover is predictor variable
Mod.S.T1 <- glm(spring_survival_t1 ~ num_shrub_t, data = CData.Transplants.s, family = "binomial")

# Shrub cover term is not significant
# Analysis of deviance test confirms this
summary(Mod.S.T1)
anova(Mod.S.T1, test = "Chisq")

# Poor residuals indicate that model is not valid
hatValues <- influence(Mod.S.T1)$hat
DevianceRes <- rstandard(Mod.S.T1)
PearsonRes <- residuals(Mod.S.T1, "pearson")/sqrt(1 - hatValues)
par(mfrow = c(1, 2))
plot(CData.Transplants.s$num_shrub_t, DevianceRes, 
     ylab = "Standardised Deviance Residual", ylim = c(-1, 4.5))
abline(h = 0, col = "red")
plot(CData.Transplants.s$num_shrub_t, PearsonRes, 
     ylab = "Standardised Pearson Residuals", ylim = c(-1, 4.5))
abline(h = 0, col = "red")

# This can also be seen here
plot(Mod.S.T1)

# Even with mixed effects model, shrub cover term is not significant
Mod.S.T2 <- glmer(spring_survival_t1 ~ num_shrub_t + (1 | unique.transect),
                  data = CData.Transplants.s, family = "binomial")
summary(Mod.S.T2)

# Also not significant if we try these regressions with bare ground as predictor
Mod.S.T3 <- glm(spring_survival_t1 ~ num_bare_t, data = CData.Transplants.s, family = "binomial")
summary(Mod.S.T3)
anova(Mod.S.T3, test = "Chisq")
Mod.S.T4 <- glmer(spring_survival_t1 ~ num_bare_t + (1 | unique.transect),
                  data = CData.Transplants.s, family = "binomial")
summary(Mod.S.T4)

# Use backwards variable selection
# Model with the best AIC still has a non-significant term!
Mod.S.TFull <- glm(spring_survival_t1 ~ num_bare_t + num_shrub_t + total.grass + total.plant, 
                   data = CData.Transplants.s, family = "binomial")
step(Mod.S.TFull, direction = "backward", data = data.new)
summary(glm(spring_survival_t1 ~ num_bare_t + num_shrub_t, family = "binomial", 
            data = CData.Transplants.s))

# Shrub cover patches may not be a good predictor of survival





##### Model survival as function of weighted density ------------------------------------------------------

# Plot of survival as a function of weighted density
plot(CData.Transplants.s$weighted.dens, CData.Transplants.s$spring_survival_t1)

# Univariate logistic regression
# Survival is response variable, weighted density is predictor variable
Mod.S.T5 <- glm(spring_survival_t1 ~ weighted.dens, data = CData.Transplants.s, family = "binomial")

# Weighted density term is not significant at 0.05 level
# However, analysis of deviance test says it is
summary(Mod.S.T5)
anova(Mod.S.T5, test = "Chisq")

# Poor residuals indicate that model is not valid
hatValues <- influence(Mod.S.T5)$hat
DevianceRes <- rstandard(Mod.S.T5)
PearsonRes <- residuals(Mod.S.T5, "pearson")/sqrt(1 - hatValues)
par(mfrow = c(1, 2))
plot(CData.Transplants.s$weighted.dens, DevianceRes, 
     ylab = "Standardised Deviance Residual", ylim = c(-1, 4.5))
abline(h = 0, col = "red")
plot(CData.Transplants.s$weighted.dens, PearsonRes, 
     ylab = "Standardised Pearson Residuals", ylim = c(-1, 4.5))
abline(h = 0, col = "red")

# This can also be seen here
plot(Mod.S.T5)

# Weighted density may not be a good predictor of survival
# Same goes for standardised weighted density





##### Model survival as function of shrub volume ----------------------------------------------------------

# Plot of survival as a function of volume
plot(CData.Transplants.s$volume_t, CData.Transplants.s$spring_survival_t1)

# Univariate logistic regression
# Survival is response variable, volume is predictor variable
Mod.S.T6 <- glm(spring_survival_t1 ~ volume_t, data = CData.Transplants.s, family = "binomial")

# Volume term is not significant
# Analysis of deviance test confirms this
summary(Mod.S.T6)
anova(Mod.S.T6, test = "Chisq")

# We have tried several different models
# None of them are significant, or residuals are so bad that the model is not valid
# The most likely case is that survival is always low for transplants and/or small shrubs
# This is what we found in S1_SupportingMaterial


# transplant growth -------------------------------------------------------
## only 20 spring survivors, so not much we can do with growth
sum(!is.na(CData.Transplants.s$logGR))

plot(CData.Transplants.s$d.stand,CData.Transplants.s$logGR)
plot(CData.Transplants.s$num_bare_t,CData.Transplants.s$logGR)


# Tom's survival analysis -------------------------------------------------
## I think all the local cover counts are too correlated with each other to include in the same model
## I will try models with bare, shrub, and grass at the smaller scale, and also\
## weighted density at the window scale to test vegetation predictors
## also going to use fall survival (more data to analyze)

## first try no random effects, since the PDC planting makes the transects so different
fall_surv_glm <- list()
## null model
fall_surv_glm[[1]] <- glm(falll_survival_t ~ 1, family="binomial", data=CData.Transplants.s)
## window-scale density effect
fall_surv_glm[[2]] <- glm(falll_survival_t ~ d.stand, family="binomial", data=CData.Transplants.s)
## plot scale bare
fall_surv_glm[[3]] <- glm(falll_survival_t ~ num_bare_t, family="binomial", data=CData.Transplants.s)
## plot scale grass
fall_surv_glm[[4]] <- glm(falll_survival_t ~ total.grass, family="binomial", data=CData.Transplants.s)
## plot scale shrub
fall_surv_glm[[5]] <- glm(falll_survival_t ~ num_shrub_t, family="binomial", data=CData.Transplants.s)

AICctab(fall_surv_glm,weights=T)
## support for model 2 - visualize with binned means

fallsurv_viz_bin <- CData.Transplants.s %>% 
  mutate(dens_bin = as.integer(cut_interval(d.stand,5))) %>% 
  group_by(dens_bin) %>% 
  summarise(mean_dens = mean(d.stand),
            mean_surv = mean(falll_survival_t),
            bin_n = n())

plot(CData.Transplants.s$d.stand,CData.Transplants.s$falll_survival_t,pch="|",col="gray",
     xlab="Weighted shrub density",ylab="Survival")
points(fallsurv_viz_bin$mean_dens,fallsurv_viz_bin$mean_surv,pch=16,cex=2)
lines(seq(min(CData.Transplants.s$d.stand),max(CData.Transplants.s$d.stand),length.out = 100),
      invlogit(coef(fall_surv_glm[[2]])[1] + coef(fall_surv_glm[[2]])[2] * seq(min(CData.Transplants.s$d.stand),max(CData.Transplants.s$d.stand),length.out = 100)),
      lwd=2)

## see if this holds up with spring survival (fewer overall survivors)
spring_surv_glm <- list()
spring_surv_glm[[1]] <- glm(spring_survival_t1 ~ 1, family="binomial", data=CData.Transplants.s)
spring_surv_glm[[2]] <- glm(spring_survival_t1 ~ d.stand, family="binomial", data=CData.Transplants.s)
spring_surv_glm[[3]] <- glm(spring_survival_t1 ~ num_bare_t, family="binomial", data=CData.Transplants.s)
spring_surv_glm[[4]] <- glm(spring_survival_t1 ~ total.grass, family="binomial", data=CData.Transplants.s)
spring_surv_glm[[5]] <- glm(spring_survival_t1 ~ num_shrub_t, family="binomial", data=CData.Transplants.s)
AICctab(spring_surv_glm,weights=T)
## still support for model 2 

springsurv_viz_bin <- CData.Transplants.s %>% 
  mutate(dens_bin = as.integer(cut_interval(d.stand,5))) %>% 
  group_by(dens_bin) %>% 
  summarise(mean_dens = mean(d.stand),
            mean_surv = mean(spring_survival_t1),
            bin_n = n())

plot(CData.Transplants.s$d.stand,CData.Transplants.s$spring_survival_t1,pch="|",col="gray",
     xlab="Weighted shrub density",ylab="Survival")
points(springsurv_viz_bin$mean_dens,springsurv_viz_bin$mean_surv,pch=16,cex=2)
lines(seq(min(CData.Transplants.s$d.stand),max(CData.Transplants.s$d.stand),length.out = 100),
      invlogit(coef(spring_surv_glm[[2]])[1] + coef(spring_surv_glm[[2]])[2] * seq(min(CData.Transplants.s$d.stand),max(CData.Transplants.s$d.stand),length.out = 100)),
      lwd=2)

## can I fit model with random effects, and are results the same?
fall_surv_rfx <- list()
fall_surv_rfx[[1]] <- glmer(falll_survival_t ~ 1 + (1|unique.transect), family="binomial", data=CData.Transplants.s)
fall_surv_rfx[[2]] <- glmer(falll_survival_t ~ d.stand + (1|unique.transect), family="binomial", data=CData.Transplants.s)
fall_surv_rfx[[3]] <- glmer(falll_survival_t ~ num_bare_t + (1|unique.transect), family="binomial", data=CData.Transplants.s)
fall_surv_rfx[[4]] <- glmer(falll_survival_t ~ total.grass + (1|unique.transect), family="binomial", data=CData.Transplants.s)
fall_surv_rfx[[5]] <- glmer(falll_survival_t ~ num_shrub_t + (1|unique.transect), family="binomial", data=CData.Transplants.s)
AICctab(fall_surv_rfx,weights=T)

plot(CData.Transplants.s$d.stand,CData.Transplants.s$falll_survival_t,pch="|",col="gray",
     xlab="Weighted shrub density",ylab="Survival")
points(fallsurv_viz_bin$mean_dens,fallsurv_viz_bin$mean_surv,pch=16,cex=2)
lines(seq(min(CData.Transplants.s$d.stand),max(CData.Transplants.s$d.stand),length.out = 100),
      invlogit(fixef(fall_surv_rfx[[2]])[1] + fixef(fall_surv_rfx[[2]])[2] * seq(min(CData.Transplants.s$d.stand),max(CData.Transplants.s$d.stand),length.out = 100)),
      lwd=2)
## note that mean survival is much lower in the mixed model, and this is good because it means
## the migh survival at PDC is being attributed to random effects, such that the overall
## survival probablility at an average site is reduced

## finally, spring survival with rfx
spring_surv_rfx <- list()
spring_surv_rfx[[1]] <- glmer(spring_survival_t1 ~ 1 + (1|unique.transect), family="binomial", data=CData.Transplants.s)
spring_surv_rfx[[2]] <- glmer(spring_survival_t1 ~ d.stand + (1|unique.transect), family="binomial", data=CData.Transplants.s)
spring_surv_rfx[[3]] <- glmer(spring_survival_t1 ~ num_bare_t + (1|unique.transect), family="binomial", data=CData.Transplants.s)
spring_surv_rfx[[4]] <- glmer(spring_survival_t1 ~ total.grass + (1|unique.transect), family="binomial", data=CData.Transplants.s)
spring_surv_rfx[[5]] <- glmer(spring_survival_t1 ~ num_shrub_t + (1|unique.transect), family="binomial", data=CData.Transplants.s)
AICctab(spring_surv_rfx,weights=T)
## interestingly, these rankings provide mixed support for models 3,4,2
## 3 has a negative effect of bare ground, 4 has a positive effect of grass,
## and 2 has a negative effect of weighted shrub density at the window scale
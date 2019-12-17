##### Initialise Data -------------------------------------------------------------------------------------

# Transform density; calculate total grass cover and plant cover
CData.Transplants.s <- CData.Transplants %>% 
  mutate(d.stand = (weighted.dens - mean(weighted.dens, na.rm = TRUE)) / sd(weighted.dens, na.rm = TRUE),
         total.grass = num_black_gramma_t + num_blue_gramma_t + num_other_grass_t,
         total.plant = num_shrub_t + total.grass + num_other_t,
         unique.transect = interaction(transect, site))





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
      mean.probs <- append(mean.probs, mean(get(plots.yvar[j], subset(CData.Transplants.s, get(plots.xvar[i]) == k))))}
    plot(mean.probs, pch = 19, xlab = paste0(plots.xlab[i], " Sections"), 
         ylab = paste0("Survival Rate (", plots.ylab[j], ")"))}
  assign(paste0("SPlot.", plots.name[i]), recordPlot())}

# To see each plot, cycle through plots in the "plots" window

# Check correlation matrix between response (survival) and predictors (cover type)
cor(CData.Transplants.s[, c("spring_survival_t1", "num_bare_t", "num_shrub_t", "total.grass", "total.plant")])

# Same as above, but more visually appealing
corrplot(cor(CData.Transplants.s[, c("spring_survival_t1", "num_bare_t", "num_shrub_t", "total.grass", "total.plant")]),
         method = "square", type = "upper")

# Multicollinearity between predictors (e.g. bare patches, shurb patches, etc.) may be a problem





##### Model survival as function of different cover types -------------------------------------------------

# Univariate logistic regression
# Survival is response variable, shrub cover is predictor variable
Mod.S.T1 <- glm(falll_survival_t ~ num_shrub_t, data = CData.Transplants.s, family = "binomial")

# Shrub cover term is not significant
# Analysis of deviance test confirms this
summary(Mod.S.T1)
anova(Mod.S.T1, test = "Chisq")

# Poor residuals indicate than model is not valid
hatValues <- influence(Mod.S.T1)$hat
DevianceRes <- rstandard(Mod.S.T1)
PearsonRes <- residuals(Mod.S.T1, "pearson")/sqrt(1 - hatValues)
par(mfrow = c(1, 2))
plot(CData.Transplants.s$num_shrub_t, DevianceRes, 
     ylab = "Standardised Deviance Residual", ylim = c(-2.5, 2.5))
abline(h = 0, col = "red")
plot(CData.Transplants.s$num_shrub_t, PearsonRes, 
     ylab = "Standardised Pearson Residuals", ylim = c(-2.5, 2.5))
abline(h = 0, col = "red")

# This can also be seen here
plot(Mod.S.T1)

# Even with mixed effects model, shrub cover term is not significant
Mod.S.T2 <- glmer(falll_survival_t ~ num_shrub_t + (1 | unique.transect),
                  data = CData.Transplants.s, family = "binomial")
summary(Mod.S.T2)

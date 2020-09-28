##### Tidy up transplant data for survival analysis -------------------------------------------------------



# Calculate conical volume; use initial volume since most plants die after a year
# Add variable indicating which plants are transplants
# Highest survival occurs at PDC, so we will use that as the best-case scenario
boot.CData.Transplants %>% 
  filter(site == "PDC") -> boot.CData.Transplants

# Combine transplants with large shrubs for later survival analysis
# Keep only location info, survival, volume, and density
select(boot.CData.Transplants, "site", "transect", "actual.window", 
       "spring_survival_t1", "volume_t", "weighted.dens", "transplant") %>% 
  rename("survival_t1" = "spring_survival_t1") %>% 
  rbind(select(boot.CData, "site", "transect", "actual.window", 
               "survival_t1", "volume_t", "weighted.dens", "transplant")) -> CData.AllSurvival





##### Create GLM for survival (transplants included) ------------------------------------------------------

# Remove entries for which there is no recorded volume or survival
# Standardise density; introduce unique transect identifier
CData.AllSurvival.s <- CData.AllSurvival %>% 
  drop_na(volume_t, survival_t1) %>% 
  mutate(d.stand = (weighted.dens - mean(weighted.dens, na.rm = TRUE)) / sd(weighted.dens, na.rm = TRUE),
         unique.transect = interaction(transect, site),
         survival_t1 = as.numeric(survival_t1))

# Create a list of possible survival models
Mod.S <- list()

# "Null" model; only random effects of site and transect within site
Mod.S[[1]] <- glmer(survival_t1 ~ (1 | unique.transect),
                    data = subset(CData.AllSurvival.s, transplant == TRUE), family = "binomial")

# Size-only model
Mod.S[[2]] <- glmer(survival_t1 ~ volume_t + (1|unique.transect),
                    data = subset(CData.AllSurvival.s, transplant == TRUE), family = "binomial")

# Density-only model
Mod.S[[3]] <- glmer(survival_t1 ~ d.stand + (1|unique.transect),
                    data = subset(CData.AllSurvival.s, transplant == TRUE), family = "binomial")

# Size and density (additive)
Mod.S[[4]] <- glmer(survival_t1 ~ volume_t + d.stand + (1|unique.transect),
                    data = subset(CData.AllSurvival.s, transplant == TRUE), family = "binomial")

# Size and density (interactive)
Mod.S[[5]] <- glmer(survival_t1 ~ volume_t * d.stand + (1|unique.transect),
                    data = subset(CData.AllSurvival.s, transplant == TRUE), family = "binomial")

# Density-only model (quadratic)
Mod.S[[6]] <- glmer(survival_t1 ~ d.stand + I(d.stand^2) + (1|unique.transect),
                    data = subset(CData.AllSurvival.s, transplant == TRUE), family = "binomial")

# Size (linear) and density (quadratic)
Mod.S[[7]] <- glmer(survival_t1 ~ volume_t + d.stand + I(d.stand^2) + (1|unique.transect),
                    data = subset(CData.AllSurvival.s, transplant == TRUE), family = "binomial")

# Size and density (interactive, quadratic)
Mod.S[[8]] <- glmer(survival_t1 ~ volume_t * d.stand + volume_t * I(d.stand^2) + (1|unique.transect),
                    data = subset(CData.AllSurvival.s, transplant == TRUE), family = "binomial")

# Calculate an AIC table, ranked from best to worst model
# Weights interpreted as the proportion of evidence in favour of each
# To do this, use AICtab(Mod.S, weights = TRUE, sort = TRUE)

# Keep all models in which the weight is greater than 0.001
Mod.S.AIC <- AICtab(Mod.S[[3]], Mod.S[[4]], Mod.S[[6]], Mod.S[[1]], 
                    Mod.S[[5]], Mod.S[[2]], Mod.S[[7]], Mod.S[[8]], 
                    weights = TRUE, sort = FALSE)

# Create vector of average GLM coefficients from the kept models
Mod.S.avg.cf <- c()

# Intercept
Mod.S.avg.cf[1] <- 
  Mod.S.AIC$weight[1]*fixef(Mod.S[[3]])["(Intercept)"] + 
  Mod.S.AIC$weight[2]*fixef(Mod.S[[4]])["(Intercept)"] + 
  Mod.S.AIC$weight[3]*fixef(Mod.S[[6]])["(Intercept)"] +
  Mod.S.AIC$weight[4]*fixef(Mod.S[[1]])["(Intercept)"] +
  Mod.S.AIC$weight[5]*fixef(Mod.S[[5]])["(Intercept)"] +
  Mod.S.AIC$weight[6]*fixef(Mod.S[[2]])["(Intercept)"] +
  Mod.S.AIC$weight[7]*fixef(Mod.S[[7]])["(Intercept)"] +
  Mod.S.AIC$weight[8]*fixef(Mod.S[[8]])["(Intercept)"]

# Volume coefficient
Mod.S.avg.cf[2] <- 
  Mod.S.AIC$weight[2]*fixef(Mod.S[[4]])["volume_t"] + 
  Mod.S.AIC$weight[5]*fixef(Mod.S[[5]])["volume_t"] +
  Mod.S.AIC$weight[6]*fixef(Mod.S[[2]])["volume_t"] + 
  Mod.S.AIC$weight[7]*fixef(Mod.S[[7]])["volume_t"] +
  Mod.S.AIC$weight[8]*fixef(Mod.S[[8]])["volume_t"]

# Density coefficient
Mod.S.avg.cf[3] <- 
  Mod.S.AIC$weight[1]*fixef(Mod.S[[3]])["d.stand"] + 
  Mod.S.AIC$weight[2]*fixef(Mod.S[[4]])["d.stand"] + 
  Mod.S.AIC$weight[3]*fixef(Mod.S[[6]])["d.stand"] +
  Mod.S.AIC$weight[5]*fixef(Mod.S[[5]])["d.stand"] +
  Mod.S.AIC$weight[7]*fixef(Mod.S[[7]])["d.stand"] + 
  Mod.S.AIC$weight[8]*fixef(Mod.S[[8]])["d.stand"]

# Volume and density interaction coefficient
Mod.S.avg.cf[4] <- 
  Mod.S.AIC$weight[5]*fixef(Mod.S[[5]])["volume_t:d.stand"] +
  Mod.S.AIC$weight[8]*fixef(Mod.S[[8]])["volume_t:d.stand"]

# Density quadratic coefficient
Mod.S.avg.cf[5] <- 
  Mod.S.AIC$weight[3]*fixef(Mod.S[[6]])["I(d.stand^2)"] +
  Mod.S.AIC$weight[7]*fixef(Mod.S[[7]])["I(d.stand^2)"] +
  Mod.S.AIC$weight[8]*fixef(Mod.S[[8]])["I(d.stand^2)"]

# Volume and quadratic density interaction coefficient
Mod.S.avg.cf[6] <- 
  Mod.S.AIC$weight[8]*fixef(Mod.S[[8]])["volume_t:I(d.stand^2)"]


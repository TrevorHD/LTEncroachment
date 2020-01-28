##### Create GLM for flowering probability ----------------------------------------------------------------

# Define inverse logit function used as link to binomial distribution
invlogit <- function(x){exp(x)/(1 + exp(x))} 

# Remove entries for which there is no recorded volume
# Then standardise density and introduce unique transect identifier
CData.s <- CData %>% 
  drop_na(volume_t) %>% 
  mutate(d.stand = (weighted.dens - mean(weighted.dens, na.rm = TRUE)) / sd(weighted.dens, na.rm = TRUE),
         unique.transect = interaction(transect, site))

# Create a list of possible flowering models
Mod.F <- list()

  # "Null" model; only random effects of site and transect within site
  Mod.F[[1]] <- glmer(did.flower ~ (1 | unique.transect),
                      data = CData.s, family = "binomial")
  
  # Size-only model
  Mod.F[[2]] <- glmer(did.flower ~ volume_t + (1 | unique.transect),
                      data = CData.s, family = "binomial")
  
  # Density-only model
  Mod.F[[3]] <- glmer(did.flower ~ d.stand + (1 | unique.transect),
                      data = CData.s, family = "binomial")
  
  # Size and density (additive)
  Mod.F[[4]] <- glmer(did.flower ~ volume_t + d.stand + (1 | unique.transect),
                      data = CData.s, family = "binomial")
  
  # Size and density (interactive)
  Mod.F[[5]] <- glmer(did.flower ~ volume_t * d.stand + (1 | unique.transect),
                      data = CData.s, family = "binomial")
  
  # Density-only model (quadratic)
  Mod.F[[6]] <- glmer(did.flower ~ d.stand + I(d.stand^2) + (1 | unique.transect),
                      data = CData.s, family = "binomial")
  
  # Size (linear) and density (quadratic)
  Mod.F[[7]] <- glmer(did.flower ~ volume_t + d.stand + I(d.stand^2) + (1 | unique.transect),
                      data = CData.s, family = "binomial")
  
  # Size and density (interactive, quadratic)
  Mod.F[[8]] <- glmer(did.flower ~ volume_t * d.stand + volume_t * I(d.stand^2) + (1 | unique.transect),
                      data = CData.s, family = "binomial")

# Calculate an AIC table, ranked from best to worst model
# Weights interpreted as the proportion of evidence in favour of each
# To do this, use AICtab(Mod.F, weights = TRUE, sort = TRUE)
  
# Model 8 has the best AIC
# For more info, use summary(Mod.F[[8]])

# Create vector of coefficients for model with the best AIC (Model 8)
Mod.F.top.cf <- c()

  # Intercept
  Mod.F.top.cf[1] <- fixef(Mod.F[[8]])["(Intercept)"]

  # Volume coefficient
  Mod.F.top.cf[2] <- fixef(Mod.F[[8]])["volume_t"]
  
  # Density coefficient
  Mod.F.top.cf[3] <- fixef(Mod.F[[8]])["d.stand"]

  # Volume and density interaction coefficient
  Mod.F.top.cf[4] <- fixef(Mod.F[[8]])["volume_t:d.stand"]

  # Density quadratic coefficient
  Mod.F.top.cf[5] <- fixef(Mod.F[[8]])["I(d.stand^2)"]

  # Volume and quadratic density interaction coefficient
  Mod.F.top.cf[6] <- fixef(Mod.F[[8]])["volume_t:I(d.stand^2)"]

  



##### Create GLM for growth ratio -------------------------------------------------------------------------

# Restore CData.s since we deleted stuff earlier
# Remove entries for which there is no recorded volume and log growth is NA
# Then standardise density and introduce unique transect identifier
CData.s <- CData %>% 
  drop_na(volume_t, logGR) %>% 
  mutate(d.stand = (weighted.dens - mean(weighted.dens, na.rm = TRUE)) / sd(weighted.dens, na.rm = TRUE),
         unique.transect = interaction(transect, site))

# Create a list of possible growth models
Mod.G <- list()

  # "Null" model; only random effects of site and transect within site
  Mod.G[[1]] <- lmer(logGR ~ (1 | unique.transect),
                     data = CData.s, REML = FALSE)

  # Size-only model
  Mod.G[[2]] <- lmer(logGR ~ volume_t + (1 | unique.transect),
                     data = CData.s, REML = FALSE)

  # Density-only model
  Mod.G[[3]] <- lmer(logGR ~ d.stand + (1 | unique.transect),
                     data = CData.s, REML = FALSE)

  # Size and density (additive)
  Mod.G[[4]] <- lmer(logGR ~ volume_t + d.stand + (1 | unique.transect),
                     data = CData.s, REML = FALSE)

  # Size and density (interactive)
  Mod.G[[5]] <- lmer(logGR ~ volume_t * d.stand + (1 | unique.transect),
                     data = CData.s, REML = FALSE)

  # Density-only model (quadratic)
  Mod.G[[6]] <- lmer(logGR ~ d.stand + I(d.stand^2) + (1 | unique.transect),
                     data = CData.s, REML = FALSE)

  # Size (linear) and density (quadratic)
  Mod.G[[7]] <- lmer(logGR ~ volume_t + d.stand + I(d.stand^2) + (1 | unique.transect),
                     data = CData.s, REML = FALSE)

  # Size and density (interactive, quadratic)
  Mod.G[[8]] <- lmer(logGR ~ volume_t * d.stand + volume_t * I(d.stand^2) + (1 | unique.transect),
                     data = CData.s, REML = FALSE)

# Calculate an AIC table, ranked from best to worst model
# Weights interpreted as the proportion of evidence in favour of each
# To do this, use AICtab(Mod.G, weights = TRUE, sort = TRUE)

# Model 7 has the best AIC
# For more info, use summary(Mod.G[[7]])
  
# Create vector of coefficients for model with the best AIC (Model 7)
Mod.G.top.cf <- c()

  # Intercept
  Mod.G.top.cf[1] <- fixef(Mod.G[[7]])["(Intercept)"]

  # Volume coefficient
  Mod.G.top.cf[2] <- fixef(Mod.G[[7]])["volume_t"]

  # Density coefficient
  Mod.G.top.cf[3] <- fixef(Mod.G[[7]])["d.stand"]
  
  # Volume and density interaction coefficient
  Mod.G.top.cf[4] <- 0
  
  # Density quadratic coefficient
  Mod.G.top.cf[5] <- fixef(Mod.G[[7]])["I(d.stand^2)"]
    
  # Volume and quadratic density interaction coefficient
  Mod.G.top.cf[6] <- 0





##### Create GLM for number of reproductive structures ----------------------------------------------------

# Restore CData.s since we deleted stuff earlier
# Remove entries for which there is no recorded volume or that did not flower
# Then standardise density and introduce unique transect identifier
CData.s <- CData %>%
  drop_na(volume_t) %>%
  filter(did.flower == 1) %>% 
  mutate(d.stand = (weighted.dens - mean(weighted.dens, na.rm = TRUE)) / sd(weighted.dens, na.rm = TRUE),
         unique.transect = interaction(transect, site))

# Create a list of possible reproduction models
Mod.R <- list()

  # "Null" model; only random effects of site and transect within site
  Mod.R[[1]] <- lmer(logTR1 ~ (1 | unique.transect),
                     data = CData.s, REML = FALSE)

  # Size-only model
  Mod.R[[2]] <- lmer(logTR1 ~ volume_t + (1 | unique.transect),
                     data = CData.s, REML = FALSE)

  # Density-only model
  Mod.R[[3]] <- lmer(logTR1 ~ d.stand + (1 | unique.transect),
                     data = CData.s, REML = FALSE)

  # Size and density (additive)
  Mod.R[[4]] <- lmer(logTR1 ~ volume_t + d.stand + (1 | unique.transect),
                     data = CData.s, REML = FALSE)

  # Size and density (interactive)
  Mod.R[[5]] <- lmer(logTR1 ~ volume_t * d.stand + (1 | unique.transect),
                     data = CData.s, REML = FALSE)

  # Density-only model (quadratic)
  Mod.R[[6]] <- lmer(logTR1 ~ d.stand + I(d.stand^2) + (1 | unique.transect),
                     data = CData.s, REML = FALSE)

  # Size (linear) and density (quadratic)
  Mod.R[[7]] <- lmer(logTR1 ~ volume_t + d.stand + I(d.stand^2) + (1 | unique.transect),
                     data = CData.s, REML = FALSE)

  # Size and density (interactive, quadratic)
  Mod.R[[8]] <- lmer(logTR1 ~ volume_t * d.stand + volume_t * I(d.stand^2) + (1 | unique.transect),
                     data = CData.s, REML = FALSE)

# Calculate an AIC table, ranked from best to worst model
# Weights interpreted as the proportion of evidence in favour of each
# To do this, use AICtab(Mod.R,weights = TRUE, sort = TRUE)

# Model 7 has the best AIC
# For more info, use summary(Mod.R[[7]])
  
# Create vector of coefficients for model with the best AIC (Model 7)
Mod.R.top.cf <- c()

# Intercept
  Mod.R.top.cf[1] <- fixef(Mod.R[[7]])["(Intercept)"]
  
  # Volume coefficient
  Mod.R.top.cf[2] <- fixef(Mod.R[[7]])["volume_t"]

  # Density coefficient
  Mod.R.top.cf[3] <- fixef(Mod.R[[7]])["d.stand"]

  # Volume and density interaction coefficient
  Mod.R.top.cf[4] <- 0

  # Density quadratic coefficient
  Mod.R.top.cf[5] <- fixef(Mod.R[[7]])["I(d.stand^2)"]
  
  # Volume and quadratic density interaction coefficient
  Mod.R.top.cf[6] <- 0





##### Create GLM for survival (transplants only) ----------------------------------------------------------

# Use models for transplants (small individuals) only
# Survival for larger individuals is pretty much guaranteed
  
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
Mod.S.AIC <- AICtab(Mod.S[[3]], Mod.S[[1]], Mod.S[[6]], Mod.S[[4]], 
                    Mod.S[[5]], Mod.S[[2]], Mod.S[[7]], Mod.S[[8]], 
                    weights = TRUE, sort = FALSE)

# Create vector of average GLM coefficients from the kept models
Mod.S.avg.cf <- c()

  # Intercept
  Mod.S.avg.cf[1] <- 
    Mod.S.AIC$weight[1]*fixef(Mod.S[[3]])["(Intercept)"] + 
    Mod.S.AIC$weight[2]*fixef(Mod.S[[1]])["(Intercept)"] + 
    Mod.S.AIC$weight[3]*fixef(Mod.S[[6]])["(Intercept)"] +
    Mod.S.AIC$weight[4]*fixef(Mod.S[[4]])["(Intercept)"] +
    Mod.S.AIC$weight[5]*fixef(Mod.S[[5]])["(Intercept)"] +
    Mod.S.AIC$weight[6]*fixef(Mod.S[[2]])["(Intercept)"] +
    Mod.S.AIC$weight[7]*fixef(Mod.S[[7]])["(Intercept)"] +
    Mod.S.AIC$weight[8]*fixef(Mod.S[[8]])["(Intercept)"]

  # Volume coefficient
  Mod.S.avg.cf[2] <- 
    Mod.S.AIC$weight[4]*fixef(Mod.S[[4]])["volume_t"] + 
    Mod.S.AIC$weight[5]*fixef(Mod.S[[5]])["volume_t"] +
    Mod.S.AIC$weight[6]*fixef(Mod.S[[2]])["volume_t"] + 
    Mod.S.AIC$weight[7]*fixef(Mod.S[[7]])["volume_t"] +
    Mod.S.AIC$weight[8]*fixef(Mod.S[[8]])["volume_t"]
    
  # Density coefficient
  Mod.S.avg.cf[3] <- 
    Mod.S.AIC$weight[1]*fixef(Mod.S[[3]])["d.stand"] + 
    Mod.S.AIC$weight[3]*fixef(Mod.S[[6]])["d.stand"] +
    Mod.S.AIC$weight[4]*fixef(Mod.S[[4]])["d.stand"] + 
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





##### Create LM for probability of recruitment from seed --------------------------------------------------

# Here we will use the baseline plants whose reproduction does not change each year
# This is much more stable and has fewer pitfalls than using yearly seed production from demography data
  
# Merge windows in baseline transect data with their respective weighted densities
CData.Recruitment <- merge(CData.Transects, Windows, 
                           by.x = c("site", "transect", "window"),
                           by.y = c("site", "transect", "window"))
  
# Use the reproduction model to calculate number of reproductive structures in a single year for each plant
mutate(CData.Recruitment,
       d.stand = (weighted.dens - mean(weighted.dens, na.rm = TRUE)) / sd(weighted.dens, na.rm = TRUE),
       seeds = exp(Mod.R.top.cf[1] + Mod.R.top.cf[2]*volume + Mod.R.top.cf[3]*d.stand + 
                   Mod.R.top.cf[5]*(d.stand^2))) -> CData.Recruitment

# Calculate number of seeds in a single year for each 5-m window
for(i in 1:nrow(CData.Recruitment)){
  CData.Recruitment$seeds.new[i] <- sum(CData.Recruitment$seeds[CData.Recruitment$window == CData.Recruitment$window[i] &
                                                                  CData.Recruitment$transect == CData.Recruitment$transect[i] &
                                                                  CData.Recruitment$site == CData.Recruitment$site[i]], na.rm = T)}

# Select only one instance of each unique combination of site, transect, and window, then merge with CData
# We're doing this because we're interested in total seeds in each window, not seeds per plant in each window
distinct(select(CData.Recruitment, "site", "transect", "window", "seeds.new")) %>% 
  merge(CData, ., by.x = c("site", "transect", "actual.window"),
        by.y = c("site", "transect", "window")) -> CData.Recruitment

CData.Recruitment <- mutate(CData.Recruitment, recruit.prob.1y = recruits.1y/(seeds.new),
                                               recruit.prob.4y = recruits.4y/(4*seeds.new))

CData.Recruitment %>% 
  mutate(d.stand = (weighted.dens - mean(weighted.dens, na.rm = TRUE)) / sd(weighted.dens, na.rm = TRUE)) %>% 
  group_by(site, transect, actual.window, year_t1) %>% 
  select(recruit.prob.1y, recruit.prob.4y, d.stand) %>% 
  summarise(recruit.prob.1y = unique(recruit.prob.1y),
            recruit.prob.4y = unique(recruit.prob.4y),
            d.stand = unique(d.stand)) -> CData3

plot(CData3$d.stand, CData3$recruit.prob.1y)
  
# Restore CData.s since we deleted stuff earlier, then standardise density
# CData.s <- CData %>%
  # mutate(d.stand = (weighted.dens - mean(weighted.dens, na.rm = TRUE)) / sd(weighted.dens, na.rm = TRUE))

# Collapse data to give unique combinations of site, transect, and window
# This helps avoid inflation of data
# Create simple linear model, using the PMod LM for the seed counts used in calculating the probabilities
# CData.s %>% 
  # group_by(site, transect, actual.window) %>% 
  # select(recruit.prob, d.stand) %>% 
  # summarise(recruit.prob = unique(recruit.prob),
  #           d.stand = unique(d.stand)) %>% 
  # lm(recruit.prob ~ d.stand, data = .) -> Mod.P

# Same as above, but with outliers removed
# subset(CData.s, recruit.prob < 0.002) %>% 
#   group_by(site, transect, actual.window, year_t) %>% 
#   select(recruit.prob, d.stand) %>% 
#   summarise(recruit.prob = unique(recruit.prob),
#             d.stand = unique(d.stand)) %>% 
#   glm(recruit.prob ~ d.stand, data = ., family = "quasibinomial") -> Mod.P

# Use logistic regression instead; will have to change 06_SIPM to make this work
# subset(CData.s, recruit.prob < 0.002) %>% 
#   group_by(site, transect, actual.window, year_t) %>% 
#   select(recruit.prob, d.stand) %>% 
#   summarise(recruit.prob = unique(recruit.prob),
#             d.stand = unique(d.stand)) %>% 
#   glm(recruit.prob ~ d.stand, data = ., family = "quasibinomial") -> Mod.P

#plot(CData.s$d.stand,CData.s$recruit.prob)
#abline(coef(Mod.P))
#anova(Mod.P)

# Create vector of coefficients for this model
# Mod.P.cf <- c()
  
  # Intercept
  # Mod.P.cf[1] <- coef(Mod.P)[1]
  
  # Volume coefficient
  # Mod.P.cf[2] <- coef(Mod.P)[2]

# Less complicated overall compared to the previous models


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

# Create vector of coefficients for model with the best AIC
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

# Create vector of coefficients for model with the best AIC
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

# Create vector of coefficients for model with the best AIC
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

# We already have the models and averaged model
# The rest of this section is simply visualising the data

# Calculate range of standardised data
CData.s.surv <- subset(CData.AllSurvival.s, transplant == TRUE)
v.range <- range(CData.s.surv$volume_t)
d.range <- range(CData.s.surv$d.stand)

# Calculate predictions for 4 "bins" of volume categories
v.cut <- CData.s.surv %>% 
  mutate(v.bin = cut_number(volume_t, n = 4)) %>% 
  group_by(v.bin) %>% 
  summarise(v.mean = mean(volume_t))
d.and.v <- as.tibble(expand.grid(v.cut$v.mean, seq(d.range[1], d.range[2], 0.5)))
names(d.and.v) <- c("v.mean", "x.d") 
S.pred <- d.and.v %>% 
  mutate(S.mean = invlogit(Mod.S.avg.cf[1] + 
                           Mod.S.avg.cf[2] * v.mean + 
                           Mod.S.avg.cf[3] * x.d +
                           Mod.S.avg.cf[4] * v.mean * x.d +
                           Mod.S.avg.cf[5] * (x.d^2) +
                           Mod.S.avg.cf[6] * v.mean * (x.d^2)))

# Get binned means of the data with respect to volume and density
S.mean.df <- CData.s.surv %>% 
  mutate(v.bin = cut_number(volume_t, n = 4),
         d.bin = cut_interval(d.stand, n = 6)) %>% 
  group_by(d.bin, v.bin) %>% 
  summarise(S.mean = mean(survival_t1),
            d.mean = mean(d.stand),
            v.mean = mean(volume_t))

# Generate visualisation; survival as function of density, with 4 volume bins
S.pred %>% 
  ggplot() +
  geom_line(aes(x = x.d, y = S.mean, colour = as.factor(v.mean)), size = 1) +
  scale_colour_manual(values = c("darkorchid2", "dodgerblue1", "navy", "firebrick2",
                                 "firebrick2", "darkorchid2", "dodgerblue1", "navy")) +
  geom_point(data = S.mean.df, aes(x = d.mean, y = S.mean, 
                                   colour = as.factor(v.bin)), size = 3) +
  labs(x = "Standardised Weighted Density (d)", y = "Survival") +
  ylim(-0.01, 0.04) +
  theme_bw() -> Mod.S.avg.dPlot

# Do same process as previous 3 steps, but switch density and volume
# This produces a graph of reproductive count as a function of volume, with 4 density bins
d.cut <- CData.s.surv %>% 
  mutate(d.bin = cut_interval(d.stand, n = 4)) %>% 
  group_by(d.bin) %>% 
  summarise(d.mean = mean(d.stand))
d.and.v.2 <- as.tibble(expand.grid(d.cut$d.mean, seq(v.range[1], v.range[2], 0.5)))
names(d.and.v.2) <- c("d.mean", "x.v") 
S.pred <- d.and.v.2 %>% 
  mutate(S.mean = invlogit(Mod.S.avg.cf[1] + 
                           Mod.S.avg.cf[2] * x.v + 
                           Mod.S.avg.cf[3] * d.mean +
                           Mod.S.avg.cf[4] * x.v * d.mean +
                           Mod.S.avg.cf[5] * (d.mean^2) +
                           Mod.S.avg.cf[6] * x.v * (d.mean^2)))
S.mean.df <- CData.s.surv %>% 
  mutate(d.bin = cut_interval(d.stand, n = 4),
         v.bin = cut_number(volume_t, n = 6)) %>% 
  group_by(v.bin, d.bin) %>% 
  summarise(S.mean = mean(survival_t1),
            d.mean = mean(d.stand),
            v.mean = mean(volume_t))
S.pred %>% 
  ggplot() +
  geom_line(aes(x = x.v, y = S.mean, colour = as.factor(d.mean)), size = 1) +
  scale_colour_manual(values = c("firebrick2", "darkorchid2", "dodgerblue1", "navy",
                                 "firebrick2", "darkorchid2", "dodgerblue1", "navy")) +
  geom_point(data = S.mean.df, aes(x = v.mean, y = S.mean, 
                                   colour = as.factor(d.bin)), size = 3) +
  labs(x = "Log-Volume (v)", y = "Survival") +
  theme_bw() -> Mod.S.avg.vPlot





##### Create LM for probability of recruitment from seed --------------------------------------------------

# Restore CData.s since we deleted stuff earlier, then standardise density
CData.s <- CData %>%
  mutate(d.stand = (weighted.dens - mean(weighted.dens, na.rm = TRUE)) / sd(weighted.dens, na.rm = TRUE))

# Collapse data to give unique combinations of year, site, transect, and window
# This helps avoid inflation of data
# Create simple linear model, using the PMod LM for the seed counts used in calculating the probabilities
CData.s %>% 
  group_by(site, transect, actual.window, year_t) %>% 
  select(recruit.prob, d.stand) %>% 
  summarise(recruit.prob = unique(recruit.prob),
            d.stand = unique(d.stand)) %>% 
  lm(recruit.prob ~ d.stand, data = .) -> Mod.P

#plot(CData.s$d.stand,CData.s$recruit.prob)
#abline(coef(Mod.P))
#anova(Mod.P)

# Create vector of coefficients for this model
Mod.P.cf <- c()
  
  # Intercept
  Mod.P.cf[1] <- coef(Mod.P)[1]
  
  # Volume coefficient
  Mod.P.cf[2] <- coef(Mod.P)[2]

# Less complicated overall compared to the previous models





##### Clean up --------------------------------------------------------------------------------------------

# Clean up variables from global environment
remove(d.and.v, d.and.v.2, d.cut, S.mean.df, S.pred, v.cut, d.range, v.range)


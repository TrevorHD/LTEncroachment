##### VFind probability of recruitment --------------------------------------------------------------------

# Deprecated because a better method of finding recruitment probabilities was implemented

# Create data frame of recruits; will be used for later calculations, but not here
CData.Recruits <- filter(CData, new.plant_t1 == 1 | seedling_t1 == 1, volume_t1 < 8)

# Calculate total amount of reproduction (empirical) in each 5-m window at time t
# Note that t is used because seeds at t will become seedlings at t1
for(i in 1:nrow(CData)){
  CData$seeds.emp[i] <- sum(CData$total.reproduction_t[CData$actual.window == CData$actual.window[i] &
                                                         CData$transect == CData$transect[i] &
                                                         CData$site == CData$site[i]], na.rm = T)}

# Create simple linear model for total reproduction in 5-m window as a function of weighted density  
seed.mod <- lm(seeds.emp ~ weighted.dens, data = CData)

# FMod: linear model fully used to estimate reproductive count for all windows
# PMod: linear model used to estimate reproductive count for windows with 0 reproduction
for(i in 1:nrow(CData)){
  CData$seeds.FMod[i] <- coef(seed.mod)[1] + coef(seed.mod)[2]*CData$weighted.dens[i]
  ifelse(CData$seeds.emp[i] == 0,
         CData$seeds.PMod[i] <- CData$seeds.FMod[i],
         CData$seeds.PMod[i] <- CData$seeds.emp[i])}

# Calculate total number of seedlings (recruits) in each 5-m window at t1
# Calculate recruitment rate per reproductive structure
for(i in 1:nrow(CData)){
  CData$recruits[i] <- sum(CData$new.plant_t1[CData$actual.window == CData$actual.window[i] &
                                                CData$transect == CData$transect[i] &
                                                CData$site == CData$site[i]], na.rm = T)
  CData$recruit.prob[i] <- CData$recruits[i]/CData$seeds.PMod[i]}





##### Flowering probability model averaging ---------------------------------------------------------------

# Deprecated because we decided to just pick the top model instead of using weighted model averages

# Keep all models in which the weight is greater than 0.001
Mod.F.AIC <- AICtab(Mod.F[[8]], Mod.F[[5]], Mod.F[[2]], Mod.F[[7]], Mod.F[[4]],
                    weights = TRUE, sort = FALSE)

# Create vector of average GLM coefficients from the kept models
Mod.F.avg.cf <- c()

# Intercept
Mod.F.avg.cf[1] <- 
  Mod.F.AIC$weight[1] * fixef(Mod.F[[8]])["(Intercept)"] + 
  Mod.F.AIC$weight[2] * fixef(Mod.F[[5]])["(Intercept)"] + 
  Mod.F.AIC$weight[3] * fixef(Mod.F[[2]])["(Intercept)"] + 
  Mod.F.AIC$weight[4] * fixef(Mod.F[[7]])["(Intercept)"] +
  Mod.F.AIC$weight[5] * fixef(Mod.F[[4]])["(Intercept)"]

# Volume coefficient
Mod.F.avg.cf[2] <- 
  Mod.F.AIC$weight[1] * fixef(Mod.F[[8]])["volume_t"] + 
  Mod.F.AIC$weight[2] * fixef(Mod.F[[5]])["volume_t"] + 
  Mod.F.AIC$weight[3] * fixef(Mod.F[[2]])["volume_t"] + 
  Mod.F.AIC$weight[4] * fixef(Mod.F[[7]])["volume_t"] +
  Mod.F.AIC$weight[5] * fixef(Mod.F[[4]])["volume_t"]

# Density coefficient
Mod.F.avg.cf[3] <- 
  Mod.F.AIC$weight[1] * fixef(Mod.F[[8]])["d.stand"] +
  Mod.F.AIC$weight[2] * fixef(Mod.F[[5]])["d.stand"] + 
  Mod.F.AIC$weight[4] * fixef(Mod.F[[7]])["d.stand"] +
  Mod.F.AIC$weight[5] * fixef(Mod.F[[4]])["d.stand"]

# Volume and density interaction coefficient
Mod.F.avg.cf[4] <- 
  Mod.F.AIC$weight[1] * fixef(Mod.F[[8]])["volume_t:d.stand"] +
  Mod.F.AIC$weight[2] * fixef(Mod.F[[5]])["volume_t:d.stand"]

# Density quadratic coefficient
Mod.F.avg.cf[5] <- 
  Mod.F.AIC$weight[1] * fixef(Mod.F[[8]])["I(d.stand^2)"] +
  Mod.F.AIC$weight[4] * fixef(Mod.F[[7]])["I(d.stand^2)"]

# Volume and quadratic density interaction coefficient
Mod.F.avg.cf[6] <-
  Mod.F.AIC$weight[1] * fixef(Mod.F[[8]])["volume_t:I(d.stand^2)"]





##### Growth model averaging ------------------------------------------------------------------------------

# Deprecated because we decided to just pick the top model instead of using weighted model averages

# Keep all models in which the weight is greater than 0.001
Mod.G.AIC <- AICtab(Mod.G[[5]], Mod.G[[4]], Mod.G[[8]], Mod.G[[7]],
                    weights = TRUE, sort = FALSE)

# Create vector of average GLM coefficients from the kept models
Mod.G.avg.cf <- c()

# Intercept
Mod.G.avg.cf[1] <- 
  Mod.G.AIC$weight[1] * fixef(Mod.G[[5]])["(Intercept)"] +
  Mod.G.AIC$weight[2] * fixef(Mod.G[[4]])["(Intercept)"] +
  Mod.G.AIC$weight[3] * fixef(Mod.G[[8]])["(Intercept)"] + 
  Mod.G.AIC$weight[4] * fixef(Mod.G[[7]])["(Intercept)"]

# Volume coefficient
Mod.G.avg.cf[2] <- 
  Mod.G.AIC$weight[1] * fixef(Mod.G[[5]])["volume_t"] +
  Mod.G.AIC$weight[2] * fixef(Mod.G[[4]])["volume_t"] +
  Mod.G.AIC$weight[3] * fixef(Mod.G[[8]])["volume_t"] + 
  Mod.G.AIC$weight[4] * fixef(Mod.G[[7]])["volume_t"]

# Density coefficient
Mod.G.avg.cf[3] <- 
  Mod.G.AIC$weight[1] * fixef(Mod.G[[5]])["d.stand"] +
  Mod.G.AIC$weight[2] * fixef(Mod.G[[4]])["d.stand"] +
  Mod.G.AIC$weight[3] * fixef(Mod.G[[8]])["d.stand"] +
  Mod.G.AIC$weight[4] * fixef(Mod.G[[7]])["d.stand"]

# Volume and density interaction coefficient
Mod.G.avg.cf[4] <- 
  Mod.G.AIC$weight[1] * fixef(Mod.G[[5]])["volume_t:d.stand"] +
  Mod.G.AIC$weight[3] * fixef(Mod.G[[8]])["volume_t:d.stand"]

# Density quadratic coefficient
Mod.G.avg.cf[5] <- 
  Mod.G.AIC$weight[3] * fixef(Mod.G[[8]])["I(d.stand^2)"] +
  Mod.G.AIC$weight[4] * fixef(Mod.G[[7]])["I(d.stand^2)"]

# Volume and quadratic density interaction coefficient
Mod.G.avg.cf[6] <- 
  Mod.G.AIC$weight[3] * fixef(Mod.G[[8]])["volume_t:I(d.stand^2)"]





##### Reproductive structure model averaging --------------------------------------------------------------

# Deprecated because we decided to just pick the top model instead of using weighted model averages

# Keep all models in which the weight is greater than 0.001
Mod.R.AIC <- AICtab(Mod.R[[7]], Mod.R[[8]], Mod.R[[5]], Mod.R[[4]], 
                    weights = TRUE, sort = FALSE)

# Create vector of average GLM coefficients from the kept models
Mod.R.avg.cf <- c()

# Intercept
Mod.R.avg.cf[1] <- 
  Mod.R.AIC$weight[1] * fixef(Mod.R[[7]])["(Intercept)"] + 
  Mod.R.AIC$weight[2] * fixef(Mod.R[[8]])["(Intercept)"] +
  Mod.R.AIC$weight[3] * fixef(Mod.R[[5]])["(Intercept)"] +
  Mod.R.AIC$weight[4] * fixef(Mod.R[[4]])["(Intercept)"]

# Volume coefficient
Mod.R.avg.cf[2] <- 
  Mod.R.AIC$weight[1] * fixef(Mod.R[[7]])["volume_t"] + 
  Mod.R.AIC$weight[2] * fixef(Mod.R[[8]])["volume_t"] + 
  Mod.R.AIC$weight[3] * fixef(Mod.R[[5]])["volume_t"] +
  Mod.R.AIC$weight[4] * fixef(Mod.R[[4]])["volume_t"]

# Density coefficient
Mod.R.avg.cf[3] <- 
  Mod.R.AIC$weight[1] * fixef(Mod.R[[7]])["d.stand"] + 
  Mod.R.AIC$weight[2] * fixef(Mod.R[[8]])["d.stand"] + 
  Mod.R.AIC$weight[3] * fixef(Mod.R[[5]])["d.stand"] +
  Mod.R.AIC$weight[4] * fixef(Mod.R[[4]])["d.stand"]

# Volume and density interaction coefficient
Mod.R.avg.cf[4] <- 
  Mod.R.AIC$weight[2] * fixef(Mod.R[[8]])["volume_t:d.stand"] +
  Mod.R.AIC$weight[3] * fixef(Mod.R[[5]])["volume_t:d.stand"]

# Density quadratic coefficient
Mod.R.avg.cf[5] <- 
  Mod.R.AIC$weight[1] * fixef(Mod.R[[7]])["I(d.stand^2)"] +
  Mod.R.AIC$weight[2] * fixef(Mod.R[[8]])["I(d.stand^2)"]

# Volume and quadratic density interaction coefficient
Mod.R.avg.cf[6] <- 
  Mod.R.AIC$weight[2] * fixef(Mod.R[[8]])["volume_t:I(d.stand^2)"]





##### Survival model averaging and visualisation ----------------------------------------------------------

# Deprecated because we decided to just pick the top model instead of using weighted model averages

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





##### Visualising survival data (NS) ----------------------------------------------------------------------

# Deprecated because visualisation is not useful and is not part of the main code

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

# Clean up variables from global environment
remove(d.and.v, d.and.v.2, d.cut, S.mean.df, S.pred, v.cut, d.range, v.range)





##### Visualising survival data (BS) ----------------------------------------------------------------------

# Deprecated because visualisation is not useful and is not part of the main code

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

# Clean up variables from global environment
remove(d.and.v, d.and.v.2, d.cut, S.mean.df, S.pred, v.cut, d.range, v.range, i, DATA)





##### Create GLM for flowering probability ----------------------------------------------------------------

# Define inverse logit function used as link to binomial distribution
invlogit <- function(x){exp(x)/(1 + exp(x))} 

# Remove entries for which there is no recorded volume
# Then standardise density and introduce unique transect identifier
boot.CData.s <- boot.CData %>% 
  drop_na(volume_t) %>% 
  mutate(d.stand = (weighted.dens - mean(weighted.dens, na.rm = TRUE)) / sd(weighted.dens, na.rm = TRUE),
         unique.transect = interaction(transect, site))

# Create a list of possible flowering models
Mod.F <- list()

# "Null" model; only random effects of site and transect within site
Mod.F[[1]] <- glmer(did.flower ~ (1 | unique.transect),
                    data = boot.CData.s, family = "binomial")

# Size-only model
Mod.F[[2]] <- glmer(did.flower ~ volume_t + (1 | unique.transect),
                    data = boot.CData.s, family = "binomial")

# Density-only model
Mod.F[[3]] <- glmer(did.flower ~ d.stand + (1 | unique.transect),
                    data = boot.CData.s, family = "binomial")

# Size and density (additive)
Mod.F[[4]] <- glmer(did.flower ~ volume_t + d.stand + (1 | unique.transect),
                    data = boot.CData.s, family = "binomial")

# Size and density (interactive)
Mod.F[[5]] <- glmer(did.flower ~ volume_t * d.stand + (1 | unique.transect),
                    data = boot.CData.s, family = "binomial")

# Density-only model (quadratic)
Mod.F[[6]] <- glmer(did.flower ~ d.stand + I(d.stand^2) + (1 | unique.transect),
                    data = boot.CData.s, family = "binomial")

# Size (linear) and density (quadratic)
Mod.F[[7]] <- glmer(did.flower ~ volume_t + d.stand + I(d.stand^2) + (1 | unique.transect),
                    data = boot.CData.s, family = "binomial")

# Size and density (interactive, quadratic)
Mod.F[[8]] <- glmer(did.flower ~ volume_t * d.stand + volume_t * I(d.stand^2) + (1 | unique.transect),
                    data = boot.CData.s, family = "binomial")

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

# Restore boot.CData.s since we deleted stuff earlier
# Remove entries for which there is no recorded volume and log growth is NA
# Then standardise density and introduce unique transect identifier
boot.CData.s <- boot.CData %>% 
  drop_na(volume_t, logGR) %>% 
  mutate(d.stand = (weighted.dens - mean(weighted.dens, na.rm = TRUE)) / sd(weighted.dens, na.rm = TRUE),
         unique.transect = interaction(transect, site))

# Create a list of possible growth models
Mod.G <- list()

# "Null" model; only random effects of site and transect within site
Mod.G[[1]] <- lmer(logGR ~ (1 | unique.transect),
                   data = boot.CData.s, REML = FALSE)

# Size-only model
Mod.G[[2]] <- lmer(logGR ~ volume_t + (1 | unique.transect),
                   data = boot.CData.s, REML = FALSE)

# Density-only model
Mod.G[[3]] <- lmer(logGR ~ d.stand + (1 | unique.transect),
                   data = boot.CData.s, REML = FALSE)

# Size and density (additive)
Mod.G[[4]] <- lmer(logGR ~ volume_t + d.stand + (1 | unique.transect),
                   data = boot.CData.s, REML = FALSE)

# Size and density (interactive)
Mod.G[[5]] <- lmer(logGR ~ volume_t * d.stand + (1 | unique.transect),
                   data = boot.CData.s, REML = FALSE)

# Density-only model (quadratic)
Mod.G[[6]] <- lmer(logGR ~ d.stand + I(d.stand^2) + (1 | unique.transect),
                   data = boot.CData.s, REML = FALSE)

# Size (linear) and density (quadratic)
Mod.G[[7]] <- lmer(logGR ~ volume_t + d.stand + I(d.stand^2) + (1 | unique.transect),
                   data = boot.CData.s, REML = FALSE)

# Size and density (interactive, quadratic)
Mod.G[[8]] <- lmer(logGR ~ volume_t * d.stand + volume_t * I(d.stand^2) + (1 | unique.transect),
                   data = boot.CData.s, REML = FALSE)

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

# Restore boot.CData.s since we deleted stuff earlier
# Remove entries for which there is no recorded volume or that did not flower
# Then standardise density and introduce unique transect identifier
boot.CData.s <- boot.CData %>%
  drop_na(volume_t) %>%
  filter(did.flower == 1) %>% 
  mutate(d.stand = (weighted.dens - mean(weighted.dens, na.rm = TRUE)) / sd(weighted.dens, na.rm = TRUE),
         unique.transect = interaction(transect, site))

# Create a list of possible reproduction models
Mod.R <- list()

# "Null" model; only random effects of site and transect within site
Mod.R[[1]] <- lmer(logTR1 ~ (1 | unique.transect),
                   data = boot.CData.s, REML = FALSE)

# Size-only model
Mod.R[[2]] <- lmer(logTR1 ~ volume_t + (1 | unique.transect),
                   data = boot.CData.s, REML = FALSE)

# Density-only model
Mod.R[[3]] <- lmer(logTR1 ~ d.stand + (1 | unique.transect),
                   data = boot.CData.s, REML = FALSE)

# Size and density (additive)
Mod.R[[4]] <- lmer(logTR1 ~ volume_t + d.stand + (1 | unique.transect),
                   data = boot.CData.s, REML = FALSE)

# Size and density (interactive)
Mod.R[[5]] <- lmer(logTR1 ~ volume_t * d.stand + (1 | unique.transect),
                   data = boot.CData.s, REML = FALSE)

# Density-only model (quadratic)
Mod.R[[6]] <- lmer(logTR1 ~ d.stand + I(d.stand^2) + (1 | unique.transect),
                   data = boot.CData.s, REML = FALSE)

# Size (linear) and density (quadratic)
Mod.R[[7]] <- lmer(logTR1 ~ volume_t + d.stand + I(d.stand^2) + (1 | unique.transect),
                   data = boot.CData.s, REML = FALSE)

# Size and density (interactive, quadratic)
Mod.R[[8]] <- lmer(logTR1 ~ volume_t * d.stand + volume_t * I(d.stand^2) + (1 | unique.transect),
                   data = boot.CData.s, REML = FALSE)

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

# Combine transplants with large shrubs for survival analysis
# Keep only location info, survival, volume, and density
select(boot.CData.Transplants, "site", "transect", "actual.window", 
       "spring_survival_t1", "volume_t", "weighted.dens", "transplant") %>% 
  rename("survival_t1" = "spring_survival_t1") %>% 
  rbind(select(boot.CData, "site", "transect", "actual.window", 
               "survival_t1", "volume_t", "weighted.dens", "transplant")) -> boot.CData.AllSurvival

# Use models for transplants (small individuals) only
# Survival for larger individuals is pretty much guaranteed

# Remove entries for which there is no recorded volume or survival
# Standardise density; introduce unique transect identifier
boot.CData.AllSurvival.s <- boot.CData.AllSurvival %>% 
  drop_na(volume_t, survival_t1) %>% 
  mutate(d.stand = (weighted.dens - mean(weighted.dens, na.rm = TRUE)) / sd(weighted.dens, na.rm = TRUE),
         unique.transect = interaction(transect, site),
         survival_t1 = as.numeric(survival_t1))

## how much mortality was there just among the naturally occurring plants?
table(boot.CData.AllSurvival.s$survival_t1[boot.CData.AllSurvival.s$transplant==F])
## why is there no natural mortality? there are natural deaths in the original data frame
table(CData.Demography$survival_t1)

# Create a list of possible survival models
Mod.S <- list()

# "Null" model; only random effects of site and transect within site
Mod.S[[1]] <- glmer(survival_t1 ~ (1 | unique.transect),
                    data = subset(boot.CData.AllSurvival.s, transplant == TRUE), family = "binomial")

# Size-only model
Mod.S[[2]] <- glmer(survival_t1 ~ volume_t + (1 | unique.transect),
                    data = subset(boot.CData.AllSurvival.s, transplant == TRUE), family = "binomial")

# Density-only model
Mod.S[[3]] <- glmer(survival_t1 ~ d.stand + (1 | unique.transect),
                    data = subset(boot.CData.AllSurvival.s, transplant == TRUE), family = "binomial")

# Size and density (additive)
Mod.S[[4]] <- glmer(survival_t1 ~ volume_t + d.stand + (1 | unique.transect),
                    data = subset(boot.CData.AllSurvival.s, transplant == TRUE), family = "binomial")

# Size and density (interactive)
Mod.S[[5]] <- glmer(survival_t1 ~ volume_t * d.stand + (1 | unique.transect),
                    data = subset(boot.CData.AllSurvival.s, transplant == TRUE), family = "binomial")

# Density-only model (quadratic)
Mod.S[[6]] <- glmer(survival_t1 ~ d.stand + I(d.stand^2) + (1 | unique.transect),
                    data = subset(boot.CData.AllSurvival.s, transplant == TRUE), family = "binomial")

# Size (linear) and density (quadratic)
Mod.S[[7]] <- glmer(survival_t1 ~ volume_t + d.stand + I(d.stand^2) + (1 | unique.transect),
                    data = subset(boot.CData.AllSurvival.s, transplant == TRUE), family = "binomial")

# Size and density (interactive, quadratic)
Mod.S[[8]] <- glmer(survival_t1 ~ volume_t * d.stand + volume_t * I(d.stand^2) + (1 | unique.transect),
                    data = subset(boot.CData.AllSurvival.s, transplant == TRUE), family = "binomial")

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
boot.CData.Recruitment <- merge(CData.Transects, Windows, 
                                by.x = c("site", "transect", "window"),
                                by.y = c("site", "transect", "window"))

# Use the reproduction model to calculate number of reproductive structures in a single year for each plant
mutate(boot.CData.Recruitment,
       d.stand = (weighted.dens - mean(weighted.dens, na.rm = TRUE)) / sd(weighted.dens, na.rm = TRUE),
       seeds = exp(Mod.R.top.cf[1] + Mod.R.top.cf[2]*log(volume) + Mod.R.top.cf[3]*d.stand + 
                     Mod.R.top.cf[5]*(d.stand^2))) -> boot.CData.Recruitment

# Calculate number of seeds in a single year for each 5-m window
for(i in 1:nrow(boot.CData.Recruitment)){
  boot.CData.Recruitment$seeds.win[i] <- sum(boot.CData.Recruitment$seeds[boot.CData.Recruitment$window == boot.CData.Recruitment$window[i] &
                                                                            boot.CData.Recruitment$transect == boot.CData.Recruitment$transect[i] &
                                                                            boot.CData.Recruitment$site == boot.CData.Recruitment$site[i]], na.rm = T)}

# Select only one instance of each unique combination of site, transect, and window, then merge with CData
# We're doing this because we're interested in total seeds in each window, not seeds per plant in each window
distinct(select(boot.CData.Recruitment, "site", "transect", "window", "seeds.win")) %>% 
  merge(boot.CData, ., by.x = c("site", "transect", "actual.window"),
        by.y = c("site", "transect", "window")) -> boot.CData.Recruitment

# Calculate recruitment rates for each 5-m window over 1- and 4-year periods
boot.CData.Recruitment <- mutate(boot.CData.Recruitment, recruit.prob.1y = recruits.1y/(seeds.win),
                                 recruit.prob.4y = recruits.4y/(4*seeds.win))

# Select only one instance of each unique combination of site, transect, window, and year
# Duplicate instances of these unique combinations will inflate the data, which we don't want
boot.CData.Recruitment %>% 
  mutate(d.stand = (weighted.dens - mean(weighted.dens, na.rm = TRUE)) / sd(weighted.dens, na.rm = TRUE)) %>% 
  group_by(site, transect, actual.window, year_t1) %>% 
  select(d.stand, seeds.win, recruits.1y, recruits.4y, recruit.prob.1y, recruit.prob.4y) %>% 
  summarise(d.stand = unique(d.stand),
            seeds.win = round(unique(seeds.win), digits = 0),
            recruits.1y = unique(recruits.1y),
            recruits.4y = unique(recruits.4y),
            recruit.prob.1y = unique(recruit.prob.1y),
            recruit.prob.4y = unique(recruit.prob.4y)) %>% 
  mutate(unique.transect = interaction(transect, site)) -> boot.CData.Recruitment

# Create a list of possible per-seed recruitment models using 1-year rates
Mod.P1 <- list()

# "Null" model; only random effects of site and transect within site
Mod.P1[[1]] <- glmer(cbind(recruits.1y, seeds.win - recruits.1y) ~ (1 | unique.transect), 
                     data = boot.CData.Recruitment, family = "binomial")

# Density-only model
Mod.P1[[2]] <- glmer(cbind(recruits.1y, seeds.win - recruits.1y) ~ d.stand  + (1 | unique.transect), 
                     data = boot.CData.Recruitment, family = "binomial")

# Density-only model (quadratic)
Mod.P1[[3]] <- glmer(cbind(recruits.1y, seeds.win - recruits.1y) ~ d.stand + I(d.stand^2) + (1 | unique.transect), 
                     data = boot.CData.Recruitment, family = "binomial")

# Create a list of possible per-seed recruitment models using 4-year rates
Mod.P4 <- list()

# "Null" model; only random effects of site and transect within site
Mod.P4[[1]] <- glmer(cbind(recruits.4y, seeds.win - recruits.4y) ~ (1 | unique.transect), 
                     data = boot.CData.Recruitment, family = "binomial")

# Density-only model
Mod.P4[[2]] <- glmer(cbind(recruits.4y, seeds.win - recruits.4y) ~ d.stand  + (1 | unique.transect), 
                     data = boot.CData.Recruitment, family = "binomial")

# Density-only model (quadratic)
Mod.P4[[3]] <- glmer(cbind(recruits.4y, seeds.win - recruits.4y) ~ d.stand + I(d.stand^2) + (1 | unique.transect), 
                     data = boot.CData.Recruitment, family = "binomial")





##### Run scripts that do not change between the two scenarios --------------------------------------------

# This will take several minutes; be patient

# "01_SeedVelocities"
# Calculate terminal velocities of seeds and fit to distributions
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/01_SeedVelocities.R")

# "02_WindSpeeds"
# Load in wind speeds and fit to distributions
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/02_WindSpeeds.R")

# "03_Dispersal"
# Construct dispersal kernels for seeds
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/03_Dispersal.R")

# "04_CDataPrep"
# Tidy up demography data before creating demography models
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/04_CDataPrep.R")





##### Wavespeeds and population growth for normal survival scenario ---------------------------------------

# This will take several minutes; be patient

# "05_CDataAnalysis_NS"
# Create demography models for growth, reproduction, survival, etc. under normal circumstances
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/05_CDataAnalysis_NS.R")

# "06_SIPM"
# Spatial integral projection model that calculates wavespeeds
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/06_SIPM.R")

# Wavespeeds as function of s; growth as function of density
c.values <- Wavespeed(200)
lambda <- c()
for(i in seq(-1.3, max(CData.s$d.stand), length.out = 100)){
  lambda.i <- TransMatrix(n = 200, d = i)
  lambda <- append(lambda, Re(eigen(lambda.i$matrix)$values[1]))}

# Clean up
remove(lambda.i, TM, i)





##### Wavespeeds and population growth for higher survival scenario ---------------------------------------

# This will take several minutes; be patient

# "05_CDataAnalysis_BS"
# Replace survival model with one for year with higher survival from above-average rainfall
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/05_CDataAnalysis_BS.R")

# "06_SIPM"
# Spatial integral projection model that calculates wavespeeds
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/06_SIPM.R")

# Wavespeeds as function of s; growth as function of density
c.values.2 <- Wavespeed(200)
lambda.2 <- c()
for(i in seq(-1.3, max(CData.s$d.stand), length.out = 100)){
  lambda.i <- TransMatrix(n = 200, d = i)
  lambda.2 <- append(lambda.2, Re(eigen(lambda.i$matrix)$values[1]))}

# Clean up
remove(lambda.i, TM, i)





# ----- Flowering as a function of volume and density (1 of 4) --------------------------------------------

# Remove entries for which there is no recorded volume
# Then standardise density and introduce unique transect identifier
CData.s <- CData %>% 
  drop_na(volume_t) %>% 
  mutate(d.stand = (weighted.dens - mean(weighted.dens, na.rm = TRUE)) / sd(weighted.dens, na.rm = TRUE),
         unique.transect = interaction(transect, site))

# Calculate range of data
v.range <- range(CData.s$volume_t)
d.range <- range(CData.s$d.stand)

# Calculate predictions for 4 "bins" of volume categories
v.cut <- CData.s %>% 
  mutate(v.bin = cut_number(volume_t, n = 4)) %>% 
  group_by(v.bin) %>% 
  summarise(v.mean = mean(volume_t))
d.and.v <- as.tibble(expand.grid(v.cut$v.mean, seq(d.range[1], d.range[2], 0.5)))
names(d.and.v) <- c("v.mean", "x.d") 
F.pred <- d.and.v %>% 
  mutate(F.mean = invlogit(Mod.F.top.cf[1] + 
                             Mod.F.top.cf[2] * v.mean + 
                             Mod.F.top.cf[3] * x.d +
                             Mod.F.top.cf[4] * v.mean * x.d +
                             Mod.F.top.cf[5] * (x.d^2)+
                             Mod.F.top.cf[6] * v.mean * (x.d^2)))

# Get binned means of the data with respect to volume and density
F.mean.df <- CData.s %>% 
  mutate(v.bin = cut_number(volume_t, n = 4),
         d.bin = cut_number(d.stand, n = 6)) %>% 
  group_by(d.bin, v.bin) %>% 
  summarise(F.mean = mean(did.flower),
            d.mean = mean(d.stand),
            v.mean = mean(volume_t))

# Generate visualisation; flowering probability as function of density, with 4 volume bins
F.pred %>% 
  ggplot() +
  geom_line(aes(x = x.d, y = F.mean, colour = as.factor(v.mean)), size = 3) +
  scale_colour_manual(values = c("darkorchid2", "dodgerblue1", "navy", "firebrick2",
                                 "darkorchid2", "dodgerblue1", "navy", "firebrick2")) +
  geom_point(data = F.mean.df, aes(x = d.mean, y = F.mean, 
                                   colour = as.factor(v.bin)), size = 7) +
  labs(x = NULL, y = "Flowering Probability") +
  ylim(0, 1) +
  xlim(-1.5, 3) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        plot.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 38, margin = margin(t = 4, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 38, margin = margin(t = 0, r = 4, b = 0, l = 0)),
        axis.title.x = element_text(size = 43, margin = margin(t = 16, r = 0, b = 0, l = 6)),
        axis.title.y = element_text(size = 43, margin = margin(t = 80, r = 16, b = 0, l = 0)),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(8, "points")) -> Mod.F.top.dPlot

# Do same process as previous 3 steps, but switch density and volume
# This produces a graph of flowering probability as a function of volume, with 4 density bins
d.cut <- CData.s %>% 
  mutate(d.bin = cut_number(d.stand, n = 4)) %>% 
  group_by(d.bin) %>% 
  summarise(d.mean = mean(d.stand))
d.and.v.2 <- as.tibble(expand.grid(d.cut$d.mean, seq(v.range[1], v.range[2], 0.5)))
names(d.and.v.2) <- c("d.mean", "x.v") 
F.pred <- d.and.v.2 %>% 
  mutate(F.mean = invlogit(Mod.F.top.cf[1] + 
                             Mod.F.top.cf[2] * x.v + 
                             Mod.F.top.cf[3] * d.mean +
                             Mod.F.top.cf[4] * x.v * d.mean +
                             Mod.F.top.cf[5] * (d.mean^2) +
                             Mod.F.top.cf[6] * x.v * (d.mean^2)))
F.mean.df <- CData.s %>% 
  mutate(d.bin = cut_number(d.stand, n = 4),
         v.bin = cut_number(volume_t, n = 6)) %>% 
  group_by(v.bin, d.bin) %>% 
  summarise(F.mean = mean(did.flower),
            d.mean = mean(d.stand),
            v.mean = mean(volume_t))
F.pred %>% 
  ggplot() +
  geom_line(aes(x = x.v, y = F.mean, colour = as.factor(d.mean)), size = 3) +
  scale_colour_manual(values = c("darkorchid2", "firebrick2", "dodgerblue1", "darkorchid2",
                                 "navy", "firebrick2", "dodgerblue1", "navy")) +
  scale_y_continuous(position = "right") +
  geom_point(data = F.mean.df, aes(x = v.mean, y = F.mean, 
                                   colour = as.factor(d.bin)), size = 7) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        plot.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 38, margin = margin(t = 4, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 38, margin = margin(t = 0, r = 4, b = 0, l = 0)),
        axis.title.x = element_text(size = 43, margin = margin(t = 16, r = 0, b = 0, l = 6)),
        axis.title.y = element_text(size = 43, margin = margin(t = 80, r = 16, b = 0, l = 0)),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(8, "points")) -> Mod.F.top.vPlot





# ----- Growth as a function of volume and density (2 of 4) -----------------------------------------------

# Restore CData.s since we deleted stuff earlier
# Remove entries for which there is no recorded volume and log growth is NA
# Then standardise density and introduce unique transect identifier
CData.s <- CData %>% 
  drop_na(volume_t, logGR) %>% 
  mutate(d.stand = (weighted.dens - mean(weighted.dens, na.rm = TRUE)) / sd(weighted.dens, na.rm = TRUE),
         unique.transect = interaction(transect, site))

# Calculate range of data
v.range <- range(CData.s$volume_t)
d.range <- range(CData.s$d.stand)

# Calculate predictions for 4 "bins" of volume categories
v.cut <- CData.s %>% 
  mutate(v.bin = cut_number(volume_t, n = 4)) %>% 
  group_by(v.bin) %>% 
  summarise(v.mean = mean(volume_t))
d.and.v <- as.tibble(expand.grid(v.cut$v.mean, seq(d.range[1], d.range[2], 0.5)))
names(d.and.v) <- c("v.mean", "x.d") 
G.pred <- d.and.v %>% 
  mutate(G.mean = Mod.G.top.cf[1] + 
           Mod.G.top.cf[2] * v.mean + 
           Mod.G.top.cf[3] * x.d +
           Mod.G.top.cf[4] * v.mean * x.d +
           Mod.G.top.cf[5] * (x.d^2)+
           Mod.G.top.cf[6] * v.mean * (x.d^2))

# Get binned means of the data with respect to volume and density
G.mean.df <- CData.s %>% 
  mutate(v.bin = cut_number(volume_t, n = 4),
         d.bin = cut_number(d.stand, n = 6)) %>% 
  group_by(d.bin, v.bin) %>% 
  summarise(G.mean = mean(logGR),
            d.mean = mean(d.stand),
            v.mean = mean(volume_t))

# Generate visualisation; growth ratio as function of density, with 4 volume bins
G.pred %>% 
  ggplot() +
  geom_line(aes(x = x.d, y = G.mean, colour = as.factor(v.mean)), size = 3) +
  scale_colour_manual(values = c("darkorchid2", "dodgerblue1", "navy", "firebrick2",
                                 "darkorchid2", "dodgerblue1", "navy", "firebrick2")) +
  scale_y_continuous(breaks = seq(0, 1, length.out = 5), limits = c(-0.1, 1)) +
  geom_point(data = G.mean.df, aes(x = d.mean, y = G.mean, 
                                   colour = as.factor(v.bin)), size = 7) +
  labs(x = NULL, y = "Log Annual Growth Ratio") +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        plot.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 38, margin = margin(t = 4, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 38, margin = margin(t = 0, r = 4, b = 0, l = 0)),
        axis.title.x = element_text(size = 43, margin = margin(t = 16, r = 0, b = 0, l = 6)),
        axis.title.y = element_text(size = 43, margin = margin(t = 80, r = 16, b = 0, l = 0)),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(8, "points")) -> Mod.G.top.dPlot

# Do same process as previous 3 steps, but switch density and volume
# This produces a graph of annual growth ratio as a function of volume, with 4 density bins
d.cut <- CData.s %>% 
  mutate(d.bin = cut_number(d.stand, n = 4)) %>% 
  group_by(d.bin) %>% 
  summarise(d.mean = mean(d.stand))
d.and.v.2 <- as.tibble(expand.grid(d.cut$d.mean, seq(v.range[1], v.range[2], 0.5)))
names(d.and.v.2) <- c("d.mean", "x.v") 
G.pred <- d.and.v.2 %>% 
  mutate(G.mean = Mod.G.top.cf[1] + 
           Mod.G.top.cf[2] * x.v + 
           Mod.G.top.cf[3] * d.mean +
           Mod.G.top.cf[4] * x.v * d.mean +
           Mod.G.top.cf[5] * (d.mean^2) +
           Mod.G.top.cf[6] * x.v * (d.mean^2))
G.mean.df <- CData.s %>% 
  mutate(d.bin = cut_number(d.stand, n = 4),
         v.bin = cut_number(volume_t, n = 6)) %>% 
  group_by(v.bin, d.bin) %>% 
  summarise(G.mean = mean(logGR),
            d.mean = mean(d.stand),
            v.mean = mean(volume_t))
G.pred %>% 
  ggplot() +
  geom_line(aes(x = x.v, y = G.mean, colour = as.factor(d.mean)), size = 3) +
  scale_colour_manual(values = c("darkorchid2", "firebrick2", "dodgerblue1", "darkorchid2",
                                 "navy", "firebrick2", "dodgerblue1", "navy")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01), position = "right",
                     breaks = seq(0, 1, length.out = 5), limits = c(-0.05, 1)) +
  geom_point(data = G.mean.df, aes(x = v.mean, y = G.mean, 
                                   colour = as.factor(d.bin)), size = 7) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        plot.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 38, margin = margin(t = 4, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 38, margin = margin(t = 0, r = 4, b = 0, l = 0)),
        axis.title.x = element_text(size = 43, margin = margin(t = 16, r = 0, b = 0, l = 6)),
        axis.title.y = element_text(size = 43, margin = margin(t = 80, r = 16, b = 0, l = 0)),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(8, "points")) -> Mod.G.top.vPlot





# ----- Reproduction as a function of volume and density (3 of 4) -----------------------------------------

# Restore CData.s since we deleted stuff earlier
# Remove entries for which there is no recorded volume or that did not flower
# Then standardise density and introduce unique transect identifier
CData.s <- CData %>%
  drop_na(volume_t) %>%
  filter(did.flower == 1) %>% 
  mutate(d.stand = (weighted.dens - mean(weighted.dens, na.rm = TRUE)) / sd(weighted.dens, na.rm = TRUE),
         unique.transect = interaction(transect, site))

# Calculate range of data
v.range <- range(CData.s$volume_t)
d.range <- range(CData.s$d.stand)

# Calculate predictions for 4 "bins" of volume categories
v.cut <- CData.s %>% 
  mutate(v.bin = cut_number(volume_t, n = 4)) %>% 
  group_by(v.bin) %>% 
  summarise(v.mean = mean(volume_t))
d.and.v <- as.tibble(expand.grid(v.cut$v.mean, seq(d.range[1], d.range[2], 0.5)))
names(d.and.v) <- c("v.mean", "x.d") 
R.pred <- d.and.v %>% 
  mutate(R.mean = Mod.R.top.cf[1] + 
           Mod.R.top.cf[2] * v.mean + 
           Mod.R.top.cf[3] * x.d +
           Mod.R.top.cf[4] * v.mean * x.d +
           Mod.R.top.cf[5] * (x.d^2) +
           Mod.R.top.cf[6] * v.mean * (x.d^2))

# Get binned means of the data with respect to volume and density
R.mean.df <- CData.s %>% 
  mutate(v.bin = cut_number(volume_t, n = 4),
         d.bin = cut_number(d.stand, n = 6)) %>% 
  group_by(d.bin, v.bin) %>% 
  summarise(R.mean = mean(logTR1),
            d.mean = mean(d.stand),
            v.mean = mean(volume_t))

# Generate visualisation; reproductive count as function of density, with 4 volume bins
R.pred %>% 
  ggplot() +
  geom_line(aes(x = x.d, y = R.mean, colour = as.factor(v.mean)), size = 3) +
  scale_colour_manual(values = c("darkorchid2", "dodgerblue1", "navy", "firebrick2",
                                 "firebrick2", "darkorchid2", "dodgerblue1", "navy")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  geom_point(data = R.mean.df, aes(x = d.mean, y = R.mean, 
                                   colour = as.factor(v.bin)), size = 7) +
  labs(x = "Standardised Weighted Density", y = "Log Reproductive Structure Count") +
  xlim(-1.5, 3) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        plot.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 38, margin = margin(t = 4, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 38, margin = margin(t = 0, r = 4, b = 0, l = 0)),
        axis.title.x = element_text(size = 43, margin = margin(t = 16, r = 0, b = 0, l = 6)),
        axis.title.y = element_text(size = 43, margin = margin(t = 80, r = 16, b = 0, l = 0)),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(8, "points")) -> Mod.R.top.dPlot

# Do same process as previous 3 steps, but switch density and volume
# This produces a graph of reproductive count as a function of volume, with 4 density bins
d.cut <- CData.s %>% 
  mutate(d.bin = cut_number(d.stand, n = 4)) %>% 
  group_by(d.bin) %>% 
  summarise(d.mean = mean(d.stand))
d.and.v.2 <- as.tibble(expand.grid(d.cut$d.mean, seq(v.range[1], v.range[2], 0.5)))
names(d.and.v.2) <- c("d.mean", "x.v") 
R.pred <- d.and.v.2 %>% 
  mutate(R.mean = Mod.R.top.cf[1] + 
           Mod.R.top.cf[2] * x.v + 
           Mod.R.top.cf[3] * d.mean +
           Mod.R.top.cf[4] * x.v * d.mean +
           Mod.R.top.cf[5] * (d.mean^2) +
           Mod.R.top.cf[6] * x.v * (d.mean^2))
R.mean.df <- CData.s %>% 
  mutate(d.bin = cut_number(d.stand, n = 4),
         v.bin = cut_number(volume_t, n = 6)) %>% 
  group_by(v.bin, d.bin) %>% 
  summarise(R.mean = mean(logTR1),
            d.mean = mean(d.stand),
            v.mean = mean(volume_t))
R.pred %>% 
  ggplot() +
  geom_line(aes(x = x.v, y = R.mean, colour = as.factor(d.mean)), size = 3) +
  scale_colour_manual(values = c("darkorchid2", "firebrick2", "dodgerblue1", "darkorchid2",
                                 "navy", "firebrick2", "dodgerblue1", "navy")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01), position = "right",
                     breaks = seq(-3, 9, length.out = 5), limits = c(-3, 9)) +
  scale_x_continuous(breaks = c(5, 8, 12, 15)) +
  geom_point(data = R.mean.df, aes(x = v.mean, y = R.mean, 
                                   colour = as.factor(d.bin)), size = 7) +
  labs(x = "Log-Volume", y = NULL) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        plot.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 38, margin = margin(t = 4, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 38, margin = margin(t = 0, r = 4, b = 0, l = 0)),
        axis.title.x = element_text(size = 43, margin = margin(t = 16, r = 0, b = 0, l = 6)),
        axis.title.y = element_text(size = 43, margin = margin(t = 80, r = 16, b = 0, l = 0)),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(8, "points")) -> Mod.R.top.vPlot





# ----- Combine figures from the 3 sections above (4 of 4) ------------------------------------------------

# Prepare graphics device
jpeg(filename = "Figure 3.jpeg", width = 2600, height = 2950, units = "px")

# Create blank page
grid.newpage()
plot.new()

# Set grid layout and activate it
gly <- grid.layout(2950, 2600)
pushViewport(viewport(layout = gly))

# Place graphs in a grid
print(Mod.F.top.dPlot, vp = viewport(layout.pos.row = 1:950, layout.pos.col = 1:1300))
print(Mod.F.top.vPlot, vp = viewport(layout.pos.row = 1:950, layout.pos.col = 1300:2600))
print(Mod.G.top.dPlot, vp = viewport(layout.pos.row = 1000:1950, layout.pos.col = 1:1300))
print(Mod.G.top.vPlot, vp = viewport(layout.pos.row = 1000:1950, layout.pos.col = 1300:2600))
print(Mod.R.top.dPlot, vp = viewport(layout.pos.row = 2000:2950, layout.pos.col = 1:1300))
print(Mod.R.top.vPlot, vp = viewport(layout.pos.row = 2000:2950, layout.pos.col = 1300:2600))

# Create legend
grid.rect(vp = viewport(layout.pos.row = 190:260, layout.pos.col = 1400:1470), gp = gpar(fill = "navy"))
grid.rect(vp = viewport(layout.pos.row = 290:360, layout.pos.col = 1400:1470), gp = gpar(fill = "dodgerblue1"))
grid.rect(vp = viewport(layout.pos.row = 390:460, layout.pos.col = 1400:1470), gp = gpar(fill = "darkorchid2"))
grid.rect(vp = viewport(layout.pos.row = 490:560, layout.pos.col = 1400:1470), gp = gpar(fill = "firebrick2"))
grid.text(label = c("[75, 100]", "[50, 75)", "[25, 50)", "[0, 25)"),
          x = rep(0.580, 4), y = c(0.925, 0.892, 0.858, 0.823), just = "left", gp = gpar(fontsize = 37))
grid.text(label = c("Volume percentile (left graphs)", "Density percentile (right graphs)"),
          x = rep(0.537, 2), y = c(0.970, 0.952), just = "left", gp = gpar(fontsize = 40))

# Deactivate grid layout; finalise graphics save
popViewport()
dev.off()

# Clean up mess from global environment
remove(d.and.v, d.and.v.2, d.cut, F.mean.df, F.pred, G.mean.df, G.pred, R.mean.df, R.pred,
       v.cut, d.range, v.range, Mod.F.top.dPlot, Mod.F.top.vPlot, Mod.G.top.dPlot, 
       Mod.G.top.vPlot, Mod.R.top.dPlot, Mod.R.top.vPlot, gly)





##### Wavespeeds and population growth for higher survival scenario ---------------------------------------

# NOTE: this part will not be used for now, until we get the other scenario working first

# This will take several minutes; be patient

# Begin bootstrapping
time.start <- Sys.time()
for(i in 1:boot.num){
  
  # "00_BootRes"
  # Run resampling subroutine for wind speeds, terminal velocities, and demography
  source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/06_BootRes.R")
  
  # "05_CDataAnalysis_NS"
  # Create demography models for growth, reproduction, survival, etc. under normal circumstances
  # Must run this before 05_CDataAnalysis_BS since it contains all of the demography models
  source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/05_CDataAnalysis_NS.R")
  
  # "05_CDataAnalysis_BS"
  # Replace survival model in 05_CDataAnalysis_NS with higher survival from above-average rainfall
  source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/05_CDataAnalysis_BS.R")
  
  # "06_SIPM"
  # Spatial integral projection model that calculates wavespeeds
  source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/07_SIPM.R")
  
  # Wavespeeds as function of s; growth as function of density
  c.values.2 <- Wavespeed(200)
  lambda.2 <- c()
  for(i in seq(-1.3, max(boot.CData.s$d.stand), length.out = 100)){
    lambda.i <- TransMatrix(n = 100, d = i)
    lambda.2 <- append(lambda.2, Re(eigen(lambda.i$matrix)$values[1]))}
  
  # Calculate minimum wavespeed
  c.min.2 <- min(c.values.2)
  
  # Append wavespeed to bootstrapped vector of estimated wavespeeds
  boot.cv2 <- append(boot.cv2, c.min.2)
  
  # Clean up
  remove(lambda.i, TM, i)}

# Get procedure time
time.end <- Sys.time()
time.end - time.start
remove(time.start, time.end)





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
               "survival_t1", "volume_t", "weighted.dens", "transplant")) -> boot.CData.AllSurvival





##### Create GLM for survival (transplants included) ------------------------------------------------------

# Remove entries for which there is no recorded volume or survival
# Standardise density; introduce unique transect identifier
boot.CData.AllSurvival.s <- boot.CData.AllSurvival %>% 
  drop_na(volume_t, survival_t1) %>% 
  mutate(d.stand = (weighted.dens - mean(weighted.dens, na.rm = TRUE)) / sd(weighted.dens, na.rm = TRUE),
         unique.transect = interaction(transect, site),
         survival_t1 = as.numeric(survival_t1))

# Create a list of possible survival models
Mod.S <- list()

# "Null" model; only random effects of site and transect within site
Mod.S[[1]] <- glmer(survival_t1 ~ (1 | unique.transect),
                    data = subset(boot.CData.AllSurvival.s, transplant == TRUE), family = "binomial")

# Size-only model
Mod.S[[2]] <- glmer(survival_t1 ~ volume_t + (1|unique.transect),
                    data = subset(boot.CData.AllSurvival.s, transplant == TRUE), family = "binomial")

# Density-only model
Mod.S[[3]] <- glmer(survival_t1 ~ d.stand + (1|unique.transect),
                    data = subset(boot.CData.AllSurvival.s, transplant == TRUE), family = "binomial")

# Size and density (additive)
Mod.S[[4]] <- glmer(survival_t1 ~ volume_t + d.stand + (1|unique.transect),
                    data = subset(boot.CData.AllSurvival.s, transplant == TRUE), family = "binomial")

# Size and density (interactive)
Mod.S[[5]] <- glmer(survival_t1 ~ volume_t * d.stand + (1|unique.transect),
                    data = subset(boot.CData.AllSurvival.s, transplant == TRUE), family = "binomial")

# Density-only model (quadratic)
Mod.S[[6]] <- glmer(survival_t1 ~ d.stand + I(d.stand^2) + (1|unique.transect),
                    data = subset(boot.CData.AllSurvival.s, transplant == TRUE), family = "binomial")

# Size (linear) and density (quadratic)
Mod.S[[7]] <- glmer(survival_t1 ~ volume_t + d.stand + I(d.stand^2) + (1|unique.transect),
                    data = subset(boot.CData.AllSurvival.s, transplant == TRUE), family = "binomial")

# Size and density (interactive, quadratic)
Mod.S[[8]] <- glmer(survival_t1 ~ volume_t * d.stand + volume_t * I(d.stand^2) + (1|unique.transect),
                    data = subset(boot.CData.AllSurvival.s, transplant == TRUE), family = "binomial")

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


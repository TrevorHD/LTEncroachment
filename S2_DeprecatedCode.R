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


##### Initialise Data -------------------------------------------------------------------------------------

# Import dataset
DATA <- "https://github.com/TrevorHD/LTEncroachment/raw/master/LT_Data.xlsx"

# Read transplant data
CData.Transplants <- read.xlsx(DATA, 5)





##### Tidy up transplant data for survival analysis -------------------------------------------------------

# Locations are multiples of 2.5 m, but our density data are in 5-m windows
# Therefore, we will round locations up to nearest 5-m window
for(i in 1:nrow(CData.Transplants)){
  if(CData.Transplants$plot_location[i] %% 2.5 == 0){
    CData.Transplants$plot_location[i] <- CData.Transplants$plot_location[i] + 2.5}}

# Calculate conical volume; use initial volume since most plants die after a year
# Add variable indicating which plants are transplants
# To use best-case scenario recruitment, uncomment the PDC line below
CData.Transplants %>% 
  filter(site == "PDC") %>% 
  mutate("volume_t" = log(vol(h = max.ht_t, w = max.w_t, p = perp.w_t)),
         "transplant" = TRUE) %>% 
  
  # Rename "plot_location" to "actual.window"
  rename("actual.window" = "plot_location") %>% 
  
  # Merge with demography data
  merge(Windows, 
        by.x = c("site", "transect", "actual.window"),
        by.y = c("site", "transect", "window")) -> CData.Transplants

# Combine transplants with large shrubs for later survival analysis
# Keep only location info, survival, volume, and density
select(CData.Transplants, "site", "transect", "actual.window", 
       "spring_survival_t1", "volume_t", "weighted.dens", "transplant") %>% 
  rename("survival_t1" = "spring_survival_t1") %>% 
  rbind(select(CData, "site", "transect", "actual.window", 
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





##### Clean up --------------------------------------------------------------------------------------------

# Clean up variables from global environment
remove(d.and.v, d.and.v.2, d.cut, S.mean.df, S.pred, v.cut, d.range, v.range, i, DATA)


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





##### Survival model averaging ----------------------------------------------------------------------------

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


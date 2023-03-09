## Purpose: trying to model non-normal growth for creosote SIPM
## Author: Tom Miller
## Last modified: 7.28.2021

library(mgcv)
library(sgt)
library(maxLik)
library(zoo)
source("C:/Users/tm9/Dropbox/github/IPM_size_transitions/Diagnostics.R")
# "04_CDataPrep"
# Tidy up demography data before creating demography models
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/04_CDataPrep.R")

  # Create unique transect as interaction of transect and site
  # First-time only; this bit of code in 06_BootRes does it all other times
  # No transplants in this data set
LATR_full <- CData %>% 
  mutate(unique.transect = interaction(transect, site))

## set the gamma argument of gam() -- gamma>1 generates smoother fits, less likely to be overfit
gamma = 1.2





##### Fit growth models -----------------------------------------------------------------------------------

# Prepare a data subset for growth that drops rows missing either t or t1 size data
# Also create log_volume as a new variable because GAM doesn't like functions of variables as variables
LATR_grow <- LATR_full  %>% drop_na(volume_t, volume_t1) %>%
  mutate(log_volume_t = log(volume_t),
         log_volume_t1 = log(volume_t1))

# Create empty list to populate with model results
LATR_gam_models <- list()

# Three candidate models for the mean: size only, size + density, or size, density, and size:density
# Three candidates for variance: size only, size + density, fitted value (all the covariates plus rfx)

# constant sigma
LATR_gam_models[[1]] <- gam(list(log_volume_t1 ~ s(log_volume_t) + s(unique.transect, bs = "re"), ~ 1),
                            data = LATR_grow, gamma = gamma, family = gaulss())
LATR_gam_models[[2]] <- gam(list(log_volume_t1 ~ s(log_volume_t) + s(weighted.dens) + s(unique.transect, bs = "re"), ~ 1), 
                            data = LATR_grow, gamma = gamma, family = gaulss())                
LATR_gam_models[[3]] <- gam(list(log_volume_t1 ~ s(log_volume_t) + s(weighted.dens) + weighted.dens:log_volume_t + s(unique.transect, bs = "re"), ~ 1), 
                            data = LATR_grow, gamma = gamma, family = gaulss()) 
# sigma as function of initial size
LATR_gam_models[[4]] <- gam(list(log_volume_t1 ~ s(log_volume_t) + s(unique.transect,bs = "re"), ~ s(log_volume_t)), 
                            data = LATR_grow, gamma = gamma, family = gaulss())
LATR_gam_models[[5]] <- gam(list(log_volume_t1 ~ s(log_volume_t) + s(weighted.dens) + s(unique.transect, bs = "re"), ~ s(log_volume_t)), 
                            data = LATR_grow, gamma = gamma, family = gaulss())                
LATR_gam_models[[6]] <- gam(list(log_volume_t1 ~ s(log_volume_t) + s(weighted.dens) + weighted.dens:log_volume_t + s(unique.transect, bs = "re"), ~ s(log_volume_t)), 
                            data = LATR_grow, gamma = gamma, family = gaulss()) 
# sigma depends on both initial size and density
LATR_gam_models[[7]] <- gam(list(log_volume_t1 ~ s(log_volume_t) + s(unique.transect, bs = "re"), ~ s(log_volume_t) + s(weighted.dens)), 
                            data = LATR_grow, gamma = gamma, family = gaulss())
LATR_gam_models[[8]] <- gam(list(log_volume_t1 ~ s(log_volume_t) + s(weighted.dens) + s(unique.transect, bs = "re"), ~ s(log_volume_t) + s(weighted.dens)), 
                            data = LATR_grow, gamma = gamma, family = gaulss())                
LATR_gam_models[[9]] <- gam(list(log_volume_t1 ~ s(log_volume_t) + s(weighted.dens) + weighted.dens:log_volume_t + s(unique.transect, bs = "re"), ~ s(log_volume_t) + s(weighted.dens)), 
                            data = LATR_grow, gamma = gamma, family = gaulss()) 

# Collect model AICs into a single table
grow_aic <- AICtab(LATR_gam_models, base = TRUE, sort = FALSE)

# Define best Gaussian model
LATR_gam_model <- LATR_gam_models[[which.min(grow_aic$AIC)]]

# Extract fitted values
LATR_grow$fitted_vals = predict(LATR_gam_model, type = "response", data = LATR_grow)

# Extract the linear predictor for the mean and sd
LATR_Xp <- predict.gam(LATR_gam_model, type = "lpmatrix")

# Fitted coefficients
LATR_beta <- coef(LATR_gam_model)





##### Extract values of the fitted splines to explore their properties ------------------------------------

# Run predict to get fitted terms
fitted_terms <- predict(LATR_gam_model, type = "terms")  

# Effect of initial size on mean of final size 
plot(LATR_grow$log_volume_t, fitted_terms[, "s(log_volume_t)"]) 

# Effect of weighted.dens on mean of final size; complicated
plot(LATR_grow$weighted.dens, fitted_terms[, "s(weighted.dens)"])

# Sigma versus size - presumably log scale
plot(LATR_grow$log_volume_t, fitted_terms[, "s.1(log_volume_t)"]) 

# Sigma versus density 
plot(LATR_grow$weighted.dens, fitted_terms[, "s.1(weighted.dens)"]) 





##### Inspect scaled residuals to evaluate the pilot model (fails!) ---------------------------------------

# Run predict to get fitted terms
fitted_all <- predict(LATR_gam_model, type = "response")  
sigma.hat <- 1/fitted_all[, 2]
scaledResids <- residuals(LATR_gam_model, type = "response")/sigma.hat
par(mfrow = c(1, 2))
plot(fitted_all[, 1], scaledResids, xlab = "Fitted values", ylab = "Scaled residuals") 

# Inspect residuals
qqPlot(scaledResids)          # really bad in both tails
jarque.test(scaledResids)     # normality test: fails, p < 0.001 
anscombe.test(scaledResids)   # kurtosis: fails, p < 0.001 
agostino.test(scaledResids)   # skewness: fails, p < 0.001 

# Rolling NP moments diagnostic; skew is small but variable, tails are fat. 
px <- LATR_grow$fitted_vals[, 1]; py <- scaledResids; 
z <- rollMomentsNP(px, py, windows = 8, smooth = TRUE, scaled = TRUE) 





##### Try to recover GAM parameters with a hand-cranked Gaussian model ------------------------------------

# Set function for Gaussian log likelihood
gausLogLik <- function(pars, response){
  val <- dnorm(x = response, 
               mean = LATR_Xp[, 1:31]%*%pars[1:31],
               sd = exp(LATR_Xp[, 32:50]%*%pars[32:50]),
               log = T) 
  return(val)}

# Set initial parameter values
p0 <- LATR_beta 
paranoid_iter <- 3
coefs <- list(paranoid_iter)
LL <- numeric(paranoid_iter)

# Get likelihood values
for(j in 1:paranoid_iter){
  out <- maxLik(logLik = gausLogLik, start = p0*exp(0.2*rnorm(length(p0))), response = LATR_grow$log_volume_t1,
                method = "BHHH", control = list(iterlim = 5000, printLevel = 2), finalHessian = FALSE)
  out <- maxLik(logLik = gausLogLik, start = out$estimate, response = LATR_grow$log_volume_t1,
                method = "NM", control = list(iterlim = 5000, printLevel = 1), finalHessian = FALSE)
  out <- maxLik(logLik = gausLogLik, start = out$estimate, response = LATR_grow$log_volume_t1,
                method = "BHHH", control = list(iterlim = 5000, printLevel = 2), finalHessian = FALSE)
  coefs[[j]] <- out$estimate
  LL[j] <- out$maximum
  cat(j, "#--------------------------------------#", out$maximum, "\n")}
j <- min(which(LL == max(LL))) 
out <- maxLik(logLik = gausLogLik, start = coefs[[j]], response = LATR_grow$log_volume_t1,
              method = "BHHH", control = list(iterlim = 5000, printLevel = 2), finalHessian = TRUE) 

# Compare to original GAM estimates
plot(LATR_beta[1:31], out$estimate[1:31])
abline(0, 1)
plot(LATR_Xp[, 1:31]%*%LATR_beta[1:31], LATR_Xp[, 1:31]%*%out$estimate[1:31])
abline(0, 1)





##### Fit SGT using design matrices from normal GAM -------------------------------------------------------

# Set function for log likelihood
# sgtLogLik = function(pars, response){
#   val <- dsgt(x = response, 
#               mu = pars[1] + pars[2]*LATR_grow$log_volume_t + pars[3]*LATR_grow$weighted.dens,
#               sigma = exp(pars[4] + pars[5]*LATR_grow$log_volume_t + pars[6]*LATR_grow$weighted.dens),
#               lambda = pars[7],
#               p = exp(pars[8]),
#               q = exp(pars[9]),
#               log = T) 
#   return(val)}

# Set up inverse logit function
invlogit <- function(x){exp(x)/(1 + exp(x))}

# Set function for log likelihood
sgtLogLik <- function(pars, response){
  val <- dsgt(x = response, 
              mu = LATR_Xp[, 1:31]%*%pars[1:31],
              sigma = exp(LATR_Xp[, 32:50]%*%pars[32:50]),
              lambda = -invlogit(pars[51] + pars[52]*LATR_grow$log_volume_t),
              p = exp(pars[53]),
              q = exp(pars[54]),
              mean.cent = T,
              var.adj = T,
              log = T) 
  return(val)}

# Set initial parameter values
p0 <- c(LATR_beta, -10, 0, 2, 2)

# Get likelihood values
out <- maxLik(logLik = sgtLogLik, start = p0*exp(0.2*rnorm(length(p0))), response = LATR_grow$log_volume_t1,
              method = "BHHH", control = list(iterlim = 5000, printLevel = 2), finalHessian = FALSE)
out <- maxLik(logLik = sgtLogLik, start = out$estimate, response = LATR_grow$log_volume_t1,
              method = "NM", control = list(iterlim = 5000, printLevel = 1), finalHessian = FALSE)
out <- maxLik(logLik = sgtLogLik, start = out$estimate, response = LATR_grow$log_volume_t1,
              method = "BHHH", control = list(iterlim = 5000, printLevel = 2), finalHessian = FALSE)
out <- maxLik(logLik = sgtLogLik, start = out$estimate, response = LATR_grow$log_volume_t1,
              method = "BHHH", control = list(iterlim = 5000, printLevel = 2), finalHessian = TRUE) 

# Compare to original GAM estimates
plot(LATR_beta, out$estimate[1:50])
abline(0, 1)
plot(LATR_Xp[, 1:31]%*%LATR_beta[1:31], LATR_Xp[, 1:31]%*%out$estimate[1:31])
abline(0, 1)

## I'm satisfied that the SGT can recover the same expected value as the gaussian gam()





##### Compare simulated and real data ---------------------------------------------------------------------

# Simulate data from normal and SGT models
n_sim <- 500
LATR_sim_NO <- LATR_sim_SGT <- matrix(NA, nrow = nrow(LATR_grow), ncol = n_sim)
for(i in 1:n_sim){
  print(i)
  LATR_sim_SGT[, i] <- rsgt(n = nrow(LATR_grow), 
                            mu = LATR_Xp[, 1:31]%*%out$estimate[1:31], 
                            sigma = exp(LATR_Xp[, 32:50]%*%out$estimate[32:50]),
                            lambda = -invlogit(out$estimate[51] + out$estimate[52]*LATR_grow$log_volume_t),
                            p = exp(out$estimate[53]),
                            q = exp(out$estimate[54]),
                            mean.cent = T,
                            var.adj = T)
  LATR_sim_NO[, i] <- rnorm(n = nrow(LATR_grow),
                            mean = fitted_all[, 1],
                            sd = 1/fitted_all[, 2])}

# Bin data for plotting
n_bins <- 12
alpha_scale <- 0.7
LATR_moments <- LATR_grow %>% 
  arrange(log_volume_t) %>% 
  mutate(size_bin = cut_interval(log_volume_t, n = n_bins)) %>% 
  group_by(size_bin) %>% 
  summarise(mean_t1 = mean(log_volume_t1),
            sd_t1 = sd(log_volume_t1),
            skew_t1 = NPskewness(log_volume_t1),
            kurt_t1 = NPkurtosis(log_volume_t1),
            bin_mean = mean(log_volume_t),
            bin_n = n()) 

par(mfrow = c(2, 2), mar = c(4, 4, 2, 1), cex.axis = 1.3, cex.lab = 1.3, mgp = c(2, 1, 0), bty = "l"); 
sim_bin_means <- sim_moment_means <- sim_moment_means_norm <- matrix(NA, n_bins, n_sim); 
for(i in 1:n_sim){
  sim_moments <- bind_cols(LATR_grow, data.frame(sim = LATR_sim_SGT[, i],
                                                 sim_norm = LATR_sim_NO[, i])) %>% 
    arrange(log_volume_t) %>% 
    mutate(size_bin = cut_interval(log_volume_t, n = n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_t1 = mean(sim),
              mean_t1_norm = mean(sim_norm),
              bin_mean = mean(log_volume_t))
  sim_bin_means[, i] <- sim_moments$bin_mean; 
  sim_moment_means[, i] <- sim_moments$mean_t1
  sim_moment_means_norm[, i] <- sim_moments$mean_t1_norm}
matplot(LATR_moments$bin_mean, sim_moment_means, col = alpha("gray", 0.5), pch = 16,
        xlab = "Mean size t0", ylab = "Mean(Size t1)", cex = 1,
        xlim = c(min(LATR_moments$bin_mean), max(LATR_moments$bin_mean) + 0.4))
matplot(LATR_moments$bin_mean + 0.4, sim_moment_means_norm,
        col = alpha("cornflowerblue", 0.5), pch = 16, add = T)
points(LATR_moments$bin_mean + 0.2, LATR_moments$mean_t1, pch = 16, lwd = 2,
       col = alpha("red", alpha_scale), cex = 1.6)
points(LATR_moments$bin_mean, apply(sim_moment_means, 1, median), pch = 1, lwd = 2,
       col = alpha("black", alpha_scale), cex = 1.6)
points(LATR_moments$bin_mean + 0.4, apply(sim_moment_means_norm, 1, median), pch = 1, lwd = 2,
       col = alpha("black", alpha_scale), cex = 1.6)
legend("topleft", legend = c("SGT", "Gaussian", "Data"),
       col = c("gray", "cornflowerblue", "red"), pch = 16, bty = "n", cex = 1.4, pt.lwd = 2, pt.cex = 1.6) 
add_panel_label("a")

for(i in 1:n_sim){
  sim_moments <- bind_cols(LATR_grow, data.frame(sim = LATR_sim_SGT[, i],
                                                 sim_norm = LATR_sim_NO[, i])) %>% 
    arrange(log_volume_t) %>% 
    mutate(size_bin = cut_interval(log_volume_t, n = n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_t1 = sd(sim),
              mean_t1_norm = sd(sim_norm),
              bin_mean = mean(log_volume_t))
  sim_bin_means[, i] <- sim_moments$bin_mean; 
  sim_moment_means[, i] <- sim_moments$mean_t1
  sim_moment_means_norm[,i] <- sim_moments$mean_t1_norm}
matplot(LATR_moments$bin_mean, sim_moment_means,col = alpha("gray", 0.5), pch = 16,
        xlab = "Mean size t0", ylab = "SD(Size t1)", cex = 1,
        xlim = c(min(LATR_moments$bin_mean), max(LATR_moments$bin_mean) + 0.4)) 
matplot(LATR_moments$bin_mean + 0.4, sim_moment_means_norm,
        col = alpha("cornflowerblue", 0.5), pch = 16, add = T)
points(LATR_moments$bin_mean + 0.2, LATR_moments$sd_t1, pch = 16, lwd = 2,
       col = alpha("red", alpha_scale), cex = 1.6)
points(LATR_moments$bin_mean, apply(sim_moment_means, 1, median), pch = 1, lwd = 2,
       col = alpha("black", alpha_scale), cex = 1.6)
points(LATR_moments$bin_mean + 0.4, apply(sim_moment_means_norm, 1, median), pch = 1, lwd = 2,
       col = alpha("black", alpha_scale), cex = 1.6)
add_panel_label("b")

for(i in 1:n_sim){
  sim_moments <- bind_cols(LATR_grow, data.frame(sim = LATR_sim_SGT[, i],
                                                 sim_norm = LATR_sim_NO[, i])) %>% 
    arrange(log_volume_t) %>% 
    mutate(size_bin = cut_interval(log_volume_t, n = n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_t1 = NPskewness(sim),
              mean_t1_norm = NPskewness(sim_norm),
              bin_mean = mean(log_volume_t))
  sim_bin_means[, i] <- sim_moments$bin_mean; 
  sim_moment_means[, i] <- sim_moments$mean_t1
  sim_moment_means_norm[, i] <- sim_moments$mean_t1_norm}
matplot(LATR_moments$bin_mean, sim_moment_means, col = alpha("gray", 0.5), pch = 16,
        xlab = "Mean size t0", ylab = "Skew(Size t1)", cex = 1,
        xlim = c(min(LATR_moments$bin_mean), max(LATR_moments$bin_mean) + 0.4))
matplot(LATR_moments$bin_mean + 0.4, sim_moment_means_norm,
        col = alpha("cornflowerblue", 0.5),pch = 16, add = T)
points(LATR_moments$bin_mean + 0.2, LATR_moments$skew_t1, pch = 16, lwd = 2,
       col = alpha("red", alpha_scale), cex = 1.6)
points(LATR_moments$bin_mean, apply(sim_moment_means, 1, median), pch = 1, lwd = 2,
       col = alpha("black", alpha_scale), cex = 1.6)
points(LATR_moments$bin_mean + 0.4, apply(sim_moment_means_norm, 1, median), pch = 1, lwd = 2,
       col = alpha("black", alpha_scale), cex = 1.6)
add_panel_label("c")

for(i in 1:n_sim){
  sim_moments <- bind_cols(LATR_grow, data.frame(sim = LATR_sim_SGT[, i],
                                                 sim_norm = LATR_sim_NO[, i])) %>% 
    arrange(log_volume_t) %>% 
    mutate(size_bin = cut_interval(log_volume_t, n = n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_t1 = NPkurtosis(sim),
              mean_t1_norm = NPkurtosis(sim_norm),
              bin_mean = mean(log_volume_t))
  sim_bin_means[, i] <- sim_moments$bin_mean; 
  sim_moment_means[, i] <- sim_moments$mean_t1
  sim_moment_means_norm[, i] <- sim_moments$mean_t1_norm}
matplot(LATR_moments$bin_mean, sim_moment_means, col = alpha("gray", 0.5), pch = 16,
        xlab = "Mean size t0", ylab = "Kurtosis(Size t1)", cex = 1,
        xlim = c(min(LATR_moments$bin_mean), max(LATR_moments$bin_mean) + 0.4))
matplot(LATR_moments$bin_mean + 0.4, sim_moment_means_norm,
        col = alpha("cornflowerblue", 0.5), pch = 16, add = T)
points(LATR_moments$bin_mean + 0.2, LATR_moments$kurt_t1, pch = 16, lwd = 2,
       col = alpha("red", alpha_scale), cex = 1.6)
points(LATR_moments$bin_mean, apply(sim_moment_means, 1, median), pch = 1, lwd = 2,
       col = alpha("black", alpha_scale), cex = 1.6)
points(LATR_moments$bin_mean + 0.4, apply(sim_moment_means_norm, 1, median), pch = 1, lwd = 2,
       col = alpha("black", alpha_scale), cex = 1.6)
add_panel_label("d")


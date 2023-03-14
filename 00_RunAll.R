##### Load required packages for all scripts --------------------------------------------------------------

# Load packages in this order
library(minpack.lm)
library(MASS)
library(tidyverse)
library(SuppDists)
library(sqldf)
library(lme4)
library(bbmle)
library(truncnorm)
library(grid)
library(gridBase)
library(corrplot)
library(mgcv)
library(popbio)
library(bbmle)
library(scales)
library(wesanderson)
library(oce)
library(sgt)
library(maxLik)
library(zoo)
library(binr)
library(plotrix)
library(ggbreak)
library(patchwork)





##### Run scripts that do not change between the two scenarios --------------------------------------------

# This will take several minutes; be patient

# "01_SeedVelocities"
# Calculate terminal velocities of seeds and fit to distributions
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/01_SeedVelocities.R")

# "02_WindSpeeds"
# Load in wind speeds and fit to distributions
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/02_WindSpeeds.R")

# "03_Dispersal"
# Construct dispersal kernel functions for seeds
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/03_Dispersal.R")

# "04_CDataPrep"
# Tidy up demography data before creating demography models
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/04_CDataPrep.R")

# "07_SIPM"
# Spatial integral projection setting up functions to calculate wavespeeds
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/07_SIPM.R")

# Pre-bootstrap switch that will flip once bootstrapping begins
# Done so that some parts of code are not re-run unnecessarily
# Should not be modified by the user
boot.switch <- FALSE

# "05_CDataAnalysis"
# Create demography models for use in SIPM
# Run once before bootstrapping to get recruit sizes and boundaries
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/05_CDataAnalysis.R")





##### Set up bootstrapping for wavespeeds -----------------------------------------------------------------

# Save bootstrapped wavespeed and lambda outputs to CSV?
boot.saveOutputs <- TRUE

# Should bootstrapping occur?
# If not, the model will be run once using full data sets
boot.on <- FALSE

# Evaluate local IPM instead of spatial IPM?
# Local IPM will not include dispersal
boot.noDisp <- FALSE

# What proportion of individuals should be resampled each bootstrap iteration?
# Should be in the interval (0, 1), exclusive
# Warning: setting this number very low may adversely affect model behaviour
# Ignore this if boot.on = FALSE
boot.prop <- 0.75

# Set number of bootstrap iterations
# Please note: one iteration may take a few minutes depending on computer and settings
# Ignore this if boot.on = FALSE
boot.num <- 500

# Create empty vectors to populate with wavespeeds
# c1 are analytic wavespeeds using mean wind speed and terminal velocity, assuming H as seed point source
# c2 are wavespeeds from simulating dispersal events over variation in heights, wind speed, and terminal velocity
boot.c1 <- c()
boot.c2 <- c()

# Should perturbation analysis be run?
pert.on <- FALSE

# Set magnitude of perturbation; 0.01 increase in vital rate function
pert <- 0.01

# Set vector of vital rates to be perturbed
elas <- c("growth.mean", "growth.sd", "survival", "flower", "fertility", "recruitment", "recruitsize.mean",
          "recruitsize.sd", "dispersal.location", "dispersal.scale")

# Store elasticity values across bootstrap replicates
boot.elas <- vector("numeric", length = length(elas))
boot.sens <- vector("numeric", length = length(elas))

# Store seed data
# seeds <- sample.int(100000, size = boot.num); write.csv(seeds, "Data/Derived/Seed_1000.csv")

# Read seed data
seeds <- read.csv("Data/Derived/Seed_1000.csv")
seed.store <- c()





##### Wavespeeds and population growth for normal survival scenario ---------------------------------------

# This takes several minutes per bootstrap replicate; be patient
# A stable internet connection is required
# Note: if boot.on = FALSE, then bootstrapping will not occur and full data will be used

# Override bootstrap replicate number if bootstrapping is turned off
if(boot.on == FALSE){
  boot.num <- 1}

# Begin bootstrapping
time.start <- Sys.time()
for(i in 1:boot.num){
  
  # Flip switch at beginning of bootstrapping
  # Failure to do so should not affect results, but WILL increase computation time
  if(i == 1){
    boot.switch <- TRUE}
  
  # "06_BootRes"
  # Run resampling subroutine for wind speeds, terminal velocities, and demography
  source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/06_BootRes.R")
  
  # "05_CDataAnalysis_NS.R"
  # Create demography models for use in SIPM, using best models from initial run of this file
  # This avoids model uncertainty; only model coefficients change, not the overall structure
  source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/05_CDataAnalysis.R")
  
  # Construct transition matrix for minimum weighted density (zero)
  TM <- TransMatrix(dens = 0, mat.size = 100)

  # Evaluate SIPM to find wavespeeds
  if(boot.noDisp == FALSE){
    
    # Get WALD parameters from this data bootstrap
    # Generates a list of length mat.size containing plant heights and corresponding WALD parameters
    params <- WALD_par()
    
    # Sample dispersal events for empirical MGF; generates a N*mat.size matrix
    D.samples <- WALD_samples(N = 10000, seed = seeds[i, 2]) 
    
    # Find the asymptotic wave speed c*(s) 
    c1star <- optimize(cs, lower = 0.05, upper = 4, emp = F)$objective
    
    # Use empirical MGF to get wavespeed from sampled dispersal events
    c2star <- optimize(cs, lower = 0.05, upper = 4, emp = T)$objective
    
    # Calculate minimum wavespeed, then append to bootstrapped vector of estimated wavespeeds
    boot.c1 <- append(boot.c1, c1star)
    boot.c2 <- append(boot.c2, c2star)
    
    # Run perturbation analysis on this data bootstrap
    if(pert.on == TRUE){
      
      # Set up elasticity vector
      # c2star.elas <- c()
      
      # Set up sensitivity vector and populate
      c2star.sens <- c()
      for(e in 1:length(elas)){
        
        # Estimate elasticities
        # TM <- TransMatrix(dens = 0, mat.size = 100, elas = elas[e])
        # D.samples <- WALD_samples(N = 10000, seed = seeds[i, 2], elas = elas[e]) 
        # c2star.elas[e] <- optimize(cs, lower = 0.05, upper = 4, emp = T)$objective
        
        # Estimate senstitivies
        TM <- TransMatrix(dens = 0, mat.size = 100,sens = elas[e])
        D.samples <- WALD_samples(N = 10000, seed = seeds[i, 2], sens = elas[e]) 
        c2star.sens[e] <- optimize(cs, lower = 0.05, upper = 4, emp = T)$objective}

      # Calculate proportional change in cstar and append to output for elasticity
      # boot.elas <- rbind(boot.elas, ((c2star.elas/c2star) - 1))
      
      # Calculate proportional change in cstar and append to output for sensitivity
      boot.sens <- rbind(boot.sens, (c2star.sens - c2star)/pert)}}

  # Uncomment code below to run analysis for lambda versus density
  # Create empty list to store lambda as a function of density; assign density values to top
  # if(i == 1){
  #  boot.lambda <- list(density = LambdaD(d.only = TRUE))}
  #
  # Calculate lambda as a function of density, then append to list
  # boot.lambda[[i + 1]] <- LambdaD()
  
  # Create empty list to store transition matrices
  if(i == 1){
    boot.TM <- list()}
  
  # Append transition matrix to list
  boot.TM[[i]] <- TM
  
  # Store seed data
  seed.store[i]<-seeds[i, 2]
  
  # Calculate elapsed time
  time.elapsed <- as.numeric(difftime(Sys.time(), time.start, units = "hours"))
  
  # Clear console (on Windows) and print bootstrapping progress to console
  shell("cls")
  cat(paste0(i, "/", boot.num, " (", round(i/boot.num, 3)*100, "%) complete; ",
               round(time.elapsed, 2), " hours elapsed."))
  
  # Flip switch back to original setting
  if(i == boot.num){
    boot.switch <- FALSE}}

# Clear console (on Windows) and print final procedure time
# Suppress any warnings or errors
suppressWarnings(if(1 == 1){
  shell("cls")
  try(ifelse(length(boot.TM) == boot.num,
             cat(paste0("Procedure 100% complete with total time of ", round(time.elapsed, 2), " hours.")),
             message(paste0("Warning: Network connection interrupted; bootstrapping could not be fully completed. \n",
                            "Procedure halted at ", round(length(boot.cv1)/boot.num, 3)*100, "% with total time of ",
                            round(time.elapsed, 2), " hours."))), silent = TRUE)})

# Save outputs to CSV if enabled
if(boot.saveOutputs == TRUE){
  
  # Write bootstrapped lambda values to csv
  write.csv(boot.lambda, "Data/Derived/Boot_Lambda.csv")
  write.csv(seed.store, "Data/Derived/Seed_Store.csv")
  
  # Write bootstrapped wavespeed values to csv
  write.csv(boot.c1, "Data/Derived/Boot_C1.csv")
  write.csv(boot.c2, "Data/Derived/Boot_C2.csv")
  # write.csv(boot.elas, "Boot_Elas.csv")
  write.csv(boot.sens, "Data/Derived/Boot_Sens.csv")}

# Remove unneeded bootstrap items from the global environment if using spatial IPM
# If running single replicate, just leave most items in global environment
# Suppress errors since some objects may not exist depending on which parts of code are re-run
if(boot.noDisp == FALSE){
  try(remove(boot.LATR_recruit_size, boot.tv.PDF, boot.ws.PDF, flower_aic, fruits_aic, grow_aic, LATR_dat_201718,
             LATR_flow_dat, LATR_flower, LATR_flower_best, LATR_flower_fitted_terms, LATR_fruits, LATR_fruits_best,
             LATR_fruits_dat, LATR_fruits_fitted_terms, LATR_gam_models, LATR_grow, LATR_grow_best, LATR_grow_fitted_terms,
             LATR_recruit, LATR_recruit_best, LATR_recruitment, LATR_recruits, LATR_surv, LATR_surv_best, LATR_surv_dat,
             LATR_surv_fitted_terms, recruit_aic, surv_aic, TM, boot.tv.raw, boot.ws.raw, c.values, grow_sd_index, i,
             time.elapsed, time.start, TM.matdim, TM.lower.extension, TM.upper.extension, gamma, seeds_per_fruit),
      silent = TRUE)}
if(boot.noDisp == TRUE){
  try(remove(boot.tv.PDF, boot.ws.PDF, boot.tv.raw, boot.ws.raw, i, time.elapsed, time.start, TM),
      silent = TRUE)}





##### Run perturbation analysis with full data ------------------------------------------------------------

# Run select data prep scripts again
boot.on <- FALSE
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/04_CDataPrep.R")
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/05_CDataAnalysis.R")
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/06_BootRes.R")
 
# Piece together the H matrix evaluated at s*
TM <- TransMatrix(dens = 0, mat.size = 100)

# Survival-growth matrix
P <- TM$Pmat 

# Fertility matrix
Fs <- TM$Fmat 

# Get dispersal pieces
params <- WALD_par()

# Sample dispersal events for empirical MGF; generates a N*mat.size matrix
D.samples <- WALD_samples(N = 10000, seed = seeds[1, 2]) 

# Get s* and c(s*)
c2star <- optimize(cs, lower = 0.05, upper = 4, emp = T)
Ms <- matrix(NA, nrow(Fs), ncol(Fs))

# Apply MGF over columns of Fs
for(j in 1:length(params$heights)){
  if(params$heights[j] < 0.15){next}
  Ms[, j] <- empiricalWALDmgf(s = c2star$minimum, D.samples[, j])
  Fs[, j] <- Fs[, j]*Ms[, j]}
Hs <- P + Fs 

# Run eigenanalysis (see Ellner et al. Ch4)
Hs.eig <- eigen(Hs)
rho1 <- Re(Hs.eig$values[1])
w.z <- Re(Hs.eig$vectors[, 1])
v.z1 <- Re(eigen(t(Hs))$vectors[, 1])
h <- diff(TM$meshpts[1:2])

# Calculate kernel sensitivity
dc.dK <- (outer(v.z1, w.z, "*")/sum(v.z1*w.z*h))*(1/(c2star$minimum*rho1))
image(TM$meshpts, TM$meshpts, t(dc.dK),
      col = grey(seq(0.6, 1, length = 100)),
      xlab = "Volume (t), z", ylab = "Volume (t+1), z\'")
contour(TM$meshpts, TM$meshpts, t(dc.dK), add = TRUE)





##### Generate main figures -------------------------------------------------------------------------------

# "08_MainFigures"
# Generate figures for wavespeeds and population growth, dispersal, and demographic data
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/08_MainFigures.R")





##### Write files for EDI data package --------------------------------------------------------------------

# 1. Observational demography census
# write.csv(CData %>% 
#   dplyr::select(-c(volume_t, volume_t1, logGR, did.flower, transplant,
#                    total.fruits_t, total.flowers_t,
#                    total.fruits_t1, total.flowers_t1,
#                    total.reproduction_t, total.reproduction_t1)), "EDI/LATR_demography_dat.csv", row.names = F)

# 2. Transplant experiment
# write.csv(CData.Transplants %>% 
#   dplyr::select(-c(name, volume_t, transplant, notes)) %>% 
#     rename(num_black_grama_t = num_black_gramma_t,
#            num_blue_grama_t = num_blue_gramma_t), "EDI/LATR_transplants_dat.csv", row.names = F)

# 3. Seed drop summary data
# write.csv(read.csv("https://github.com/TrevorHD/LTEncroachment/raw/master/Data/SD_Summary.csv"), "EDI/SD_Summary.csv")

# 4. Seed drop position data
# write.csv(read.csv("https://github.com/TrevorHD/LTEncroachment/raw/master/Data/SD_Trials.csv"), "EDI/SD_Trials.csv")


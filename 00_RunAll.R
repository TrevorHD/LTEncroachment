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





##### Set up bootstrapping for wavespeeds -----------------------------------------------------------------

# Pre-bootstrap switch that will flip once bootstrapping begins
# Done so that some parts of code are not re-run unnecessarily
# Should not be modified by the user
boot.switch <- FALSE

# Save bootstrapped wavespeed and lambda outputs to CSV?
boot.saveOutputs <- FALSE

# "05_CDataAnalysis_NS.R"
# Create demography models for use in SIPM
# Run once before bootstrapping to get recruit sizes and boundaries
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/05_CDataAnalysis.R")

# Should bootstrapping occur?
# If not, the model will be run once using full data sets
boot.on <- TRUE

# Evaluate local IPM instead of spatial IPM?
# Local IPM will not include dispersal
boot.noDisp <- FALSE

# What proportion of individuals should be resampled each bootstrap interation?
# Should be in the interval (0, 1), exclusive
# Warning: setting this number very low may adversely affect model behaviour
# Ignore this if boot.on = FALSE
boot.prop <- 0.75

# Set number of bootstrap iterations
# Please note: one iteration takes some time (5-15 minutes) depending on computer and settings
# Ignore this if boot.on = FALSE
boot.num <- 2

# Create empty vectors to populate with wavespeeds
boot.cv1 <- c()





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
  
  # "07_SIPM"
  # Spatial integral projection setting up functions to calculate wavespeeds
  source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/07_SIPM.R")
  
  # Evaluate SIPM to find wavespeeds
  if(boot.noDisp == FALSE){
    
    # Create empty vector to store wavespeeds
    c.values <- Wavespeed()
    
    # Calculate minimum wavespeed, then append to bootstrapped vector of estimated wavespeeds
    boot.cv1 <- append(boot.cv1, min(c.values))}
  
  # Create empty list to store lambda as a function of density; assign density values to top
  if(i == 1){
    boot.lambda <- list(density = LambdaD(d.only = TRUE))}
  
  # Calculate lambda as a function of density, then append to list
  boot.lambda[[i + 1]] <- LambdaD()
  
  # Create empty list to store transition matrices
  if(i == 1){
    boot.TM <- list()}
  
  # Append transition matrix to list
  boot.TM[[i]] <- TM
  
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
  write.csv(boot.lambda, "BootLambda.csv")
  
  # Write bootstrapped wavespeed values to csv
  write.csv(boot.cv1, "BootCV.csv")}

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





##### Generate main figures -------------------------------------------------------------------------------

# "08_MainFigures"
# Generate figures for wavespeeds and population growth, dispersal, and demographic data
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/08_MainFigures.R")


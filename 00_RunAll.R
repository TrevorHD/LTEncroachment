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
# Done so that computationally-expensive parts of 05_CDataAnalysis_NS are not re-run unnecessarily
boot.switch <- FALSE

# "05_CDataAnalysis_NS.R"
# Create demography models for use in SIPM
# Run once before bootstrapping, using full data set to find the best models
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/05_CDataAnalysis.R")

# Should bootstrapping occur?
# If not, the model will be run once using full data sets
boot.on <- TRUE

# Determine resampling proportion; that is, the proportion of individuals selected each interation
# Should be in the interval (0, 1), exclusive
# Warning: setting this number very low may adversely affect model behaviour
boot.prop <- 0.80

# Set number of bootstrap iterations
# Please note: one iteration takes some time (5-10 minutes), so choose this number wisely
# Ignore this if boot.on = FALSE
boot.num <- 2

# Create empty vectors to populate with wavespeeds for normal and higher survival scenarios
boot.cv1 <- c()





##### Wavespeeds and population growth for normal survival scenario ---------------------------------------

# This takes 5-10 minutes per bootstrap replicate; be patient
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
  
  # Wavespeeds as function of s; growth as function of density
  c.values <- Wavespeed()
  
  # Calculate minimum wavespeed
  c.min <- min(c.values)
  
  # Append wavespeed to bootstrapped vector of estimated wavespeeds
  boot.cv1 <- append(boot.cv1, c.min)
  
  # Calculate elapsed time
  time.elapsed <- as.numeric(difftime(Sys.time(), time.start, units = "hours"))
  
  # Clear console (on Windows) and print bootstrapping progress to console
  shell("cls")
  print(paste0(i, "/", boot.num, " (", round(i/boot.num, 3)*100, "%) complete; ",
               round(time.elapsed, 2), " hours elapsed."))
  
  # Flip switch back to original setting
  if(i == boot.num){
    boot.switch <- FALSE}}

# Clear console (on Windows) and print final procedure time
shell("cls")
print(paste0("Procedure complete with total time of ", time.elapsed, " hours." ))
remove(time.start, time.elapsed, time.end)

# Remove other unneeded items from the global environment
# Don't worry if this returns errors; it just means the items were already cleared
remove(fitGAU, fitted_all, err, fitted_vals, i, k, n_cuts_dens, new_fitted_vals, new_weights,
       weights, boot.tv.raw, boot.tv.PDF, boot.ws.raw, boot.ws.PDF, c.values, c.min, mod)





##### Generate main figures -------------------------------------------------------------------------------

# "07_MainFigures"
# Generate figures for wavespeeds and population growth, dispersal, and demographic data
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/08_MainFigures.R")


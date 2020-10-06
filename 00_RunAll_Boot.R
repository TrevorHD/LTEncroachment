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





# ----- Set up bootstrapping for wavespeeds ---------------------------------------------------------------

# Determine resampling proportion; that is, the proportion of individuals selected each interation
# Should be in the interval (0, 1), exclusive
# Warning: setting this number very low may adversely affect model behaviour
boot.prop <- 0.75

# Set number of bootstrap iterations
# Please note: one iteration takes a long time (~30 minutes), so choose this number wisely
boot.num <- 5

# Create empty vectors to populate with wavespeeds for normal and higher survival scenarios
boot.cv1 <- c()
boot.cv2 <- c()





# ----- Wavespeeds and population growth for normal survival scenario -------------------------------------

# This will take several minutes; be patient

# Begin bootstrapping
time.start <- Sys.time()
for(i in 1:boot.num){
  
  # "00_BootRes"
  # Run resampling subroutine for wind speeds, terminal velocities, and demography
  source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/00_BootRes.R")
  
  # "05_CDataAnalysis_NS"
  # Create demography models for growth, reproduction, survival, etc. under normal circumstances
  source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/05_CDataAnalysis_NS.R")
  
  # "06_SIPM"
  # Spatial integral projection setting up functions to calculate wavespeeds
  source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/06_SIPM.R")
  
  # Wavespeeds as function of s; growth as function of density
  c.values <- Wavespeed(200)
  lambda <- c()
  for(i in seq(-1.3, max(boot.CData.s$d.stand), length.out = 100)){
    lambda.i <- TransMatrix(n = 200, d = i)
    lambda <- append(lambda, Re(eigen(lambda.i$matrix)$values[1]))}
  
  # Calculate minimum wavespeed
  c.min <- min(c.values)
  
  # Append wavespeed to bootstrapped vector of estimated wavespeeds
  boot.cv1 <- append(boot.cv1, c.min)
  
  # Clean up
  remove(lambda.i, TM, i)}

# Get procedure time
time.end <- Sys.time()
time.start - time.end
remove(time.start, time.end)





# ----- Wavespeeds and population growth for higher survival scenario -------------------------------------

# This will take several minutes; be patient

# Begin bootstrapping
time.start <- Sys.time()
for(i in 1:boot.num){
  
  # "00_BootRes"
  # Run resampling subroutine for wind speeds, terminal velocities, and demography
  source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/00_BootRes.R")
  
  # "05_CDataAnalysis_NS"
  # Create demography models for growth, reproduction, survival, etc. under normal circumstances
  # Must run this before 05_CDataAnalysis_BS since it contains all of the demography models
  source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/05_CDataAnalysis_NS.R")
  
  # "05_CDataAnalysis_BS"
  # Replace survival model in 05_CDataAnalysis_NS with higher survival from above-average rainfall
  source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/05_CDataAnalysis_BS.R")
  
  # "06_SIPM"
  # Spatial integral projection model that calculates wavespeeds
  source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/06_SIPM.R")
  
  # Wavespeeds as function of s; growth as function of density
  c.values.2 <- Wavespeed(200)
  lambda.2 <- c()
  for(i in seq(-1.3, max(boot.CData.s$d.stand), length.out = 100)){
    lambda.i <- TransMatrix(n = 200, d = i)
    lambda.2 <- append(lambda.2, Re(eigen(lambda.i$matrix)$values[1]))}
  
  # Calculate minimum wavespeed
  c.min.2 <- min(c.values.2)
  
  # Append wavespeed to bootstrapped vector of estimated wavespeeds
  boot.cv2 <- append(boot.cv2, c.min.2)
  
  # Clean up
  remove(lambda.i, TM, i)}

# Get procedure time
time.end <- Sys.time()
time.start - time.end
remove(time.start, time.end)


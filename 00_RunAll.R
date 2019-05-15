##### Load required packages for all scripts --------------------------------------------------------------

# Load packages in this order
library(xlsx)
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





##### Run scripts that do not change between the two scenarios --------------------------------------------

# This will take several minutes; be patient

# "01_SeedVelocities"
# Calculate terminal velocities of seeds and fit to distributions
source("C:\\Users/Trevor Drees/Dropbox/Trevor Creosote analysis/Code/01_SeedVelocities.R")

# "02_WindSpeeds"
# Load in wind speeds and fit to distributions
source("C:\\Users/Trevor Drees/Dropbox/Trevor Creosote analysis/Code/02_WindSpeeds.R")

# "03_Dispersal"
# Construct dispersal kernels for seeds
source("C:\\Users/Trevor Drees/Dropbox/Trevor Creosote analysis/Code/03_Dispersal.R")

# "04_CDataPrep"
# Tidy up demography data before creating demography models
source("C:\\Users/Trevor Drees/Dropbox/Trevor Creosote analysis/Code/04_CDataPrep.R")





# ----- Wavespeeds and population growth for normal survival scenario -------------------------------------

# This will take several minutes; be patient

# "05_CDataAnalysis_NS"
# Create demography models for growth, reproduction, survival, etc. under normal circumstances
source("C:\\Users/Trevor Drees/Dropbox/Trevor Creosote analysis/Code/05_CDataAnalysis_NS.R")

# "06_SIPM"
# Spatial integral projection model that calculates wavespeeds
source("C:\\Users/Trevor Drees/Dropbox/Trevor Creosote analysis/Code/06_SIPM.R")

# Wavespeeds as function of s; growth as function of density
c.values <- Wavespeed(200)
lambda <- c()
for(i in seq(-1.3, max(CData.s$d.stand), length.out = 100)){
  lambda.i <- TransMatrix(n = 200, d = i)
  lambda <- append(lambda, Re(eigen(lambda.i$matrix)$values[1]))}

# Clean up
remove(lambda.i, TM, i)





# ----- Wavespeeds and population growth for higher survival scenario -------------------------------------

# This will take several minutes; be patient

# "05_CDataAnalysis_BS"
# Replace survival model with one for year with higher survival from above-average rainfall
source("C:\\Users/Trevor Drees/Dropbox/Trevor Creosote analysis/Code/05_CDataAnalysis_BS.R")

# "06_SIPM"
# Spatial integral projection model that calculates wavespeeds
source("C:\\Users/Trevor Drees/Dropbox/Trevor Creosote analysis/Code/06_SIPM.R")

# Wavespeeds as function of s; growth as function of density
c.values.2 <- Wavespeed(200)
lambda.2 <- c()
for(i in seq(-1.3, max(CData.s$d.stand), length.out = 100)){
  lambda.i <- TransMatrix(n = 200, d = i)
  lambda.2 <- append(lambda.2, Re(eigen(lambda.i$matrix)$values[1]))}

# Clean up
remove(lambda.i, TM, i)





# ----- Generate main figures -----------------------------------------------------------------------------

# "07_MainFigures"
# Generate figures for wavespeeds and population growth, dispersal, and demographic data
source("C:\\Users/Trevor Drees/Dropbox/Trevor Creosote analysis/Code/07_MainFigures.R")


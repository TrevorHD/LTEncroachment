##### Initialise data -------------------------------------------------------------------------------------------------------------------------------

# Code author(s): Trevor, Tom

# Goal of this script is to quantify how wavespeed responds to parameter perturbations

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

# Run necessary scripts to prep for wavespeed calculations
# See "RunALL" for descriptions
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/01_SeedVelocities.R")
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/02_WindSpeeds.R")
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/03_Dispersal.R")
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/04_CDataPrep.R")
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/05_CDataAnalysis.R")

# Turn off bootstrapping
boot.on <- FALSE
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/06_BootRes.R")

# Spatial integral projection setting up functions to calculate wavespeeds
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/07_SIPM.R")

# How much variability is there across calculations of wave speed?
# for(i in 1:10){
#   c.values <- Wavespeed()
#   Calculate minimum wavespeed
#   print(min(c.values))}





##### Calculate recruitment sensitivity -------------------------------------------------------------------------------------------------------------

# Loop over variation in recruitment rate (a constant in the SIPM) and recalculate wavespeed

# Initialise vectors
recruitment <- 10:0
lambda.recruitment <- numeric(length = length(recruitment))
wavespeed.recruitment <- numeric(length = length(recruitment))

# Start loop
for(i in 1:length(recruitment)){
  
  # Redefine recruitment function - should overwrite the function defined in 07_SIPM
  TM.recruitment <- function(d = NULL){
    return(10^(-recruitment[i]))}
  
  # Combined flowering, fertility, and recruitment
  TM.fertrecruit <- function(x, y, d){
    TM.flower(x, d) * TM.seeds(x, d) * TM.recruitment(d) * TM.recruitsize(y)}
  
  # Put it all together; projection matrix is a function of weighted density (dens)
  # We need a large lower extension because growth variance (Gaussian) is greater for smaller plants
  TransMatrix <- function(dens, ext.lower = TM.lower.extension, ext.upper = TM.upper.extension,
                          min.size = LATR_size_bounds$min_size, max.size = LATR_size_bounds$max_size,
                          mat.size = TM.matdim){
    
    # Matrix size and size extensions (upper and lower integration limits)
    n <- mat.size
    L <- min.size + ext.lower
    U <- max.size + ext.upper
    
    # Bin size for n bins
    h <- (U - L)/n
    
    # Lower boundaries of bins 
    b <- L + c(0:n)*h
    
    # Bin midpoints
    y <- 0.5*(b[1:n] + b[2:(n + 1)])
    
    # Growth/Survival matrix
    Pmat <- t(outer(y, y, TM.growsurv, d = dens))*h 
    
    # Fertility/Recruiment matrix
    Fmat <- t(outer(y, y, TM.fertrecruit, d = dens))*h 
    
    # Put it all together
    IPMmat <- Pmat + Fmat
    
    # And transition matrix
    return(list(IPMmat = IPMmat, Fmat = Fmat, Pmat = Pmat, meshpts = y))}
  
  # Construct transition matrix for minimum weighted density (zero)
  TM <- TransMatrix(dens = 0)
  lambda.recruitment[i] <- lambda(TM$IPMmat)
  
  # Function to calculate the minimum wavespeed across a range of s
  Wavespeed <- function(n = TM.matdim){
    
    # Fit equation to convert volume to height for dispersal kernel use
    LATR_full %>%
      select(max.ht_t, volume_t) %>% 
      drop_na(max.ht_t, volume_t) %>% 
      rename("h" = max.ht_t, "v" = volume_t) %>% 
      arrange(v) %>% 
      nlsLM(h ~ A*v^(1/3),
            start = list(A = 0), data = .) %>% 
      coef() %>% 
      as.numeric() -> A
    
    # To see data, store data frame as test and then use plot(test$v, test$h)
    # To see fit line, store the model as fit and then use lines(sort(test$v), fitted(fit), col = "red")
    
    # Function converting volume to height
    vol.to.height <- function(v){
      h <- A*v^(1/3)
      return(h)}
    
    # Vector of n heights across which dispersal kernel will be evaluated
    z.list <- sapply(exp(TM$meshpts), vol.to.height)/100
    
    # List of simulated dispersal distances for each height
    r.list <- as.list(sapply(z.list[z.list >= 0.15], WALD.f.e.h))
    
    # Define modified bessel function for product of s and dispersal distance
    bessel <- function(r, t){
      return(besselI(t*r, 0))}
    
    # Function to evaluate MGF for each height; returns MGFs at each height for a specific value of s
    MGF.s <- function(s){
      
      # No dispersal when z < 0.15; use dirac delta function, with MGF of 1
      mgf.values.a <- rep(1, length(z.list[z.list < 0.15]))
      
      # For all other heights, evaluate bessel at each dispersal distance; find mean across all distances
      sapply(r.list[1:length(r.list)], bessel, t = s) %>% 
        sapply(., mean) -> mgf.values.b
      
      # Return MGF values for all heights
      mgf.values <- c(mgf.values.a, mgf.values.b)
      return(mgf.values)}
    
    # Set up range of s values over which to calculate wavespeeds
    s.seq <- c(0.0001, 0.0005, 0.001, 0.005, seq(0.01, 2, length.out = 96))
    
    # Apply MGF for each value of s
    mgf.over.s <- mapply(MGF.s, s = s.seq)
    
    # Define function to calculate wavespeed for each value of s
    ws.calc <- function(m, s){
      H <- TM$Fmat %*% diag(as.vector(m)) + TM$Pmat
      rho <- Re(eigen(H)$values[1])
      (1/s)*log(rho)}
    
    # Create empty vector on which to add wave speeds
    vec <- c()
    
    # Calculate wavespeed for each s and add it to the vector
    for(i in 1:100){
      val <- ws.calc(m = mgf.over.s[, i], s = s.seq[i])
      vec <- append(vec, val)}
    
    # Return vector of wavespeeds
    return(vec)}
  
  # Take the average of wave wavespeeds at this recruitment rate
  # Wavespeeds as function of s; growth as function of density
  hold.speed <- c()
  for(j in 1:3){
    c.values <- Wavespeed()
    hold.speed <- c(hold.speed, min(c.values))}
  wavespeed.recruitment[i] <- mean(hold.speed, na.rm = T)
  
  # Calculate minimum wavespeed
  # c.values <- Wavespeed()
  # wavespeed.recruitment[i] <- min(c.values)
  
  # Print wavespeeds
  print(i); print(recruitment[i]); print(wavespeed.recruitment[i])}

# Plot wavespeeds
xaxis <- 10^(-recruitment)
plot(log10(xaxis), wavespeed.recruitment, type = "b", pch = 16,
     ylab = "Wavespeed (m/yr)", xlab = "log(Recruitment probability)", cex.lab = 1.4)
plot(log10(10^(-recruitment)), wavespeed.recruitment, type = "b", pch = 16)
plot(lambda.recruitment, wavespeed.recruitment, type = "b", pch = 16)


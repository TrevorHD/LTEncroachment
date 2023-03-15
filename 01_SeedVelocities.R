##### Calculate terminal velocities -----------------------------------------------------------------------------------------------------------------

# Code author(s): Trevor

# Load data file with seed positions and times
file <- "https://github.com/TrevorHD/LTEncroachment/raw/master/Data/SD_Trials.csv"

# Create empty vector to populate with terminal velocities
tv.raw <- c()

# Define gravity (m/s)
g <- 9.81

# Go through velocities for each trial and find terminal velocities
for(i in 13:62){
  
  # Load each column of positions
  TX <- paste0("T", i, ".y")
  data <- read.csv(file) %>% 
    dplyr::select("Time", TX) %>% 
    rename(t = "Time", y = TX) %>% 
    na.omit()
  
  # Find initial position
  y.0 <- min(data$y)
  
  # Fit to equation of motion for quadratic drag using Levenberg-Marquardt algorithm
  nlslrc <- nlsLM(data$y ~ y.0 + ((v^2)/g)*log(cosh(g*t/v)), 
                  start = list(v = 2.5), data = data)
  
  # Extract terminal velocity
  v <- as.numeric(coef(nlslrc))
  
  # Add terminal velocity to vector
  tv.raw <- append(tv.raw, v)}
  
# Remove outliers; set index to [tv.raw < 5] to further restrict outliers
tv.raw <- tv.raw[tv.raw < 10]





##### Create terminal velocity distributions --------------------------------------------------------------------------------------------------------

# Empirical probability distribution
tv.PDF <- density(tv.raw, from = 0, to = 5, bw = 0.2)

# Fit lognormal probability distribution to data
tv.ln <- fitdistr(tv.raw, "lognormal")

# Fit gamma probability distribution to data
tv.gm <- fitdistr(tv.raw, "gamma")

# Store coefficients in one vector
tv.fits <- c()

  # Lognormal coefficients
  tv.fits[1] <- as.numeric(tv.ln[[1]][1])       # Mean log
  tv.fits[2] <- as.numeric(tv.ln[[1]][2])       # SD log

  # Gamma coefficients
  tv.fits[3] <- as.numeric(tv.gm[[1]][1])       # Shape
  tv.fits[4] <- as.numeric(tv.gm[[1]][2])       # Rate





##### Clean up --------------------------------------------------------------------------------------------------------------------------------------

# Clean up variables from global environment
remove(v, y.0, g, TX, i, data, nlslrc, tv.ln, tv.gm, file)


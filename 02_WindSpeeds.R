##### Initialise data -------------------------------------------------------------------------------------

# TM update 6 July 2022: grabbing 2015-19 met data from station 49 only (speeds things up)
# Get wind data from repo

# Create URL to read from; replace with local file path for faster performance
url <- "https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/Data/Weather/Sevilleta_LTER_Hourly_Meteorological_Data_2015_2019_station49.csv"
    
# Using read.csv is too slow; SQL will be much more efficient
# Read comma-separated data from text file, and select relevant columns
df <- read.csv.sql(url, sql = "SELECT Mean_WindSpeed, Wind_Dir FROM file", eol = "\n")
    
# Close SQL database connection
sqldf()
  
# Coerce data frame to vector
# Must use matrix as intermediate since it won't work directly for some reason
ws.raw <- as.vector(as.matrix(df$Mean_WindSpeed))
  
# Remove (negative) markers for missing values
# Also remove zero wind speeds, as we will assume no seed release occurs in absence of wind
ws.raw <- ws.raw[ws.raw > 0]





##### Create wind speed distributions ---------------------------------------------------------------------

# Create empirical PDF
ws.PDF <- density(ws.raw, from = 0, to = 15)

# Fit lognormal probability distribution to data
ws.ln <- fitdistr(ws.raw, "lognormal")

# Fit gamma probability distribution to data
ws.gm <- fitdistr(ws.raw, "gamma")

# Store coefficients in one vector
ws.fits <- c()
  
  # Lognormal coefficients
  ws.fits[1] <- as.numeric(ws.ln[[1]][1])       # Mean log
  ws.fits[2] <- as.numeric(ws.ln[[1]][2])       # SD log

  # Gamma coefficients
  ws.fits[3] <- as.numeric(ws.gm[[1]][1])       # Shape
  ws.fits[4] <- as.numeric(ws.gm[[1]][2])       # Rate





##### Clean up --------------------------------------------------------------------------------------------

# Clean up variables from global environment
remove(ws.ln, ws.gm, df, url)


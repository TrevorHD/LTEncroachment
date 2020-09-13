##### Initialise data -------------------------------------------------------------------------------------

# Create vector of year ranges used in file URLs
years <- c("1988-1995_0", "1996-2000_1", "2001-2005", "2006-2010")

# Get wind data from online
for(i in 1:4){
  
  # Create URL
  url <- paste0("http://sevlter.unm.edu/sites/default/files/data/sev1_meteorology_",
                years[i], ".txt")
    
  # Base R does not play nice with read.table and these files
  # SQL will be much more efficient
    
  # Read comma-separated data from text file; keep all columns
  # Can't use SELECT to extract only wind speeds since 3rd file has unnamed columns
  df <- read.csv.sql(url, sql = "SELECT * FROM file", sep = ",", skip = 1)
    
  # Close SQL database connection
  sqldf()
  
  # Use data from the Five Points site (ID = 49)
  df <- subset(df, df[1] == 49)
    
  # Only keep the mean hourly wind speeds (column 9)
  df <- df[9]
  
  # Coerce data frame to vector
  # Must use matrix as intermediate since it won't work directly for some reason
  df <- as.vector(as.matrix(df))
  
  # Assign name to vector
  assign(paste0("ws.", i), df)}
  
# Combine individual vectors
ws.raw <- c(ws.1, ws.2, ws.3, ws.4)
  
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
remove(ws.ln, ws.gm, df, i, url, years, ws.1, ws.2, ws.3, ws.4)


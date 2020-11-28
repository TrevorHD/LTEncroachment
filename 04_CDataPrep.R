##### Initialise Data -------------------------------------------------------------------------------------

# Read transect densities
CData.Transects <- read.csv("https://github.com/TrevorHD/LTEncroachment/raw/master/Data/LT_TransectData.csv")
names(CData.Transects)[1] <- "site"

# Read longitudinal demography data
CData.Demography <- read.csv("https://github.com/TrevorHD/LTEncroachment/raw/master/Data/LT_DemographyData.csv")
names(CData.Demography)[1] <- "site"

# Read transect lengths
CData.Lengths <- read.csv("https://github.com/TrevorHD/LTEncroachment/raw/master/Data/LT_TransectLengths.csv")

# Read transplant data
CData.Transplants <- read.csv("https://github.com/TrevorHD/LTEncroachment/raw/master/Data/LT_TransplantExp.csv")
names(CData.Transplants)[1] <- "site"





##### Calculate densities for each 5-m window -------------------------------------------------------------

# Calculate plant volume (elliptical cone)
vol <- function(h, w, p){(1/3)*pi*h*(sqrt((w*p))/2)^2}

# Apply calculations to transects data frame
CData.Transects$volume <- vol(h = CData.Transects$max.ht, w = CData.Transects$max.w,
                              p = CData.Transects$perp.w)
# Define study sites
site <- c("FPS", "SLP", "PDC", "MOD")

# Set up windows for each site
for(i in 1:4){
  
  # Create lists of 5-m windows
  exprs <- list(seq(0, CData.Lengths$Length_m[CData.Lengths$Site == site[i] & CData.Lengths$Transect == 1], 5),
                seq(0, CData.Lengths$Length_m[CData.Lengths$Site == site[i] & CData.Lengths$Transect == 2], 5),
                seq(0, CData.Lengths$Length_m[CData.Lengths$Site == site[i] & CData.Lengths$Transect == 3], 5))
  
  # Assign names to these window lists
  assign(paste0("Windows.", site[i]), exprs)}

# Prepare window lists for binding
for(i in 1:4){
  
  # Get data for window lists
  window.i <- paste0("Windows.", site[i]) %>% get()
  
  # Bind window lists for each transect into 3 columns
  exprs <- cbind(rep(site[i], length(c(window.i[[1]], window.i[[2]], window.i[[3]]))),
                 c(rep(1, length(window.i[[1]])), rep(2,length(window.i[[2]])), rep(3,length(window.i[[3]]))),
                 c(window.i[[1]], window.i[[2]], window.i[[3]]))
  
  # Assign names to bind-ready window lists
  assign(paste0("WBind.", site[i]), exprs)}
  
# Bind all window list columns into a data frame
Windows <- data.frame(rbind(WBind.FPS, WBind.SLP, WBind.PDC, WBind.MOD))

# Name Windows and convert characters to integers
names(Windows) <- c("site", "transect", "window")
Windows$window <- as.integer(as.character(Windows$window))
Windows$transect <- as.integer(as.character(Windows$transect))

# Loop over data frame to add densities by window
for(i in 1:nrow(Windows)){
  
  # Density (# of individuals per window)
  plants <- na.omit(CData.Transects$window[CData.Transects$window == Windows$window[i] &
                                           CData.Transects$transect == Windows$transect[i] &
                                           CData.Transects$site == Windows$site[i]])
  
  # Add density to data frame
  Windows$density[i] <- length(plants)
  
  # Add weighted density (sum of log volume per window) to data frame
  Windows$weighted.dens[i] <- sum(log(CData.Transects$volume[CData.Transects$window == Windows$window[i] &
                                                             CData.Transects$transect == Windows$transect[i] &
                                                             CData.Transects$site == Windows$site[i]]), na.rm = T)}

# Final result is a list of 5-m windows and their densities for each transect at each site





##### Demography data QA/QC -------------------------------------------------------------------------------

# Remove extra site level, since it is just a duplicate with a space after it
CData.Demography$site[CData.Demography$site == "FPS "] <- "FPS"

# Rename "reproductive_fraction" to "reproductive_fraction_t1"
# This should be done because it is recorded at the end of the transition year
CData.Demography %>% 
  rename(reproductive_fraction_t1 = reproductive_fraction) %>% 
  mutate(transect = as.integer(transect),
         site = droplevels(site)) -> CData.Demography

# Blanks in survival should be converted to NA
CData.Demography$survival_t1[CData.Demography$survival_t1 == ""] <- NA

# There is a "." for one survival entry; original data indicate this plant was dead
CData.Demography[which(CData.Demography$survival_t1 == "."), "survival_t1"] <- 0

# Clean up survival factor
CData.Demography$survival_t1 <- as.integer(as.character(droplevels.factor(CData.Demography$survival_t1)))

# Fix entries whose reproductive fraction should be 1
CData.Demography$reproductive_fraction_t[CData.Demography$reproductive_fraction_t > 1] <- 1
CData.Demography$reproductive_fraction_t1[CData.Demography$reproductive_fraction_t1 > 1] <- 1

# Find implausible/incorrect size transitions and remove these observations
# Check height changes below the 2.5th and about the 97.5th percentile
CData.Demography %>% mutate(height_change = log(max.ht_t1/max.ht_t)) %>% 
  filter(height_change > quantile(height_change, 0.975, na.rm = T) | height_change < quantile(height_change, 0.025, na.rm = T)) %>% 
  select(site, transect, designated.window, plant, year_t, max.ht_t, 
         max.w_t, perp.w_t, max.ht_t1, max.w_t1, perp.w_t1) %>% 
  arrange(site, transect, designated.window, plant, year_t)

# Go through line by line and pull out the plant-years are problems and should be dropped
# These plants have inexplicable/unbelievable size changes that cannot be verified or corrected with raw data
# Note: height changes (esp. reductions) are sensitive to single branches dying back, which we count as "real"
problems <- data.frame(site = factor(NA, levels = c("FPS", "MOD", "PDC", "SLP")), 
                       transect = NA, designated.window = NA, plant = NA, year_t = NA)

problems[1,] <- c("FPS", 1, 150, 2, 2013)

# FPS 2-0 plants 4,5,6 were likely in a bunch and mixed up
problems[2,] <- c("FPS", 2, 0, 4, 2014)
problems[3,] <- c("FPS", 2, 0, 5, 2014)
problems[4,] <- c("FPS", 2, 0, 5, 2015)

# FPS 2-150-12 was recorded using data from 12s instead of 12
problems[5,] <- c("FPS", 2, 150, 12, 2016)

# FPS 3-100-7 has a data entry problem (confirmed after checking data sheets)
CData.Demography$max.ht_t1[CData.Demography$site == "FPS" & CData.Demography$transect == 3 & CData.Demography$designated.window == 100 & CData.Demography$plant == 7 & CData.Demography$year_t == 2015] <- 38
CData.Demography$max.ht_t[CData.Demography$site == "FPS" & CData.Demography$transect == 3 & CData.Demography$designated.window == 100 & CData.Demography$plant == 7 & CData.Demography$year_t == 2016] <- 38

# FPS 3-100 plants 4 and 6 seem like they were mixed up or perhaps growing on top of each other
# Dropping these for all years
problems[6,] <- c("FPS", 3, 100, 4, 2013)
problems[7,] <- c("FPS", 3, 100, 4, 2014)
problems[8,] <- c("FPS", 3, 100, 4, 2015)
problems[9,] <- c("FPS", 3, 100, 4, 2016)
problems[10,] <- c("FPS", 3, 100, 6, 2013)
problems[11,] <- c("FPS", 3, 100, 6, 2014)
problems[12,] <- c("FPS", 3, 100, 6, 2015)
problems[13,] <- c("FPS", 3, 100, 6, 2016)

# This one can't be right; 2013 measurement of 126 should be 26 but dropping due to uncertainty
problems[14,] <- c("FPS", 3, 100, 9, 2013)

# FPS 3-200-8 has some strange size changes, might have gotten mixed up with another plant
# Dropping these observations
problems[15,] <- c("FPS", 3, 200, 8, 2013)
problems[16,] <- c("FPS", 3, 200, 8, 2014)
problems[17,] <- c("FPS", 3, 200, 8, 2015)
problems[18,] <- c("FPS", 3, 200, 8, 2016)

# FPS-3-500 1 and 2 likely got mixed up in 2013-2014
problems[19,] <- c("FPS", 3, 500, 1, 2013)
problems[20,] <- c("FPS", 3, 500, 1, 2014)
problems[21,] <- c("FPS", 3, 500, 2, 2013)
problems[22,] <- c("FPS", 3, 500, 2, 2014)

# This change is not believable
problems[23,] <- c("MOD", 1, 150, 1, 2016)

# This change is not believable; no apparent problem in raw data
problems[24,] <- c("MOD", 2, 50, 3, 2014)

# These dimensions are not believale... someone probably said "60" and someone wrote "16"; dropping
problems[25,] <- c("MOD", 3, 0, 9, 2014)
problems[26,] <- c("MOD", 3, 0 ,9, 2015)

# This one is not believable; can't find the 2015 data to verify
problems[27,] <- c("PDC", 1, 200, 2, 2014)
problems[28,] <- c("PDC", 1, 200, 2, 2015)

# This one has a note that it was untagged and we thought it was right... it wasn't
problems[29,] <- c("SLP", 3, 100, 8, 2016)

# Ensure data types are correct before removing problem entries
problems %>% 
  mutate(transect = as.integer(transect),
         designated.window = as.integer(designated.window),
         plant = as.factor(plant),
         year_t = as.integer(year_t)) -> problems

# Remove problematic entries and create new CData set
CData <- anti_join(CData.Demography,problems, by = c("site", "transect", "designated.window", "plant", "year_t"))

# The new set should have 29 fewer rows; compare nrow(problems) to nrow(CData)-nrow(CData.Demography)





##### Fix window-related issues ---------------------------------------------------------------------------

# Identify entries that lack an actual window (5-m resolution)
subset(CData, is.na(actual.window)) %>% 
  select("site", "transect", "designated.window", "plant", "year_t") %>% 
  unique() -> CData.MissingWindow

# Check why actual Windows are not present
# See "Missing Windows" file for more information

# check whether any designated windows are missing
CData %>% filter(is.na(designated.window))

# Two SLP recruits are missing designated windows, but the notes place them under specified plants
CData$designated.window[CData$site == "SLP"&CData$transect == 3&CData$plant == "10s"&CData$year_t==2016] <- 
  CData$designated.window[CData$site == "SLP"&CData$transect == 3 & CData$designated.window == 150 & CData$plant == "7"][1]
CData$designated.window[CData$site == "SLP"&CData$transect == 3&CData$plant == "11s"&CData$year_t==2016] <- 
  CData$designated.window[CData$site == "SLP"&CData$transect == 3 & CData$designated.window == 150 & CData$plant == "5"][1]

# Use designated window for sites that don't have an actual window
CData$actual.window[is.na(CData$actual.window)] <- 
  CData$designated.window[is.na(CData$actual.window)]





##### Remove unnecessary columns and merge densities with demography data ---------------------------------

# Keep only useful columns
select(CData, "site", "transect", "designated.window", "actual.window", "plant", "year_t",
       "max.ht_t", "max.w_t", "perp.w_t", "flowers_t", "fruits_t", "reproductive_fraction_t", 
       "year_t1", "new.plant_t1", "seedling_t1", "survival_t1", "max.ht_t1", "max.w_t1", 
       "perp.w_t1", "flowers_t1", "fruits_t1", "reproductive_fraction_t1") %>% 

# Merge with demography data
merge(Windows, 
      by.x = c("site", "transect", "actual.window"),
      by.y = c("site", "transect", "window")) -> CData

# Final result is each plant and its demography, marked with its 5-m window weighted density





##### Calculate quantities that will be used in analyses --------------------------------------------------

# Add additional columns to data, starting with log initial volume of plant before year has elapsed
CData %>%
  mutate("volume_t" = vol(h = max.ht_t, w = max.w_t, p = perp.w_t),
       
         # Final log conical volume of plant after year of growth
         "volume_t1" = vol(h = max.ht_t1, w = max.w_t1, p = perp.w_t1),
       
         # Logarithmic annual growth ratio
         "logGR" = log(volume_t1) - log(volume_t),

         # Total number of fruits on a given plant before year has elapsed
         "total.fruits_t" = fruits_t*(1/reproductive_fraction_t),
       
         # Total number of flowers on a given plant before year has elapsed
         "total.flowers_t" = flowers_t*(1/reproductive_fraction_t),
       
         # Total number of fruits on a given plant after year has elapsed
         "total.fruits_t1" = fruits_t1*(1/reproductive_fraction_t1),
       
         # Total number of flowers on a given plant after year has elapsed
         "total.flowers_t1" = flowers_t1*(1/reproductive_fraction_t1),
       
         # Total number of reproductive structures on a given plant before year has elapsed
         "total.reproduction_t" = total.fruits_t + total.flowers_t,
         
         # Total number of reproductive structures on a given plant after year has elapsed
         "total.reproduction_t1" = total.fruits_t1 + total.flowers_t1,
       
         # Boolean stating whether or not the plant flowered in t1
         "did.flower" = total.reproduction_t1 > 0,
       
         # Variable indicating these are not transplants
         "transplant" = FALSE) -> CData





##### Quantify total recruitment  -------------------------------------------------------------------------

# Create data frame of recruits; will be used for later calculations, but not here
# We'll need to report this criterion for designating a recruit (log vol < 8)
CData.Recruits <- filter(CData, new.plant_t1 == 1 | seedling_t1 == 1, log(volume_t1) < 8)

# The code below will likely be deprecated
# Calculate total number of seedlings (recruits) in a single year for each 5-m window
# for(i in 1:nrow(CData)){
#  CData$recruits.1y[i] <- sum(CData$new.plant_t1[CData$actual.window == CData$actual.window[i] &
#                                                 CData$transect == CData$transect[i] &
#                                                 CData$site == CData$site[i] &
#                                                 CData$year_t1 == CData$year_t1[i]], na.rm = T)}

# Calculate total number of seedlings (recruits) across all years for each 5-m window
# for(i in 1:nrow(CData)){
#  CData$recruits.4y[i] <- sum(CData$new.plant_t1[CData$actual.window == CData$actual.window[i] &
#                                                 CData$transect == CData$transect[i] &
#                                                 CData$site == CData$site[i]], na.rm = T)}

# Write transect data to use elsewhere for estimating recruitment per seed; merge Windows and Cdata.Transects
# Create data frame with plant sizes and window densities for all transects and windows
left_join(CData.Transects,Windows,by=c("site","transect","window")) %>% 
  select(site,transect,window,volume,weighted.dens) -> Cdata.Transects.Windows





##### Tidy up transplant data for survival analysis -------------------------------------------------------

# Locations are multiples of 2.5 m, but our density data are in 5-m windows
# Therefore, we will round locations up to nearest 5-m window
for(i in 1:nrow(CData.Transplants)){
  if(CData.Transplants$plot_location[i] %% 2.5 == 0){
    CData.Transplants$plot_location[i] <- CData.Transplants$plot_location[i] + 2.5}}

# Calculate conical volume; use initial volume since most plants die after a year
# Add variable indicating which plants are transplants
CData.Transplants %>% 
  mutate("volume_t" = vol(h = max.ht_t, w = max.w_t, p = perp.w_t),
         "transplant" = TRUE) %>% 
    
# Rename "plot_location" to "actual.window"
rename("actual.window" = "plot_location") %>% 
  
# Merge with demography data
merge(Windows, 
      by.x = c("site", "transect", "actual.window"),
      by.y = c("site", "transect", "window")) -> CData.Transplants





##### Clean up global environment -------------------------------------------------------------------------

# Remove variables that will no longer be used
remove(site, plants, i, Windows.FPS, Windows.MOD, Windows.PDC, Windows.SLP,
       window.i, WBind.FPS, WBind.MOD, WBind.PDC, WBind.SLP, exprs)


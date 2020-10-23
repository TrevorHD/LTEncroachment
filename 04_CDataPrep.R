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
  
  # Add weighted density (amount of volume per window) to data frame
  Windows$weighted.dens[i] <- sum(CData.Transects$volume[CData.Transects$window == Windows$window[i] &
                                                         CData.Transects$transect == Windows$transect[i] &
                                                         CData.Transects$site == Windows$site[i]], na.rm = T)}

# Final result is a list of 5-m windows and their densities for each transect at each site




# Data QA/QC --------------------------------------------------------------

# Fix entry where "." in survival should be a "0"; drop "." factor
CData[which(CData$survival_t1 == "."), "survival_t1"] <- 0
CData$survival_t1 <- droplevels.factor(CData$survival_t1)

# Remove other invalid entries
# FPS 1-0-4, t1=2016, CDataV2 entry 651 (2 plants were accidentally measured as 1)
# FPS 1-0-6, t1=2016, CDataV2 entry 653 (2 plants were accidentally measured as 1)
# FPS 2-150-12, t1=2017, CDataV2 entry 1245 (large loss likely indicates measurement error)
# SLP 3-100-8, t1=2017, CDataV2 entry 1443 (large loss likely indicates measurement error)
# PDC 3-0-8, t1=2017, CDataV2 entry 1519 (large gain likely indicates measurement error)
# MOD 1-150-1, t1=2017, CDataV2 entry 1574 (large loss likely indicates measurement error)
# MOD 2-50-3, t1=2017, CDataV2 entry 1593 (large loss likely indicates measurement error)
# MOD 3-200-1, t1=2017, CDataV2 entry 1633 (large loss likely indicates measurement error)
CData <- CData[-c(1610, 1570, 1552, 1497, 1421, 1225, 643, 641), ]

# Note that entry position =! row number (due to previous deletions)
# Might change if original spreadsheet or previous code is modified!

## Find implausible / incorrect size transitions and remove these observations
## check height changes below the 2.5th and about the 97.5th percentile
CData %>% mutate(height_change = log(max.ht_t1/max.ht_t)) %>% 
  filter(height_change > quantile(height_change,0.975,na.rm=T) | height_change < quantile(height_change,0.025,na.rm=T)) %>% 
  select(site,transect,designated.window,plant,year_t,max.ht_t,max.ht_t1)
drop_row <-c()
#     site transect designated.window plant year_t max.ht_t max.ht_t1
# 17  FPS        2               150    12   2016     47.0       7.0 
## checked data sheets and this appears to be a field error -- DROP
drop_row <- c(drop_row,17)
# 29  FPS        3               100     9   2013    126.0      42.0
## there are field notes that FPS-3-100 8,9, and 10 may be connected, so I am going to drop these since there was likely confusion
drop_row <- c(drop_row,29)
# 32  FPS        3               100     4   2014     58.0     109.0
## field data were entered correctly but this size change can't be right, dropping
drop_row <- c(drop_row,32)
# 34  FPS        3               100     7   2016      8.0      36.0
# 35  FPS        3               100     7   2015     35.0       8.0
## data entry error: the 8's should be 38's
CData$max.ht_t[CData$site=="FPS"&CData$transect==3&CData$designated.window==100&CData$plant==7&CData$year_t==2016]<-
CData$max.ht_t1[CData$site=="FPS"&CData$transect==3&CData$designated.window==100&CData$plant==7&CData$year_t==2015]<-38
# 45  FPS        3               500     1   2014     95.0      45.0
# 46  FPS        3               500     1   2013     38.0      95.0
# 47  FPS        3               500     2   2014     37.0     105.0
# 48  FPS        3               500     2   2013     83.0      37.0
## It looks like FPS-3-500 1 and 2 were mixed up for field data collection in 2014 only - swap these observations

# 51  MOD        3                 0     9   2015     16.0      67.0
# 52  MOD        3                 0     9   2014     60.0      16.0
# 58  PDC        1               200     2   2014     30.0      21.0
# 59  PDC        1               200     2   2015     21.0      35.0
# 80  SLP        3               100     8   2016     83.0      41.0
##### Fix window-related issues ---------------------------------------------------------------------------

# Identify entries that lack an actual window (5-m resolution)
subset(CData, is.na(actual.window)) %>% 
  select("site", "transect", "designated.window", "plant") %>% 
  unique() -> CData.MissingWindow

# Check why actual Windows are not present
# See "Missing Windows" file for more information

# Eliminate recruits missing a designated window (50-m resolution)
CData <- CData[complete.cases(CData[ , 3]), ]

# Use designated window for sites that don't have an actual window
CData$actual.window[is.na(CData$actual.window)] <- 
  CData$designated.window[is.na(CData$actual.window)]





##### Remove unnecessary columns and merge densities with demography data ---------------------------------

# Keep only useful columns
select(CData, "site", "transect", "designated.window", "actual.window", "plant", "year_t",
       "max.ht_t", "max.w_t", "perp.w_t", "flowers_t", "fruits_t", "reproductive_fraction_t", 
       "year_t1", "new.plant_t1", "seedling_t1", "survival_t1", "max.ht_t1", "max.w_t1", 
       "perp.w_t1", "flowers_t1", "fruits_t1", "reproductive_fraction") %>% 

# Merge with demography data
merge(Windows, 
      by.x = c("site", "transect", "actual.window"),
      by.y = c("site", "transect", "window")) -> CData

# Final result is each plant and its demography, marked with its 5-m window weighted density





##### Calculate quantities that will be used in analyses --------------------------------------------------

CData %>%
  
  # Rename "reproductive_fraction" to "reproductive_fraction_t1" to avoid confusion
  rename("reproductive_fraction_t1" = "reproductive_fraction") %>% 

  # Add additional columns to data, starting with log initial volume of plant before year has elapsed
  mutate("volume_t" = log(vol(h = max.ht_t, w = max.w_t, p = perp.w_t)),
       
         # Final log conical volume of plant after year of growth
         "volume_t1" = log(vol(h = max.ht_t1, w = max.w_t1, p = perp.w_t1)),
       
         # Logarithmic annual growth ratio
         "logGR" = ifelse(is.nan(volume_t1 - volume_t) == TRUE, NA, log(volume_t1) - log(volume_t)),

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
       
         # Log of total t1 reproduction
         # Plants with no reproduction get -Inf, so we simply treat this as a NA
         "logTR1" = ifelse(log(total.reproduction_t1) == -Inf, NA, log(total.reproduction_t1)),
       
         # Boolean stating whether or not the plant flowered in t1
         "did.flower" = ifelse(total.fruits_t1 > 0 | total.flowers_t1 > 0, 1, 0),
       
         # Variable indicating these are not transplants
         "transplant" = FALSE) -> CData

## write out CData for Tom to use elsewhere
#write.csv(CData,"C:/Users/tm9/Desktop/git local/IPM_size_transitions/creosote/CData.csv")



##### Quantify total recruitment  -------------------------------------------------------------------------

# Create data frame of recruits; will be used for later calculations, but not here
## tom: we'll need to report this criterion for designating a recruit (log vol < 8)
## also, 8 is pretty big in log(vol). did you mean log(8)?
CData.Recruits <- filter(CData, new.plant_t1 == 1 | seedling_t1 == 1, log(volume_t1) < 8)

# Calculate total number of seedlings (recruits) in a single year for each 5-m window
for(i in 1:nrow(CData)){
  CData$recruits.1y[i] <- sum(CData$new.plant_t1[CData$actual.window == CData$actual.window[i] &
                                                 CData$transect == CData$transect[i] &
                                                 CData$site == CData$site[i] &
                                                 CData$year_t1 == CData$year_t1[i]], na.rm = T)}

# Calculate total number of seedlings (recruits) across all years for each 5-m window
for(i in 1:nrow(CData)){
  CData$recruits.4y[i] <- sum(CData$new.plant_t1[CData$actual.window == CData$actual.window[i] &
                                                 CData$transect == CData$transect[i] &
                                                 CData$site == CData$site[i]], na.rm = T)}

# Per-seed recruitment rates are calculated in 05_CDataAnalysis





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


## write out transplant data for Tom to use elsewhere
#write.csv(CData.Transplants,"C:/Users/tm9/Desktop/git local/IPM_size_transitions/creosote/CData.Transplants.csv")


##### Clean up global environment -------------------------------------------------------------------------

# Remove variables that will no longer be used
remove(site, plants, i, Windows.FPS, Windows.MOD, Windows.PDC, Windows.SLP,
       window.i, WBind.FPS, WBind.MOD, WBind.PDC, WBind.SLP, exprs)


##### Load and plot monsoon precipitation -----------------------------------------------------------------------------------------------------------

# Code author(s): Tom

# Data from Kris Hall: here is annual July-September monsoon precipitation by calendar year (not by water year)
# for each of the met stations, with units in millimeters. 

# Load data
monsoon <- read.csv("https://github.com/TrevorHD/LTEncroachment/raw/master/Data/Weather/SEVMonsoonPrecip.csv")
str(monsoon)

# Calculate seedlings per area by getting transect lengths
transect_lengths <- read.csv("https://github.com/TrevorHD/LTEncroachment/raw/master/Data/LT_TransectLengths.csv") %>% 
  mutate(area = Length_m*2)
CData %>% 
  group_by(site, transect,year_t1) %>% 
  summarise(seedlings = sum(new.plant_t1 == 1 & seedling_t1 == 1, na.rm = T)) %>% 
  rename(Site = site, Transect = transect) %>% 
  left_join(., transect_lengths, by = c("Site", "Transect"))%>% 
  group_by(year_t1) %>% 
  summarise(tot_recruits = sum(seedlings),
            tot_area = sum(area)) %>% 
  mutate(recruits_per_area = tot_recruits/tot_area) -> recruits

# Prepare plotting to PDF
pdf("Manuscript/Figures/monsoon_seedlings.pdf", useDingbats = F, height = 4, width = 8)
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))

# Generate plots
plot(monsoon$Year[monsoon$StationID == 49],
     monsoon$Precipitation[monsoon$StationID == 49], type = "b",
     ylab = "Monsoon (July-September) precip. (mm)", xlab = "Year")
points(2013:2017, monsoon$Precipitation[monsoon$StationID == 49 & monsoon$Year %in% 2013:2017], pch = 16)
title("A", adj = 0, font = 3)
plot(monsoon$Precipitation[monsoon$StationID == 49 & monsoon$Year %in% 2013:2016],
     recruits$recruits_per_area[recruits$year_t1 %in% 2014:2017],
     xlab = "Monsoon precip. (mm)", ylab = "Recruit density per m2", cex = 2, pch = 16)
title("B", adj = 0, font = 3)

# Turn off graphics device
dev.off()


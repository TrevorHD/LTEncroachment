##### Get monsoon precipitation data ----------------------------------------------------------------------------------------------------------------

# Code author(s): Tom

# Data from Kris Hall: here is annual July-September monsoon precipitation by calendar year (not by water year)
# for each of the met stations, with units in millimeters.

# Load packages
library(tidyverse)
library(lubridate)

# Set folder where the output file will be written to
path_for_data_export <- "~/Documents/SEV/Requests/Miller_Tom_monsoon_precip_20220611/"

# Get meteorological data from EDI
m8894 <- read_csv("https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-sev.1.15&entityid=87ff32c8e179e69743c6514cbbc8b31c",
                  guess_max = 1000000)
m9599 <- read_csv("https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-sev.1.15&entityid=dbe548a38a8dc4c1c829657cee190c4d",
                  guess_max = 1000000)
m0004 <- read_csv("https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-sev.1.15&entityid=371109f8068b35cf65edc8ba4237c8bd",
                  guess_max = 1000000)
m0509 <- read_csv("https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-sev.1.15&entityid=e326dbe48c0cdc5b91496a469a50e36d",
                  guess_max = 1000000)
m1014 <- read_csv("https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-sev.1.15&entityid=011fd6eb9726321cace6c72b50cb8056",
                  guess_max = 1000000)
m1519 <- read_csv("https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-sev.1.15&entityid=76922a0b041ac5ab05be6132ff7f90d7",
                  guess_max = 1000000)
m2021 <- read_csv("https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-sev.1.15&entityid=bd8e2163c866a22a72136db5b52dd3b9",
                  guess_max = 1000000)

# Combine data, view, and tidy
m <- rbind(m8894, m9599, m0004, m0509, m1014, m1519, m2021)
glimpse(m)
m <- m %>% 
  mutate(StationID = as.factor(StationID)) %>% 
  select(StationID, Year, Month, Precipitation) %>% 
  filter(Month %in% c(7, 8, 9)) %>% 
  select(-Month)

# Calculate monsoon precip (July-September) by calendar year
ppt <- m %>% 
  group_by(StationID, Year) %>% 
  summarize(Precipitation = round(sum(Precipitation, na.rm = TRUE), 3)) %>% 
  filter(StationID != "46" & StationID != "47")

# Plot precipitation data
ggplot(ppt, aes(x = Year, y = Precipitation, color = StationID)) +
  geom_line() +
  facet_wrap(~ StationID)

# Write to file
write_csv(ppt, paste0(path_for_data_export, "SEV_monsoon_precipitation.csv"))


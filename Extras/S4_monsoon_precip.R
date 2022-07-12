## read in monsoon data
## From kris hall: Here is annual July-September monsoon precipitation by calendar year (not by water year) for each of the met stations. Units are millimeters. 

monsoon <- read.csv("Data/Weather/SEV_monsoon_precipitation.csv")
str(monsoon)

plot(monsoon$Year[monsoon$StationID==49],
     monsoon$Precipitation[monsoon$StationID==49],type="b")
points(2013:2016,
       monsoon$Precipitation[monsoon$StationID==49 & monsoon$Year%in%2013:2016],
       pch=16)

CData %>% 
  filter(new.plant_t1==1 & seedling_t1==1) %>% 
  group_by(year_t1) %>% 
  summarise(n())

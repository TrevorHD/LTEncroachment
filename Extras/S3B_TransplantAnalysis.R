## This is Tom's analysis of the transplant data

## create a new data frame from the existing transplant data
## for analysis that will appear in an appendix
transplants.app <- CData.Transplants %>%
  # # Add additional columns to data, starting with standardised weighted density
  mutate("d.stand" = (weighted.dens - mean(weighted.dens, na.rm = TRUE)) / sd(weighted.dens, na.rm = TRUE),
         # Total number of patches that contain any type of grass
         "total.grass" = num_black_gramma_t + num_blue_gramma_t + num_other_grass_t,
         # Total number of patches that contain any type of plant (including shrubs)
         "total.plant" = num_shrub_t + total.grass + num_other_t,
         # Unique transect identifier
         "unique.transect" = interaction(transect, site),
         "volume_t1" = log(vol(h = max.ht_t1, w = max.w_t1, p = perp.w_t1)),
         "logGR" = volume_t1 - volume_t)

## Density here is measured at two scales with two different metrics
## See how tightly those scales are correlated by summarising over plots

## First, some bookkeeping. There should be four reps of each plot number on each transect
table(transplants.app$unique.transect,transplants.app$plot)
## Why are some transects missing some plot numbers?

transplant.plots <- transplants.app %>% 
  group_by()
  summarise()
  

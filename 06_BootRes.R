##### Get sample of seed terminal velocities  -------------------------------------------------------------------------------------------------------

# Sample percentage of empirical seed terminal velocity distribution and create PDF
if(boot.on == TRUE){
  set.seed(seeds[i, 2])
  boot.tv.raw <- sample(tv.raw, round(length(tv.raw)*boot.prop), replace = TRUE)
  boot.tv.PDF <- density(boot.tv.raw, from = 0, to = 5, bw = 0.2)}

# Use full terminal velocity data if bootstrapping is turned off
if(boot.on == FALSE){
  boot.tv.raw <- tv.raw
  boot.tv.PDF <- tv.PDF}





##### Get sample of wind speeds  --------------------------------------------------------------------------------------------------------------------

# Sample percentage of empirical wind speed distribution and create PDF
if(boot.on == TRUE){
  set.seed(seeds[i, 2])
  boot.ws.raw <- sample(ws.raw, round(length(ws.raw)*boot.prop), replace = TRUE)
  boot.ws.PDF <- density(boot.ws.raw, from = 0, to = 15)}

# Use full wind speed data if bootstrapping is turned off
if(boot.on == FALSE){
  boot.ws.raw <- ws.raw
  boot.ws.PDF <- ws.PDF}





##### Get sample from demographic data  -------------------------------------------------------------------------------------------------------------

# Sample percentage of demographic data
if(boot.on == TRUE){
  LATR_full <- CData %>% 
    mutate(unique.transect = interaction(transect, site))
  set.seed(seeds[i, 2])
  LATR_full <- LATR_full[sample(1:nrow(LATR_full), round(nrow(LATR_full)*boot.prop), replace = FALSE), ]}





##### Get sample of recruits  -----------------------------------------------------------------------------------------------------------------------

# Sample percentages of recruits to calculate mean and sd recruit size
if(boot.on == TRUE){
  set.seed(seeds[i, 2])
  boot.LATR_recruit_size <- LATR_recruit_size[sample(1:nrow(LATR_recruit_size), round(nrow(LATR_recruit_size)*boot.prop), 
                                              replace = FALSE), ]
  boot.LATR_recruit_size <- data.frame(recruit_mean = mean(boot.LATR_recruit_size$log_volume),
                                       recruit_sd = sd(boot.LATR_recruit_size$log_volume))}

# Use full recruit data if bootstrapping is turned off
if(boot.on == FALSE){
  boot.LATR_recruit_size <- data.frame(recruit_mean = mean(LATR_recruit_size$log_volume),
                                       recruit_sd = sd(LATR_recruit_size$log_volume))}


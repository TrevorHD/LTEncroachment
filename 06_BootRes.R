##### Set up resampling for wind speed, terminal velocitiy, and demography  -------------------------------

# Sample percentage of empirical wind speed distribution and create PDF
boot.tv.raw <- sample(tv.raw, round(length(tv.raw)*boot.prop), replace = TRUE)
boot.tv.PDF <- density(boot.tv.raw, from = 0, to = 5, bw = 0.2)





##### Get sample of wind speeds  --------------------------------------------------------------------------

# Sample percentage of empirical wind speed distribution and create PDF
boot.ws.raw <- sample(ws.raw, round(length(ws.raw)*boot.prop), replace = TRUE)
boot.ws.PDF <- density(boot.ws.raw, from = 0, to = 15)





##### Get sample of recruits  -----------------------------------------------------------------------------

# Sample individuals (sans transplants) after calculating weighted density
boot.LATR_full <- LATR_full[sample(1:nrow(LATR_full), round(nrow(LATR_full)*boot.prop), replace = FALSE), ]

# Filter out seedlings and get their sizes
LATR_recruit_size <- boot.LATR_full %>% 
  filter(seedling_t1 == 1) %>% 
  mutate(log_volume = log(volume_t1))

# Plot distribution of recruit sizes
hist(LATR_recruit_size$log_volume)

# Create df of recruit sizes
LATR_recruit_size <- data.frame(recruit_mean = mean(LATR_recruit_size$log_volume),
                                recruit_sd = sd(LATR_recruit_size$log_volume))


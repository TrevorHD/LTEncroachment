##### Set up resampling for wind speed, terminal velocitiy, and demography  -------------------------------

# Sample percentage of empirical wind speed distribution and create PDF
boot.tv.raw <- sample(tv.raw, round(length(tv.raw)*boot.prop), replace = TRUE)
boot.tv.PDF <- density(boot.tv.raw, from = 0, to = 5, bw = 0.2)

# Note: in scripts 03-07, replace "tv.raw" with "boot.tv.raw" and "tv.PDF" with "boot.tv.PDF"





##### Get sample of wind speeds  --------------------------------------------------------------------------

# Sample percentage of empirical wind speed distribution and create PDF
boot.ws.raw <- sample(ws.raw, round(length(ws.raw)*boot.prop), replace = TRUE)
boot.ws.PDF <- density(boot.ws.raw, from = 0, to = 15)

# Note: in scripts 03-07, replace "ws.raw" with "boot.ws.raw" and "ws.PDF" with "boot.ws.PDF"





##### Get sample of individuals  --------------------------------------------------------------------------

# Do we bootstrap individuals before or after we've calculated weighted density?
# Do we bootstrap transplants as well?

# Will leave this section blank for now
##### Set up resampling for wind speed, terminal velocitiy, and demography  -------------------------------

# Sample percentage of empirical wind speed distribution and create PDF
boot.tv.raw <- sample(tv.raw, round(length(tv.raw)*boot.prop), replace = TRUE)
boot.tv.PDF <- density(boot.tv.raw, from = 0, to = 5, bw = 0.2)





##### Get sample of wind speeds  --------------------------------------------------------------------------

# Sample percentage of empirical wind speed distribution and create PDF
boot.ws.raw <- sample(ws.raw, round(length(ws.raw)*boot.prop), replace = TRUE)
boot.ws.PDF <- density(boot.ws.raw, from = 0, to = 15)





##### Get sample of individuals  --------------------------------------------------------------------------

# Sample individuals (sans transplants) after calculating weighted density
boot.CData <- CData[sample(1:nrow(CData), round(nrow(CData)*boot.prop), replace = FALSE), ]

# Sample transplants
boot.CData.Transplants <- CData.Transplants[sample(1:nrow(CData.Transplants), round(nrow(CData.Transplants)*boot.prop), 
                                            replace = FALSE), ]


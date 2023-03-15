##### How are recruits distributed across size classes? ---------------------------------------------------------------------------------------------

# Code author(s): Trevor

# Plot distribution of recruits across volume
CData.Recruits %>%
  select(volume_t1, seedling_t1) %>% 
  drop_na(seedling_t1) %>% 
  mutate(seedling_t1 = factor(seedling_t1, levels = c(0, 1))) %>% 
  ggplot() + aes(x = volume_t1, fill = seedling_t1) +
  geom_histogram(binwidth = 1, center = -0.5, col = "black") +
  scale_y_continuous(breaks = seq(0, 18, 2)) +
  scale_x_continuous(breaks = seq(0, 12, 2)) +
  scale_fill_manual(values = c("0" = "firebrick2", "1" = "dodgerblue1"),
                    breaks = c("0", "1"))

# new.plant_t1 is a more inclusive category for recruits than seedling_t1 is
# All seedlings are new plants, but not vice-versa





##### How are recruits distributed across densities? ------------------------------------------------------------------------------------------------

# Get only densities at which recruits are found
CData.Recruits %>% 
  select(weighted.dens) %>% 
  mutate(Type = "Recruit Densities") %>% 
  rename("Density" = "weighted.dens") -> RD

# Combine transect data with window data
TWD <- merge(CData.Transects, Windows, by = c("site", "transect", "window"))

# Number of different densities is far greater than number of recruits
# Get around this by bootstrapping, drawing 58 densities so number is equal to number of recruits
TWD.BS <- as.data.frame(replicate(1000000, sample(TWD$weighted.dens, nrow(RD), 
                                                  replace = FALSE), simplify = TRUE))

# Sort and find mean density for each of 58 entry positions
# Then bind with recruit densities and plot
TWD.BS %>% sapply(function(x) sort(x)) %>% 
  rowMeans() %>% as.data.frame() %>% 
  mutate(Type = "All Densities") %>% 
  rename("Density" = !!names(.[1])) %>% 
  rbind(RD) %>% 
  ggplot() + aes(Density, fill = Type) +
    geom_density(alpha = 0.2) + 
    scale_fill_manual(values = c("yellow", "darkgreen")) + 
    labs(x = "Weighted Density", y = "Probability Density") + 
    scale_x_continuous(limits = c(0, 280), expand = c(0, 0), breaks = seq(0, 280, 40)) + 
    scale_y_continuous(limits = c(0, 0.02), expand = c(0, 0))

# Note that this graph does NOT show the recruitment probability at a given density
# It only shows the probability of a recruit occuring at a density GIVEN the plant is a recruit

# Clean up
remove(RD, TWD, TWD.BS)





##### Choosing a model for seed recruitment probability ---------------------------------------------------------------------------------------------

# We choose to use the partial linear model to calculate seed numbers
plot(CData.s$d.stand, CData.s$recruit.prob)

# Calculating the probability with FMod instead of PMOD highlights more artefacts from the linear model
# Thus, we use PMod to reduce the error of using this model and inform the recruitment probability model

# Recruitment will depend on local density only, for the most part
# Here, we try both a linear model and a quadratic model
# Quadratic model has much better p-value, but still low R^2 value
summary(lm(recruit.prob ~ d.stand, data = CData.s))
summary(lm(recruit.prob ~ d.stand + I(d.stand^2), data = CData.s))

# Quadratic model has much better p-values; both models have low R^2
# Next, check AIC values
AIC(lm(recruit.prob ~ d.stand, data = CData.s))
AIC(lm(recruit.prob ~ d.stand + I(d.stand^2), data = CData.s))

# Linear model has better AIC, so we will use that





##### Determining the transition between low-survival and high-survival regimes ---------------------------------------------------------------------

# After a certain volume "cutoff" point, almost no plants die
# This creates a quasi-perfect separation that linear models can't really handle
# Where is the volume cutoff?
subset(CData.AllSurvival.s) %>% 
  mutate(v.bin = cut_number(volume_t, n = 100)) %>% 
  group_by(v.bin) %>% 
  summarise(s.mean = mean(survival_t1)) %>% 
  print(n = Inf) %>% 
  ggplot() +
  geom_point(aes(x = v.bin, y = s.mean))

# Seems that volume cutoff should be at volume_t = 7.294?
# Check this against distribution of transplants and non-transplants
# Why? Transplants almost always die, while non-transplants almost always live
arrange(CData.AllSurvival.s, desc(transplant)) %>% 
  mutate(transplant = factor(transplant, levels = c(TRUE, FALSE))) %>% 
  ggplot() + aes(x = volume_t, fill = transplant) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = 7.294) +
  scale_fill_manual(values = c("TRUE" = "dodgerblue1", "FALSE" = "firebrick2"),
                    breaks = c("FALSE", "TRUE"))

# Seems to be a good cutoff, but there are still transplants above it and non-transplants below it
# Check percentages of each category on "correct" side of cutoff
CData.AllSurvival.s %>% 
  filter(volume_t < 7.294 & transplant == TRUE) %>% 
  nrow()/nrow(subset(CData.AllSurvival.s, transplant == TRUE))
CData.AllSurvival.s %>% 
  filter(volume_t >= 7.294 & transplant == FALSE) %>% 
  nrow()/nrow(subset(CData.AllSurvival.s, transplant == FALSE))

# Thus, we use linear models for any volume below cutoff, and guaranteed survival for anything above





##### Graphs for fitting lognormal distribution to terminal velocities ------------------------------------------------------------------------------

# Set up function to create Q-Q plot and PDF
tv.ln.graphs <- function(){
  
  # Create Q-Q plot for lognormal distribution
  t <- seq(0, 1, length.out = length(tv.raw) + 2)[-c(1, length(tv.raw))]
  qt <- qlnorm(t, meanlog = tv.fits[1], sdlog = tv.fits[2])
  qs <- quantile(tv.raw, t)
  
  # Set graph layout
  layout(mat = matrix(c(1, 2,
                        1, 2), 2, byrow = TRUE))
  
  # Graph Q-Q plot on left, with red line for theoretical lognormal distribution
  plot(qt, qs, main = "Lognormal Q-Q Plot", xlab = "Theoretical Quantiles",
       ylab = "Sample Quantiles", pch = 20); abline(0, 1, col = 'red')
  
  # Create dummy variable for plotting
  x <- seq(0, 6, length.out = 500)
  
  # Graph wind speed distribution on right
  # Red line is density estimated directly from data, black is from fit lognormal distribution
  plot(x, dlnorm(x, meanlog = tv.fits[1], sdlog = tv.fits[2]),
       xlab = "Terminal Velocity (m/s)", ylab = "Probability Density", pch = 20, type = "l",
       main = "Terminal Velocity Distribution (Lognormal)", ylim = c(0, 1), xlim = c(0, 6))
  lines(tv.PDF, col = 'red')
  
  # Reset graphical parameters
  par(mfrow = c(1, 1))
  
  # Check to make sure area under curve is 1
  tv.dist.lognormal <- function(u){
    ((u*tv.fits[2]*(sqrt(2*pi)))^(-1))*exp(-((log(u) - tv.fits[1])^2)/(2*tv.fits[2]^2))}
  print(integrate(tv.dist.lognormal, 0, Inf))}

# Call function to display graphs; then clean it up after done
tv.ln.graphs(); remove(tv.ln.graphs)





##### Graphs for fitting lognormal distribution to wind speeds --------------------------------------------------------------------------------------

# Set up function to create Q-Q plot and PDF
ws.ln.graphs <- function(){
  
  # Create Q-Q plot for lognormal distribution
  t <- seq(0, 1, length.out = length(ws.raw) + 2)[-c(1, length(ws.raw))]
  qt <- qlnorm(t, meanlog = ws.fits[1], sdlog = ws.fits[2])
  qs <- quantile(ws.raw, t)
  
  # Set graph layout
  layout(mat = matrix(c(1, 2,
                        1, 2), 2, byrow = TRUE))
  
  # Graph Q-Q plot on left, with red line for theoretical lognormal distribution
  plot(qt, qs, main = "Lognormal Q-Q Plot", xlab = "Theoretical Quantiles",
       ylab = "Sample Quantiles", pch = 20); abline(0, 1, col = 'red')
  
  # Create dummy variable for plotting
  x <- seq(0, 20, length.out = 2000)
  
  # Graph wind speed distribution on right
  # Red line is density estimated directly from data, black is from fit lognormal distribution
  plot(x, dlnorm(x, meanlog = ws.fits[1], sdlog = ws.fits[2]),
       xlab = "Wind Speed (m/s)", ylab = "Probability Density", pch = 20, type = "l",
       main = "Wind Speed Distribution (Lognormal)", ylim = c(0, 0.4), xlim = c(0, 20))
  lines(ws.PDF, col = 'red')
  
  # Reset graphical parameters
  par(mfrow = c(1, 1))
  
  # Check to make sure area under curve is 1
  ws.dist.lognormal <- function(u){
    ((u*ws.fits[2]*(sqrt(2*pi)))^(-1))*exp(-((log(u) - ws.fits[1])^2)/(2*ws.fits[2]^2))}
  print(integrate(ws.dist.lognormal, 0, Inf))}

# Call function to display graphs; then clean it up after done
ws.ln.graphs(); remove(ws.ln.graphs)





##### Graphs for fitting gamma distribution to terminal velocities ----------------------------------------------------------------------------------

# Set up function to create Q-Q plot and PDF
tv.gm.graphs <- function(){
  
  # Create Q-Q plot for gamma distribution
  t <- seq(0, 1, length.out = length(tv.raw) + 2)[-c(1, length(tv.raw))]
  qt <- qgamma(t, shape = tv.fits[3], rate = tv.fits[4])
  qs <- qgamma(pgamma(sort(tv.raw), shape = tv.fits[3], rate = tv.fits[4]),
               shape = tv.fits[3], rate = tv.fits[4])
  
  # Set graph layout
  layout(mat = matrix(c(1, 2,
                        1, 2), 2, byrow = TRUE))
  
  # Graph Q-Q plot on left, with red line for theoretical gamma distribution
  plot(qt, qs, main = "Gamma Q-Q Plot", xlab = "Theoretical Quantiles",
       ylab = "Sample Quantiles", pch = 20); abline(0, 1, col = 'red')
  
  # Create dummy variable for plotting
  x <- seq(0, 6, length.out = 500)
  
  # Graph wind speed distribution on right
  # Red line is density estimated directly from data, black is from fit gamma distribution
  plot(x, dgamma(x, shape = tv.fits[3], rate = tv.fits[4]),
       xlab = "Terminal Velocity (m/s)", ylab = "Probability Density", pch = 20, type = "l",
       main = "Terminal Velocity Distribution (Gamma)", ylim = c(0, 1), xlim = c(0, 6))
  lines(tv.PDF, col = 'red')
  
  # Reset graphical parameters
  par(mfrow = c(1, 1))
  
  # Check to make sure area under curve is 1
  tv.dist.gamma <- function(u){
    (tv.fits[4]^tv.fits[3])*(u^(tv.fits[3]-1))*exp(-tv.fits[4]*u)/gamma(tv.fits[3])}
  print(integrate(tv.dist.gamma, 0, Inf))}

# Call function to display graphs; then clean it up after done
tv.gm.graphs(); remove(tv.gm.graphs)





##### Graphs for fitting gamma distribution to wind speeds ------------------------------------------------------------------------------------------

# Set up function to create Q-Q plot and PDF
ws.gm.graphs <- function(){
  
  # Create Q-Q plot for gamma distribution
  t <- seq(0, 1, length.out = length(ws.raw) + 2)[-c(1, length(ws.raw))]
  qt <- qgamma(t, shape = ws.fits[3], rate = ws.fits[4])
  qs <- qgamma(pgamma(sort(ws.raw), shape = ws.fits[3], rate = ws.fits[4]),
               shape = ws.fits[3], rate = ws.fits[4])
  
  # Set graph layout
  layout(mat = matrix(c(1, 2,
                        1, 2), 2, byrow = TRUE))
  
  # Graph Q-Q plot on left, with red line for theoretical gamma distribution
  plot(qt, qs, main = "Gamma Q-Q Plot", xlab = "Theoretical Quantiles",
       ylab = "Sample Quantiles", pch = 20); abline(0, 1, col = 'red')
  
  # Create dummy variable for plotting
  x <- seq(0, 20, length.out = 2000)
  
  # Graph wind speed distribution on right
  # Red line is density estimated directly from data, black is from fit gamma distribution
  plot(x, dgamma(x, shape = ws.fits[3], rate = ws.fits[4]),
       xlab = "Wind Speed (m/s)", ylab = "Probability Density", pch = 20, type = "l",
       main = "Wind Speed Distribution (Gamma)", ylim = c(0, 0.4), xlim = c(0, 20))
  lines(ws.PDF, col = 'red')
  
  # Reset graphical parameters
  par(mfrow = c(1, 1))
  
  # Check to make sure area under curve is 1
  ws.dist.gamma <- function(u){
    (ws.fits[4]^ws.fits[3])*(u^(ws.fits[3]-1))*exp(-ws.fits[4]*u)/gamma(ws.fits[3])}
  print(integrate(ws.dist.gamma, 0, Inf))}

# Call function to display graphs; then clean it up after done
ws.gm.graphs(); remove(ws.gm.graphs)





##### Numerically comparing dispersal between heights -----------------------------------------------------------------------------------------------

# This may take a few minutes, so please be patient
# Also, results from the same code may differ slightly since random sampling occurs

# How likely is it that propagules exceed 5 and 10 m for the tallest shrub?
probs1 <- c()
probs2 <- c()
for(i in 1:50){
  sapply(max(CData$max.ht_t1, na.rm = TRUE)/100, WALD.f.e.h.2) -> dist
  probs1 <- append(probs1, length(dist[dist > 5])/length(dist))
  probs2 <- append(probs2, length(dist[dist > 10])/length(dist))}
mean(probs1); mean(probs2)

# What proportion of seeds fall within 1 m of the shrub?
probs1 <- c()
for(i in 1:50){
  sapply(max(CData$max.ht_t1, na.rm = TRUE)/100, WALD.f.e.h.2) -> dist
  probs1 <- append(probs1, length(dist[dist < 1])/length(dist))}
mean(probs1)

# How likely is it that propagules exceed 5 and 10 m for a 1 m shrub?
probs1 <- c()
probs2 <- c()
for(i in 1:50){
  sapply(1, WALD.f.e.h.2) -> dist
  probs1 <- append(probs1, length(dist[dist > 5])/length(dist))
  probs2 <- append(probs2, length(dist[dist > 10])/length(dist))}
mean(probs1); mean(probs2)

# How likely is it that propagules exceed 5 and 10 m for the median shrub height?
probs1 <- c()
probs2 <- c()
for(i in 1:100){
  sapply(median(CData$max.ht_t1, na.rm = TRUE)/100, WALD.f.e.h.2) -> dist
  probs1 <- append(probs1, length(dist[dist > 5])/length(dist))
  probs2 <- append(probs2, length(dist[dist > 10])/length(dist))}
mean(probs1); mean(probs2)

# Clean up mess from for loop
remove(i, dist, probs1, probs2)


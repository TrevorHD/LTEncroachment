##### Plot vital rates ------------------------------------------------------------------------------------------------------------------------------

# Code author(s): Tom

# Set bin splits for size and density
size_breaks <- 4
density_breaks <- 5
LATR_cols <- wes_palette("Zissou1", size_breaks, type = "continuous")

# Output graphics as PDF
pdf("Manuscript/Figures/VitalRates.pdf", useDingbats = F, height = 9, width = 7)
par(mfrow = c(3, 2), mar = c(5, 5, 1, 1))

# Plot survival (mature plants)
LATR_surv_dat %>% 
  filter(transplant == F) %>% 
  mutate(size_bin = as.numeric(cut(log_volume_t, breaks = size_breaks)),
         density_bin = as.numeric(cut(weighted.dens, breaks = density_breaks))) %>% 
  group_by(size_bin, density_bin) %>% 
  summarise(mean_size = mean(log_volume_t),
            mean_dens = mean(weighted.dens),
            mean_surv = mean(survival_t1),
            bin_n = n()) -> LATR_surv_plot_transF
LATR_surv_plot_transF$size_col <- LATR_cols[LATR_surv_plot_transF$size_bin]
mean_size <- LATR_surv_plot_transF %>% group_by(size_bin) %>% summarise(mean_size = mean(mean_size))
surv_predict <- expand.grid(
  weighted.dens = seq(min(LATR_surv_dat$weighted.dens), max(LATR_surv_dat$weighted.dens), 1),
  size_bin = mean_size$size_bin,
  transplant = F,
  unique.transect = LATR_flow_dat$unique.transect[1])
surv_predict$log_volume_t <- mean_size$mean_size[surv_predict$size_bin]
surv_predict$pred <- predict.gam(LATR_surv_best, newdata = surv_predict, type = "response",
                                 exclude = "s(unique.transect)")
surv_predict$size_col <- LATR_cols[surv_predict$size_bin]
plot(LATR_surv_plot_transF$mean_dens, LATR_surv_plot_transF$mean_surv,
     col = alpha(LATR_surv_plot_transF$size_col, 0.5), pch = 16,
     ylim = c(0, 1), xlab = "Weighted density", ylab = "Survival probability", cex.lab = 1.4,
     cex = LATR_surv_plot_transF$bin_n/max(LATR_surv_plot_transF$bin_n)*2 + 1,
     xlim = c(min(LATR_surv_dat$weighted.dens), max(LATR_surv_dat$weighted.dens)))
points(surv_predict$weighted.dens, surv_predict$pred, col = surv_predict$size_col, pch = 16, cex = 0.75)
title("A", adj = 0)

# Plot survival (transplants)
LATR_surv_dat %>% 
  filter(transplant == T) %>% 
  mutate(density_bin = as.numeric(cut(weighted.dens, breaks = density_breaks))) %>% 
  group_by(density_bin) %>% 
  summarise(mean_dens = mean(weighted.dens),
            mean_surv = mean(survival_t1),
            bin_n = n()) -> LATR_surv_plot_transT
surv_trans_predict <- expand.grid(
  weighted.dens = seq(min(LATR_surv_dat$weighted.dens), max(LATR_surv_dat$weighted.dens), 1),
  log_volume_t = mean(LATR_surv_dat$log_volume_t[LATR_surv_dat$transplant == T]),
  transplant = T,
  unique.transect = LATR_flow_dat$unique.transect[1])
surv_trans_predict$pred <- predict.gam(LATR_surv_best, newdata = surv_trans_predict, type = "response",
                                       exclude = "s(unique.transect)")
points(LATR_surv_plot_transT$mean_dens, LATR_surv_plot_transT$mean_surv, col = alpha("black", 0.5), pch = 16,
       cex.lab = 1.4, cex = LATR_surv_plot_transT$bin_n/max(LATR_surv_plot_transT$bin_n)*2 + 1)
lines(surv_trans_predict$weighted.dens, surv_trans_predict$pred, lwd = 3)
# points(LATR_surv_dat$weighted.dens[LATR_surv_dat$transplant == T],
#        LATR_surv_dat$survival_t1[LATR_surv_dat$transplant == T],
#        col = alph)
text(30, 0.1, "Transplants", font = 3)

# First pass at plotting growth, starting with raw data visualization, binning sizes and densities
# Odd at first since there is clear positive density dependence at smallest size, regardless of binning
LATR_grow %>% 
  mutate(size_bin = as.numeric(cut(log_volume_t, breaks = size_breaks)),
         density_bin = as.numeric(cut(weighted.dens, breaks = density_breaks))) %>% 
  group_by(size_bin, density_bin) %>% 
  summarise(mean_size = mean(log_volume_t),
            mean_dens = mean(weighted.dens),
            mean_sizet1 = mean(log_volume_t1),
            sd_sizet1 = sd(log_volume_t1),
            bin_n = n()) -> LATR_grow_plot
LATR_grow_plot$size_col <- LATR_cols[LATR_grow_plot$size_bin]
# plot(LATR_grow_plot$mean_dens, LATR_grow_plot$mean_sizet1, col = LATR_grow_plot$size_bin,
#      pch = 16, cex = 2)
# plot(LATR_grow_plot$mean_dens, (LATR_grow_plot$mean_sizet1 - LATR_grow_plot$mean_size),
#      col = LATR_grow_plot$size_bin, pch = 16, cex = 2)

# Bins actually differ in initial size, which may explain the apparent positive DD
# See commented code above for plotting the change in size to control for size differences

# Thus, we avoid binning the growth data for this read and create dummy data frame for GAM prediction
mean_size <- LATR_grow %>% 
  mutate(size_bin = as.numeric(cut(log_volume_t, breaks = size_breaks))) %>% 
  group_by(size_bin) %>% 
  summarise(mean_size=mean(log_volume_t))
grow_predict <- expand.grid(
  weighted.dens = seq(min(LATR_grow$weighted.dens), max(LATR_grow$weighted.dens), 1),
  size_bin = mean_size$size_bin,
  unique.transect = LATR_grow$unique.transect[1])
grow_predict$log_volume_t <- mean_size$mean_size[grow_predict$size_bin]
grow_predict[, c("pred_mu", "pred_sigma")] <- predict.gam(LATR_grow_best, newdata = grow_predict,
                                                          type = "response", exclude = "s(unique.transect)")
LATR_grow$size_col <- LATR_cols[as.numeric(cut(LATR_grow$log_volume_t, breaks = size_breaks))]
grow_predict$size_col <- LATR_cols[grow_predict$size_bin]

# Plot mean growth with SD inset
plot(LATR_grow$weighted.dens, LATR_grow$log_volume_t1, pch = 16, cex = 0.5, cex.lab = 1.4,
     col = alpha(LATR_grow$size_col, 0.5), xlab = "Weighted density",
     ylab = expression(paste(Size[t + 1], " (log ", cm^3, ")")))
points(grow_predict$weighted.dens, grow_predict$pred_mu, col = grow_predict$size_col, pch = 16, cex = 0.75)
rect(100, -2.5, 220, 6, col = "white", lwd = 0.5)
plotInset(100, -1, 220, 6,
          expr = plot(grow_predict$weighted.dens, 1/(grow_predict$pred_sigma^2),
                      col = grow_predict$size_col, pch = 16, cex = 0.25, xlab = "",
                      ylab = expression(paste("SD(", Size[t + 1], ")")),
                      cex.axis = 0.5, mgp = c(3/2, 1/2, 0)),
          mar = c(0, 3, 0, 0))
title("B", adj = 0)

# Plot flowering
LATR_flow_dat %>% 
  mutate(size_bin = as.numeric(cut(log_volume_t, breaks = size_breaks)),
         density_bin = as.numeric(cut(weighted.dens, breaks = density_breaks))) %>% 
  group_by(size_bin, density_bin) %>% 
  summarise(mean_size = mean(log_volume_t),
            mean_dens = mean(weighted.dens),
            mean_flow = mean(total.reproduction_t > 0),
            bin_n = n()) -> LATR_flow_plot
LATR_flow_plot$size_col <- LATR_cols[LATR_flow_plot$size_bin]
mean_size <- LATR_flow_plot %>% group_by(size_bin) %>% summarise(mean_size = mean(mean_size))
flow_predict <- expand.grid(
  weighted.dens = seq(min(LATR_flow_dat$weighted.dens), max(LATR_flow_dat$weighted.dens), 1),
  size_bin = mean_size$size_bin,
  unique.transect = LATR_flow_dat$unique.transect[1])
flow_predict$log_volume_t <- mean_size$mean_size[flow_predict$size_bin]
flow_predict$pred <- predict.gam(LATR_flower_best, newdata = flow_predict, type = "response",
                                 exclude = "s(unique.transect)")
flow_predict$size_col <- LATR_cols[flow_predict$size_bin]
plot(LATR_flow_plot$mean_dens, LATR_flow_plot$mean_flow, col = alpha(LATR_flow_plot$size_col, 0.5), pch = 16,
     ylim = c(0, 1), xlab = "Weighted density", ylab = "Flowering probability",
     cex = LATR_flow_plot$bin_n/max(LATR_flow_plot$bin_n)*2 + 1, cex.lab = 1.4,
     xlim = c(min(LATR_flow_dat$weighted.dens), max(LATR_flow_dat$weighted.dens)))
points(flow_predict$weighted.dens, flow_predict$pred, col = flow_predict$size_col, pch = 16, cex = 0.75)
title("C", adj = 0)

# Plot fruits
# Note that these size breaks are *different* compared to other vital rates (b/c this is the flowering subset)
mean_size <- LATR_fruits_dat %>% 
  mutate(size_bin = as.numeric(cut(log_volume_t, breaks = size_breaks))) %>% 
  group_by(size_bin) %>% 
  summarise(mean_size = mean(log_volume_t))
fruits_predict <- expand.grid(
  weighted.dens = seq(min(LATR_fruits_dat$weighted.dens), max(LATR_fruits_dat$weighted.dens), 1),
  size_bin = mean_size$size_bin,
  unique.transect = LATR_fruits_dat$unique.transect[1])
LATR_fruits_dat$size_col <- LATR_cols[as.numeric(cut(LATR_fruits_dat$log_volume_t, breaks = size_breaks))]
fruits_predict$log_volume_t <- mean_size$mean_size[fruits_predict$size_bin]
fruits_predict$pred <- predict.gam(LATR_fruits_best, newdata = fruits_predict, type = "response",
                                   exclude = "s(unique.transect)")
fruits_predict$size_col <- LATR_cols[fruits_predict$size_bin]
plot(LATR_fruits_dat$weighted.dens, LATR_fruits_dat$total.reproduction_t, pch = 16, cex = 0.5,
     ylim = c(0, 4000), cex.lab = 1.4, col = alpha(LATR_fruits_dat$size_col, 0.5),
     xlab = "Weighted density", ylab = "Flowers and fruits")
points(fruits_predict$weighted.dens, fruits_predict$pred, col = fruits_predict$size_col, pch = 16, cex = 0.75)
title("D", adj = 0)

# Plot recruitment
recruit_predict <- expand.grid(
  weighted.dens = seq(min(LATR_recruitment$weighted.dens), max(LATR_recruitment$weighted.dens), 1),
  unique.transect = LATR_fruits_dat$unique.transect[1])
recruit_predict$pred <- predict.gam(LATR_recruit_best, newdata = recruit_predict, type = "response",
                                    exclude = "s(unique.transect)")
plot(LATR_recruitment$weighted.dens, LATR_recruitment$recruits/LATR_recruitment$total_seeds, cex.lab = 1.4,
     col = alpha("black", 0.5), pch = 1, xlab = "Weighted density", ylab = "Per-seed recruitment probability",
     ylim = c(0, 0.001))
lines(recruit_predict$weighted.dens, recruit_predict$pred, lwd = 3)
title("E", adj = 0)
plot(LATR_recruit_size$weighted.dens, LATR_recruit_size$log_volume, col = alpha("black", 0.5), cex.lab = 1.4,
     pch = 16, xlab = "Weighted density", ylab = expression(paste("Recruit size (log ", cm^3, ")")))
lines(LATR_recruit_size$weighted.dens, LATR_recruit_size$pred[, 1], lwd = 3)
title("F", adj = 0)

# Turn off graphics device
dev.off()





##### Plot growth rate as function of density -------------------------------------------------------------------------------------------------------

# Get bootstrapped results
boot.lambda <- read.csv("https://github.com/TrevorHD/LTEncroachment/raw/master/Data/Derived/Boot_Lambda.csv")

# Compute mean result with full data
boot.prop <- 1
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/06_BootRes.R")
source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/05_CDataAnalysis.R")
mean.lambda <- LambdaD()

# Output graphics as PDF
pdf("Manuscript/Figures/LambdaDensity.pdf", useDingbats = F, height = 5, width = 6)
par(mar = c(5, 5, 1, 1))

# Plot bootstrapped lambda
plot(boot.lambda$density, boot.lambda[, 3], type = "n", xlab = "Weighted density",
     ylab = expression(paste(lambda)), cex.lab = 1.2, ylim = c(1, 1.06))
for(i in 3:dim(boot.lambda)[2]){
  lines(boot.lambda$density, boot.lambda[, i], col = alpha("darkgray", 0.15))}
lines(boot.lambda$density, mean.lambda, lwd = 4)

# Turn off graphics device
dev.off()

# Calculate fraction of samples with peak fitness at zero density
max.fit <- c()
for(i in 3:dim(boot.lambda)[2]){
  max.fit[i - 2] <- which.max(boot.lambda[, i])}
sum(max.fit == 1)/length(max.fit)

# Find "weird" samples
# Here are 3 weird ones from the first 100 samples
which(boot.lambda[25, ] > 1.005)
plot(boot.lambda$density, boot.lambda[, 3], type = "n", xlab = "Weighted density",
     ylab = expression(paste(lambda)), cex.lab = 1.2, ylim = c(1, 1.06))
for(i in 3:dim(boot.lambda)[2]){
  lines(boot.lambda$density, boot.lambda[, i], col = alpha("black", 0.15))}
lines(boot.lambda$density, boot.lambda[, 30], col = "red", lwd = 4)
lines(boot.lambda$density, boot.lambda[, 52], col = "blue", lwd = 4)
lines(boot.lambda$density, boot.lambda[, 100], col = "darkgreen", lwd = 4)

# Find reps with max lambda at density > 0
allee <- which(apply(as.matrix(boot.lambda[, 3:102]), 2, which.max) > 1)
lines(boot.lambda$density, boot.lambda[, 28 + 3], lwd = 4, lty = 2)





##### Plot dispersal kernels ------------------------------------------------------------------------------------------------------------------------

# Get hights and set colour palette
heights <- quantile(LATR_full$max.ht_t, na.rm = T)
D_cols <- wes_palette("Zissou1", 4, type = "continuous")

# Output graphics as PDF
pdf("Manuscript/Figures/DispersalKernels.pdf", useDingbats = F, height = 5, width = 6)
par(mar = c(5, 5, 1, 1))

# Plot kernels
plot(density(WALD.f.e.h.tom(10000, H = 0.6, elas = "none", sens = "none"), from = 0, to = 10,
             bw = 0.05), main = " ", type = "n", cex.lab = 1.4,
     lwd = 4, xlim = c(0, 4),
     ylim = c(0, 5), xlab = "Seed dispersal distance (m)", 
     ylab = "Probability Density")
for(i in 2:length(heights)){
  lines(density(WALD.f.e.h.tom(10000, H = heights[i]/100, elas = "none", sens = "none"),
                from = 0, to = 10, bw = 0.05), lwd = 4, col = D_cols[i - 1])}
legend("topright", title = "Shrub height", bty = "n",
       legend = c("0.45m (25%)", "0.62m (50%)", "0.83m (75%)", "1.96m (100%)"),
       lwd = 3, col = D_cols, cex = 1.2)

# Turn off graphics device
dev.off()





##### Plot wavespeeds and sensitivity ---------------------------------------------------------------------------------------------------------------

# Load in wavespeed and sensitivity data
boot.wave <- read.csv("https://github.com/TrevorHD/LTEncroachment/raw/master/Data/Derived/Boot_C2.csv")
boot.sens <- read.csv("https://github.com/TrevorHD/LTEncroachment/raw/master/Data/Derived/Boot_Sens.csv")
#boot.elas <- read.csv("https://github.com/TrevorHD/LTEncroachment/raw/master/Data/Derived/Boot_Elas.csv")

# Drop first row (garbage)
boot.sens <- boot.sens[-1, ]
names(boot.sens)[2:11] <- elas

# Output graphics as PDF
pdf("Manuscript/Figures/Wavespeeds.pdf", useDingbats = F, height = 5, width = 6)
par(mar = c(5, 5, 1, 1))

# Plot wavespeeds
plot(density(boot.wave[, 2]), type = "n", main = " ", cex.lab = 1.2,
     xlab = "Wavespeed (m/yr)", 
     ylab = "Probability Density")
polygon(density(boot.wave[, 2]), col = alpha(D_cols[4], .5))

# Turn off graphics device
dev.off()

# Plot alternative wavespeed/sensitivity combo
ggplot(boot.wave, aes(x = x)) +
  geom_density(color = "darkgrey",
               fill = "darkgrey",
               alpha = 0.5) +
  xlab("Wavespeed (m/yr)") +
  theme_classic() + 
  ggtitle('A') -> wave.plot
sens.dat <- data.frame(colMeans(boot.sens[, 2:11]))
names(sens.dat) <- "Wavespeed sensitivity"
sens.dat$vital_rate <- c("Growth mean", "Growth SD", "Survival",
                         "Flowering", "Fertility", "Recruitment",
                         "Recruit size mean", "Recruit size SD",
                         "Dispersal location", "Dispersal scale")
ggplot(data = sens.dat,
       aes(x = vital_rate, y = `Wavespeed sensitivity`))+
  geom_bar(stat = "identity") + 
  scale_y_break(c(0.05, 0.16)) +
  scale_y_break(c(0.182, 1.91)) +
  scale_y_break(c(1.925, 446.1)) +
  coord_flip() + 
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 45, size = 6),
        axis.title.y = element_blank()) + 
  ggtitle('B') -> sens.plot

# Output graphics as PDF
pdf("Manuscript/Figures/WavespeedSens.pdf", useDingbats = F, height = 10, width = 6)
wave.plot/sens.plot

# Turn off graphics device
dev.off()





##### Plot transect resurvey data -------------------------------------------------------------------------------------------------------------------

# Load data
resurv <- read.csv("https://github.com/TrevorHD/LTEncroachment/raw/master/Data/LT_TransectResurvey.csv")

# Output graphics as PDF
pdf("Manuscript/Figures/TransectResurveys.pdf", useDingbats = F, height = 4, width = 8)

# Plot resurvey data
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
plot(resurv$Quad[resurv$Year == 2001 & resurv$Transect == "East"],
     resurv$LATR_percent[resurv$Year == 2001 & resurv$Transect == "East"],
     type = "l", lwd = 2, col = "gray", ylim = c(0, 50),
     xlab = "Meter location", ylab = "Shrub cover (%)")
points(resurv$Quad[resurv$Year == 2013 & resurv$Transect == "East"],
       resurv$LATR_percent[resurv$Year == 2013 & resurv$Transect == "East"],
       type = "l", lwd = 2, col = "black")
title("A", adj = 0, font = 3)
plot(resurv$Quad[resurv$Year == 2001 & resurv$Transect == "West"],
     resurv$LATR_percent[resurv$Year == 2001 & resurv$Transect == "West"],
     type = "l", lwd = 2, col = "gray", ylim = c(0, 50),
     xlab = "Meter location", ylab = "Shrub cover (%)")
points(resurv$Quad[resurv$Year == 2013 & resurv$Transect == "West"],
       resurv$LATR_percent[resurv$Year == 2013 & resurv$Transect == "West"],
       type = "l", lwd = 2, col = "black")
title("B", adj = 0, font = 3)
legend("topright", legend = c("2001", "2013"), bty = "n", lwd = 2, col = c("gray", "black"))

# Turn off graphics device
dev.off()

# Find fraction of quads with zero creosote in 2001 and non-zero in 2013
resurv$trans_quad <- interaction(resurv$Transect, resurv$Quad)
zero2001 <- resurv[which(resurv$LATR_percent[resurv$Year == 2001] == 0), ]$trans_quad
nonzero2013 <- resurv[resurv$trans_quad %in% zero2001 & resurv$Year == 2013, ]$LATR_percent > 0
zero_to_nonzero <- sum(nonzero2013, na.rm = T)/length(na.omit(nonzero2013))

# Reverse of above: fraction of quads with nonzero in 2001 had zero in 2013?
nonzero2001 <- resurv[which(resurv$LATR_percent[resurv$Year == 2001] > 0), ]$trans_quad
zero2013 <- resurv[resurv$trans_quad %in% nonzero2001 & resurv$Year == 2013, ]$LATR_percent == 0
nonzero_to_zero <- sum(zero2013, na.rm = T)/length(na.omit(zero2013))

# Calculate overall mean cover in both years
mean_cover <- resurv %>% 
  group_by(Transect, as.factor(Year)) %>% 
  summarise(mean(LATR_percent, na.rm = T))
resurv %>% 
  group_by(as.factor(Year)) %>% 
  summarise(mean(LATR_percent, na.rm = T))
resurv %>% 
  group_by(as.factor(Year)) %>% 
  summarise(sum(LATR_percent > 0, na.rm = T))
resurv %>% 
  group_by(Transect,as.factor(Year)) %>% 
  summarise(sum(LATR_percent > 0, na.rm = T))
resurv %>% 
  group_by(Transect, as.factor(Year)) %>% 
  summarise(n())


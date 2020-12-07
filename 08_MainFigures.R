# ----- Plot bootstrapped wavespeeds ----------------------------------------------------------------------

# If you need wave speed data, uncomment the placeholder below
# This data came from a previous bootstrapped run
# boot.cv1 <- c(0.09300555, 0.07136389, 0.07973939, 0.06934206, 0.08390636, 0.10090147, 0.07331644,
#               0.07617637, 0.09215193, 0.07756776, 0.08905111, 0.09805242, 0.08950228 ,0.06447266,
#               0.06261184, 0.09056981, 0.07074311, 0.07019130, 0.07587839, 0.09814358, 0.07025593,
#               0.10539437, 0.11545569, 0.09206888, 0.03625982, 0.08895308, 0.07050688, 0.09903403, 
#               0.06102468, 0.11949837, 0.07716352, 0.08338977, 0.08842232, 0.10004509, 0.08942910,
#               0.08462161, 0.08994115, 0.09529132, 0.09029985, 0.13229806, 0.05929692, 0.04649399,
#               0.06152422, 0.07297820, 0.11533840, 0.10979438, 0.07682500, 0.09560812, 0.07322340,
#               0.08183116, 0.12628294, 0.10143668, 0.05289771, 0.10925509, 0.07852553, 0.13390366,
#               0.08429845, 0.08730562, 0.05310306, 0.07233902, 0.05901759, 0.07969898, 0.11591684,
#               0.08551957)

# Prepare graphics device
jpeg(filename = "Figure 1.jpeg", width = 750, height = 500, units = "px")

# Create blank page
grid.newpage()
plot.new()

# Set grid layout and activate it
gly <- grid.layout(500, 750)
pushViewport(viewport(layout = gly))

# Setup to plot data with equal point scales
pushViewport(viewport(layout = gly, layout.pos.row = 1:500, layout.pos.col = 1:750))
par(fig = gridFIG())
par(new = TRUE)
par(mar = c(5, 8, 4, 2))

# Plot PDF of bootstrapped wavespeeds
plot(density(boot.cv1), lwd = 3, xlim = c(0, 0.2), ylim = c(0, 25), 
     axes = FALSE, ann = FALSE)
axis(1, at = seq(0, 0.2, length.out = 5), cex.axis = 1.5, mgp = c(1, 1, 0))
axis(2, at = seq(0, 25, length.out = 6), cex.axis = 1.5, mgp = c(1, 1, 0), las = 1)
mtext("Wave Speed (m/yr)", side = 1, cex = 2, line = 3.5)
mtext("Probability Density", side = 2, cex = 2, line = 4)
box()
popViewport()

# Deactivate grid layout; finalise graphics save
popViewport()
dev.off()





# ----- Plot growth rate as function of density -----------------------------------------------------------

# Calculate growth rate lambda across a range of densities
d.test <- seq(min(LATR_full$weighted.dens, na.rm = TRUE), 
              max(LATR_full$weighted.dens, na.rm = TRUE), length.out = 25)
lambda_density <- c()
for(d in 1:length(d.test)){lambda_density[d] <- lambda(TransMatrix(dens = d.test[d], mat.size = 200)$IPMmat)}

# Prepare graphics device
jpeg(filename = "Figure 2.jpeg", width = 750, height = 500, units = "px")

# Create blank page
grid.newpage()
plot.new()

# Set grid layout and activate it
gly <- grid.layout(500, 750)
pushViewport(viewport(layout = gly))

# Setup to plot data with equal point scales
pushViewport(viewport(layout = gly, layout.pos.row = 1:500, layout.pos.col = 1:750))
par(fig = gridFIG())
par(new = TRUE)
par(mar = c(5, 8, 4, 2))

# Plot growth rate lambda across a range of densities
plot(d.test, lambda_density, type = "l", lwd = 3, ylim = c(1, 1.04), xlim = c(0, 200), 
     axes = FALSE, ann = FALSE)
axis(1, at = seq(0, 200, length.out = 5), cex.axis = 1.5, mgp = c(1, 1, 0))
axis(2, at = seq(1, 1.04, length.out = 5), cex.axis = 1.5, mgp = c(1, 1, 0), las = 1)
mtext("Weighted Density", side = 1, cex = 2, line = 3.5)
mtext("Growth Rate", side = 2, cex = 2, line = 5)
box()
popViewport()

# Deactivate grid layout; finalise graphics save
popViewport()
dev.off()

# Clean up mess from global environment
remove(gly, d.test, lambda_density)





# ----- Plot dispersal kernels ----------------------------------------------------------------------------

# Prepare graphics device
jpeg(filename = "Figure 3.jpeg", width = 1600, height = 900, units = "px")

# Create blank page
grid.newpage()
plot.new()

# Set grid layout and activate it
gly <- grid.layout(900, 1600)
pushViewport(viewport(layout = gly))

# Setup to plot dispersal kernels (large plot with NO log axis)
par(fig = gridFIG())
par(new = TRUE)
par(mar = c(7, 7, 2.1, 2.1))

# Create plot of dispersal kernels, with various drop heights (0.6 m, 0.8 m, 1 m, 1.2 m, 1.4 m)
plot(density(na.omit(WALD.f.e.h(0.6)), from = 0, to = 10, bw = 0.05), lwd = 4, xlim = c(0, 2),
     ylim = c(0, 5), xlab = "Distance (m)", ylab = "Probability Density", axes = FALSE, ann = FALSE)
axis(1, at = seq(0, 2, length.out = 5), cex.axis = 2, mgp = c(5, 2, 0))
axis(2, at = seq(0, 5, length.out = 6), cex.axis = 2, mgp = c(5, 2, 0), las = 1)
box()
mtext("Distance (m)", side = 1, cex = 3, line = 5)
mtext("Probability Density", side = 2, cex = 3, line = 4)
for(i in 1:4){
  colours <- c("black", "blue", "forestgreen", "red", "orange")
  lines(density(na.omit(WALD.f.e.h(0.6 + 0.2*i)), from = 0, to = 10, bw = 0.05), 
        col = colours[i + 1], lwd = 4)}
text(x = c(0.22, 0.31, 0.36, 0.42, 0.49), y = c(4, 2.25, 1.75, 1.45, 1.2), col = colours,
     labels = c("H = 0.6 m", "H = 0.8 m", "H = 1.0 m", "H = 1.2 m", "H = 1.4 m"), cex = 2)

# Setup to plot dispersal kernels (inset plot with log axis)
pushViewport(viewport(layout = gly, layout.pos.row = 32:600, layout.pos.col = 500:1569))
par(fig = gridFIG())
par(new = TRUE)
par(mar = c(7, 14, 0, 0))

# Create plot of dispersal kernels, with various drop heights (0.6 m, 0.8 m, 1 m, 1.2 m, 1.4 m)
plot(density(na.omit(WALD.f.e.h(0.6)), from = 0, to = 10, bw = 0.05), lwd = 4, cex.axis = 2,
     log = "x", xlab = "Distance (m)", ylab = "Probability Density", ann = FALSE, las = 1)
for(i in 1:4){
  colours <- c("black", "blue", "forestgreen", "red", "orange")
  lines(density(na.omit(WALD.f.e.h(0.6 + 0.2*i)), from = 0, to = 10, bw = 0.05), 
        col = colours[i + 1], lwd = 4)}
popViewport()

# Deactivate grid layout; finalise graphics save
popViewport()
dev.off()

# Clean up mess from for loop
remove(colours, i, gly)





##### Plot density effects on flowering and fruit production ----------------------------------------------

# Create vector of colours for plotting
PlotCol <- c("firebrick2", "darkorchid2", "dodgerblue1", "navy")

# Determine size of density and volume bins
n_cuts_dens <- 8
n_cuts_size <- 4

# Create bins for flowering data and model
LATR_flow_dat %>% 
  mutate(size_bin = as.integer(cut_number(log_volume_t, n_cuts_size)),
         dens_bin = as.integer(cut_number(weighted.dens, n_cuts_dens))) %>% 
  group_by(size_bin,dens_bin) %>% 
  summarise(mean_size = mean(log_volume_t),
            mean_density = mean(weighted.dens),
            mean_flower = mean(total.reproduction_t > 0),
            pred_flower = mean(pred),
            bin_n = n()) -> LATR_flow_dat_plot

# Generate flowering predictions for plotting
size_means_flow <- LATR_flow_dat_plot %>% group_by(size_bin) %>% summarise(mean_size=mean(mean_size))
LATR_flow_pred <- data.frame(
  weighted.dens = rep(seq(min(LATR_flow_dat$weighted.dens), max(LATR_flow_dat$weighted.dens), length.out = 20), times = n_cuts_size),
  log_volume_t = rep(size_means_flow$mean_size, each = 20),
  unique.transect = "1.FPS",
  size_bin = rep(size_means_flow$size_bin, each = 20))
LATR_flow_pred$pred <- predict.gam(LATR_flower_best, newdata = LATR_flow_pred, exclude = "s(unique.transect)")

# Create bins for fruits data and model
LATR_fruits_dat %>% 
  mutate(size_bin = as.integer(cut_number(log_volume_t, n_cuts_size)),
         dens_bin = as.integer(cut_number(weighted.dens, n_cuts_dens))) %>% 
  group_by(size_bin, dens_bin) %>% 
  summarise(mean_size = mean(log_volume_t),
            mean_density = mean(weighted.dens),
            mean_fruits = mean(total.reproduction_t),
            pred_fruits = mean(pred),
            bin_n = n()) -> LATR_fruits_dat_plot

# Generate fruits data for prediction
size_means_fruit <- LATR_fruits_dat_plot %>% group_by(size_bin) %>% summarise(mean_size = mean(mean_size))
LATR_fruit_pred <- data.frame(
  weighted.dens = rep(seq(min(LATR_fruits_dat$weighted.dens), max(LATR_fruits_dat$weighted.dens), length.out = 20), times = n_cuts_size),
  log_volume_t = rep(size_means_fruit$mean_size, each = 20),
  unique.transect = "1.FPS",
  size_bin = rep(size_means_fruit$size_bin, each = 20))
LATR_fruit_pred$pred <- predict.gam(LATR_fruits_best, newdata = LATR_fruit_pred, exclude = "s(unique.transect)")

# Prepare graphics device
jpeg(filename = "Figure 4.jpeg", width = 1000, height = 1200, units = "px")

# Create blank page
grid.newpage()
plot.new()

# Set grid layout and activate it
gly <- grid.layout(1200, 1000)
pushViewport(viewport(layout = gly))

# Setup to plot fruits data
pushViewport(viewport(layout = gly, layout.pos.row = 1:635, layout.pos.col = 25:1000))
par(fig = gridFIG())
par(new = TRUE)
par(mar = c(5, 6, 4, 2))

# Plot fruits data
plot(LATR_fruits_dat_plot$mean_density, LATR_fruits_dat_plot$mean_fruits, type = "n", 
     xlim = c(0, 200), ylim = c(0, 2500), axes = FALSE, ann = FALSE)
axis(1, at = seq(0, 200, length.out = 5), labels = FALSE, mgp = c(1, 1, 0))
axis(2, at = seq(0, 2500, length.out = 6), cex.axis = 2, mgp = c(1, 1, 0), las = 1)
mtext("Flowers and Fruits", side = 2, cex = 2.5, line = 5)
box()
for(i in 1:n_cuts_size){
  points(LATR_fruits_dat_plot$mean_density[LATR_fruits_dat_plot$size_bin == i],
         LATR_fruits_dat_plot$mean_fruits[LATR_fruits_dat_plot$size_bin == i], pch = 16, col = PlotCol[i],
         cex = (LATR_fruits_dat_plot$bin_n[LATR_fruits_dat_plot$size_bin == i]/max(LATR_fruits_dat_plot$bin_n))*4)
  lines(LATR_fruit_pred$weighted.dens[LATR_fruit_pred$size_bin == i],
        exp(LATR_fruit_pred$pred[LATR_fruit_pred$size_bin == i]), col = PlotCol[i], lwd = 3)}
popViewport()

# Setup to plot flowering data
pushViewport(viewport(layout = gly, layout.pos.row = 540:1175, layout.pos.col = 25:1000))
par(fig = gridFIG())
par(new = TRUE)
par(mar = c(5, 6, 4, 2))

# Plot flowering data
plot(LATR_flow_dat_plot$mean_density, LATR_flow_dat_plot$mean_flower, type = "n", 
     xlim = c(0, 200), ylim = c(0, 1), axes = FALSE, ann = FALSE)
axis(1, at = seq(0, 200, length.out = 5), cex.axis = 2, mgp = c(1, 1, 0))
axis(2, at = seq(0, 1, length.out = 6), cex.axis = 2, mgp = c(1, 1, 0), las = 1)
mtext("Weighted density", side = 1, cex = 2.5, line = 3.5)
mtext("Pr(Flowering)", side = 2, cex = 2.5, line = 5)
box()
for(i in 1:n_cuts_size){
  points(LATR_flow_dat_plot$mean_density[LATR_flow_dat_plot$size_bin == i],
         LATR_flow_dat_plot$mean_flower[LATR_flow_dat_plot$size_bin == i], pch = 16, col = PlotCol[i],
         cex = (LATR_flow_dat_plot$bin_n[LATR_flow_dat_plot$size_bin == i]/max(LATR_flow_dat_plot$bin_n))*4)
  lines(LATR_flow_pred$weighted.dens[LATR_flow_pred$size_bin == i],
        invlogit(LATR_flow_pred$pred[LATR_flow_pred$size_bin == i]), col = PlotCol[i], lwd = 3)}
popViewport()

# Create legend
grid.rect(vp = viewport(layout.pos.row = 100:125, layout.pos.col = 930:955), gp = gpar(fill = "navy"))
grid.rect(vp = viewport(layout.pos.row = 135:160, layout.pos.col = 930:955), gp = gpar(fill = "dodgerblue1"))
grid.rect(vp = viewport(layout.pos.row = 170:195, layout.pos.col = 930:955), gp = gpar(fill = "darkorchid2"))
grid.rect(vp = viewport(layout.pos.row = 205:230, layout.pos.col = 930:955), gp = gpar(fill = "firebrick2"))
grid.text(label = c("(75, 100]", "(50, 75]", "(25, 50]", "[0, 25]"),
          x = rep(0.923, 4), y = c(0.909, 0.880, 0.851, 0.822), just = "right", gp = gpar(fontsize = 20))
grid.text(label = bquote(underline("Size percentile")),
          x = 0.957, y = 0.936, just = "right", gp = gpar(fontsize = 26))

# Deactivate grid layout; finalise graphics save
popViewport()
dev.off()





##### Plot density effects on survival --------------------------------------------------------------------

# Visualize data and model for the natural census
n_cuts_dens <- 8
n_cuts_size <- 4
LATR_surv_dat %>% 
  filter(transplant == FALSE) %>% 
  mutate(size_bin = as.integer(cut_interval(log_volume_t, n_cuts_size)),
         dens_bin = as.integer(cut_interval(weighted.dens, n_cuts_dens))) %>% 
  group_by(size_bin, dens_bin) %>% 
  summarise(mean_size = mean(log_volume_t),
            mean_density = mean(weighted.dens),
            mean_surv = mean(survival_t1),
            pred_surv = mean(pred),
            bin_n = n()) -> LATR_surv_nat_plot

# Generate predictions for natural census plotting
size_means_surv_nat <- LATR_surv_nat_plot %>% group_by(size_bin) %>% summarise(mean_size = mean(mean_size))
LATR_surv_nat_pred <- data.frame(
  weighted.dens = rep(seq(0, max(LATR_surv_nat_plot$mean_density), length.out = 20), times = n_cuts_size),
  log_volume_t = rep(size_means_surv_nat$mean_size, each = 20),
  unique.transect = "1.FPS",
  transplant = FALSE,
  size_bin = rep(size_means_surv_nat$size_bin, each = 20))
LATR_surv_nat_pred$pred <- predict.gam(LATR_surv_best,newdata = LATR_surv_nat_pred, exclude = "s(unique.transect)")

# Generate predictions for transplant plotting
LATR_surv_dat %>% 
  filter(transplant == TRUE) %>% 
  mutate(dens_bin = as.integer(cut_interval(weighted.dens, n_cuts_dens)),
         mean_size = mean(log_volume_t)) %>% 
  group_by(dens_bin) %>% 
  summarise(mean_size = unique(mean_size),
            mean_density = mean(weighted.dens),
            mean_surv = mean(survival_t1),
            bin_n = n()) -> LATR_surv_exp_plot
LATR_surv_exp_pred <- data.frame(
  weighted.dens = seq(0, max(LATR_surv_exp_plot$mean_density), length.out = 20),
  log_volume_t = LATR_surv_exp_plot$mean_size[1],
  unique.transect = "1.FPS",
  transplant = TRUE)
LATR_surv_exp_pred$pred <- predict.gam(LATR_surv_best, newdata = LATR_surv_exp_pred, exclude = "s(unique.transect)")

# Prepare graphics device
jpeg(filename = "Figure 5.jpeg", width = 1000, height = 1200, units = "px")

# Create blank page
grid.newpage()
plot.new()

# Set grid layout and activate it
gly <- grid.layout(1200, 1000)
pushViewport(viewport(layout = gly))

# Setup to plot data with equal point scales
pushViewport(viewport(layout = gly, layout.pos.row = 1:635, layout.pos.col = 25:1000))
par(fig = gridFIG())
par(new = TRUE)
par(mar = c(5, 6, 4, 2))

# Plot natural census and transplant survival data with equal point scales
plot(LATR_surv_nat_plot$mean_density, LATR_surv_nat_plot$mean_surv, type = "n",
     xlim = c(0, 200), ylim = c(0, 1), axes = FALSE, ann = FALSE)
for(i in 1:n_cuts_size){
  points(LATR_surv_nat_plot$mean_density[LATR_surv_nat_plot$size_bin == i],
         LATR_surv_nat_plot$mean_surv[LATR_surv_nat_plot$size_bin == i], pch = 16, col = PlotCol[i], cex = 2.5)
  lines(LATR_surv_nat_pred$weighted.dens[LATR_surv_nat_pred$size_bin == i],
        invlogit(LATR_surv_nat_pred$pred[LATR_surv_nat_pred$size_bin == i]), col = PlotCol[i], lwd = 3)}
points(LATR_surv_exp_plot$mean_density, LATR_surv_exp_plot$mean_surv, ylim = c(0, 1), pch = 2, cex = 2.5)
lines(LATR_surv_exp_pred$weighted.dens, invlogit(LATR_surv_exp_pred$pred), lty = 2)
axis(1, at = seq(0, 200, length.out = 5), labels = FALSE, mgp = c(1, 1, 0))
axis(2, at = seq(0, 1, length.out = 6), cex.axis = 2, mgp = c(1, 1, 0), las = 1)
mtext("Pr(Survival)", side = 2, cex = 2.5, line = 5)
box()
popViewport()

# Setup to plot data with points scaled on number of observations
pushViewport(viewport(layout = gly, layout.pos.row = 540:1175, layout.pos.col = 25:1000))
par(fig = gridFIG())
par(new = TRUE)
par(mar = c(5, 6, 4, 2))

# Plot natural census and transplant survival data with points scaled on number of observations
plot(LATR_surv_nat_plot$mean_density, LATR_surv_nat_plot$mean_surv, type = "n",
     xlim = c(0, 200), ylim = c(0, 1), axes = FALSE, ann = FALSE)
for(i in 1:n_cuts_size){
  points(LATR_surv_nat_plot$mean_density[LATR_surv_nat_plot$size_bin == i],
         LATR_surv_nat_plot$mean_surv[LATR_surv_nat_plot$size_bin == i], pch = 16, col = PlotCol[i],
  cex = (LATR_surv_nat_plot$bin_n[LATR_surv_nat_plot$size_bin == i]/max(LATR_surv_nat_plot$bin_n))*7)
  lines(LATR_surv_nat_pred$weighted.dens[LATR_surv_nat_pred$size_bin == i],
        invlogit(LATR_surv_nat_pred$pred[LATR_surv_nat_pred$size_bin == i]), col = PlotCol[i], lwd = 3)}
points(LATR_surv_exp_plot$mean_density, LATR_surv_exp_plot$mean_surv, ylim = c(0, 1), pch = 2,
       cex = (LATR_surv_exp_plot$bin_n/max(LATR_surv_nat_plot$bin_n))*7)
lines(LATR_surv_exp_pred$weighted.dens, invlogit(LATR_surv_exp_pred$pred), lty = 2)
axis(1, at = seq(0, 200, length.out = 5), cex.axis = 2, mgp = c(1, 1, 0))
axis(2, at = seq(0, 1, length.out = 6), cex.axis = 2, mgp = c(1, 1, 0), las = 1)
mtext("Weighted density", side = 1, cex = 2.5, line = 3.5)
mtext("Pr(Survival)", side = 2, cex = 2.5, line = 5)
box()
popViewport()

# Block out lines to reduce clash with legend
grid.rect(vp = viewport(layout.pos.row = 65:90, layout.pos.col = 765:960), 
          gp = gpar(fill = "white", col = "white"))
grid.rect(vp = viewport(layout.pos.row = 130:190, layout.pos.col = 815:960), 
          gp = gpar(fill = "white", col = "white"))

# Create legend
grid.rect(vp = viewport(layout.pos.row = 100:125, layout.pos.col = 930:955), gp = gpar(fill = "navy"))
grid.rect(vp = viewport(layout.pos.row = 135:160, layout.pos.col = 930:955), gp = gpar(fill = "dodgerblue1"))
grid.rect(vp = viewport(layout.pos.row = 170:195, layout.pos.col = 930:955), gp = gpar(fill = "darkorchid2"))
grid.rect(vp = viewport(layout.pos.row = 205:230, layout.pos.col = 930:955), gp = gpar(fill = "firebrick2"))
grid.text(label = c("(10.7, 14.9]", "(6.45, 10.7]", "(2.21, 6.45]", "[-2.03, 2.21]"),
          x = rep(0.923, 4), y = c(0.909, 0.880, 0.851, 0.822), just = "right", gp = gpar(fontsize = 20))
grid.text(label = bquote(underline("Log-size interval")),
          x = 0.957, y = 0.936, just = "right", gp = gpar(fontsize = 26))

# Deactivate grid layout; finalise graphics save
popViewport()
dev.off()





# ----- Plot per-seed recruitment rate --------------------------------------------------------------------

# Add points at max and min density to predict recruitment rate
LATR_recruitment_line <- LATR_recruitment
LATR_recruitment_line[1324:1325, 1:8] <- LATR_recruitment_line[1322:1323, 1:8]
LATR_recruitment_line[1324, 4] <- 0; LATR_recruitment_line[1325, 4] <- 300
LATR_recruitment_line$pred <- predict.gam(LATR_recruit_best, newdata = LATR_recruitment_line, exclude = "s(unique.transect)")

# Prepare graphics device
jpeg(filename = "Figure 6.jpeg", width = 750, height = 500, units = "px")

# Create blank page
grid.newpage()
plot.new()

# Set grid layout and activate it
gly <- grid.layout(500, 750)
pushViewport(viewport(layout = gly))

# Setup to plot data with equal point scales
pushViewport(viewport(layout = gly, layout.pos.row = 1:500, layout.pos.col = 1:750))
par(fig = gridFIG())
par(new = TRUE)
par(mar = c(5, 8, 4, 2))

# Plot null model for per-seed recruitment
# No evidence for density dependence in recruitment, just a really low overall recruitment rate
plot(LATR_recruitment$weighted.dens, LATR_recruitment$recruits/LATR_recruitment$total_seeds, pch = 16,
     col = rgb(r = 0, g = 0, b = 0, alpha = 0.15), xlim = c(0, 300), ylim = c(0, 0.008), cex = 2,
     axes = FALSE, ann = FALSE)
axis(1, at = seq(0, 300, length.out = 4), cex.axis = 1.5, mgp = c(1, 1, 0))
axis(2, at = seq(0, 0.008, length.out = 5), cex.axis = 1.5, mgp = c(1, 1, 0), las = 1)
mtext("Weighted density", side = 1, cex = 2, line = 3.5)
mtext("Pr(Recruitment)", side = 2, cex = 2, line = 5)
lines(LATR_recruitment_line$weighted.dens, invlogit(LATR_recruitment_line$pred), col = "red")
box()
popViewport()

# Deactivate grid layout; finalise graphics save
popViewport()
dev.off()


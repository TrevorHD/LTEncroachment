# ----- Plot bootstrapped wavespeeds and growth rates -----------------------------------------------------

# If you need wave speed data, uncomment the placeholder below
#boot.cv1 <- c(0.185, 0.156, 0.158, 0.164, 0.175, 0.161, 0.17)

# Calculate growth rate lambda across a range of densities
d.test <- seq(min(LATR_full$weighted.dens, na.rm = TRUE), 
              max(LATR_full$weighted.dens, na.rm = TRUE), length.out = 25)
lambda_density <- c()
for(d in 1:length(d.test)){lambda_density[d] <- lambda(TransMatrix(dens = d.test[d], mat.size = 200)$IPMmat)}

# Prepare graphics device
jpeg(filename = "Figure 1.jpeg", width = 2600, height = 900, units = "px")

# Create blank page
grid.newpage()
plot.new()

# Set grid layout and activate it
gly <- grid.layout(900, 2600)
pushViewport(viewport(layout = gly))

# Setup to plot bootstrapped wavespeeds
pushViewport(viewport(layout = gly, layout.pos.row = 1:900, layout.pos.col = 1:1300))
par(fig = gridFIG())
par(new = TRUE)
par(mar = c(7, 9, 2.1, 2.1))

# Plot PDF of bootstrapped wavespeeds
plot(density(boot.cv1), col = "yellow3", lwd = 6, xlim = c(0.125, 0.2), ylim = c(0, 40), 
     axes = FALSE, ann = FALSE)
axis(1, at = seq(0.125, 0.2, length.out = 4), cex.axis = 2, mgp = c(5, 2, 0))
axis(2, at = seq(0, 40, length.out = 5), cex.axis = 2, mgp = c(5, 1, 0), las = 1)
mtext("Wave Speed (m/yr)", side = 1, cex = 3, line = 5)
mtext("Probability Density", side = 2, cex = 3, line = 5)
box()
#lines(density(boot.cv2), col = "green4", lwd = 6)
# Need to update code below to reflect new values once BS scenario is ready
#segments(x0 = c(-0.2, -0.2, s.c.min, s.c.min.2), x1 = c(s.c.min, s.c.min.2, s.c.min, s.c.min.2),
#         y0 = c(c.min, c.min.2, -0.05, -0.05), y1 = c(c.min, c.min.2, c.min, c.min.2), lty = 2)
#text(labels = c(paste0("~ ", round(c.min.2, digits = 4), " m/yr"),
#                paste0("~ ", round(c.min, digits = 4), " m/yr")),
#     x = rep(0.110, 2), y = c(c.min.2, c.min) + 0.0007, cex = 2.5)
popViewport()

# Setup to plot growth rates
pushViewport(viewport(layout.pos.row = 1:900, layout.pos.col = 1300:2600))
par(fig = gridFIG())
par(new = TRUE)
par(mar = c(7, 2.1, 2.1, 9))

# Plot growth rate lambda across a range of densities
plot(d.test, lambda_density, type = "l", lwd = 6, col = "yellow3", ylim = c(1, 1.04), xlim = c(0, 200), 
     axes = FALSE, ann = FALSE)
axis(1, at = seq(0, 200, length.out = 5), cex.axis = 2, mgp = c(5, 2, 0))
axis(4, at = seq(1, 1.04, length.out = 5), cex.axis = 2, mgp = c(5, 1, 0), las = 1)
mtext("Weighted Density", side = 1, cex = 3, line = 5)
mtext("Growth Rate", side = 4, cex = 3, line = 6.5)
box()
popViewport()

# Create legend
grid.rect(vp = viewport(layout.pos.row = 60:100, layout.pos.col = 2420:2460), gp = gpar(fill = "green4"))
grid.rect(vp = viewport(layout.pos.row = 120:160, layout.pos.col = 2420:2460), gp = gpar(fill = "yellow3"))
grid.text(label = c("Higher seedling survival", "Normal conditions"),
          x = rep(0.927, 2), y = c(0.912, 0.849), just = "right", gp = gpar(fontsize = 27))

# Deactivate grid layout; finalise graphics save
popViewport()
dev.off()

# Clean up mess from global environment
remove(gly, d.test, lambda_density)





# ----- Plot dispersal kernels ----------------------------------------------------------------------------

# Prepare graphics device
jpeg(filename = "Figure 2.jpeg", width = 1600, height = 900, units = "px")

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
jpeg(filename = "Figure 3.jpeg", width = 1000, height = 1200, units = "px")

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
     xlim = c(0, 200), ylim = c(0, 2500), cex.lab = 2, axes = FALSE, ann = FALSE)
axis(1, at = seq(0, 200, length.out = 5), labels = FALSE, mgp = c(1, 1, 0))
axis(2, at = seq(0, 2500, length.out = 6), cex.axis = 1.5, mgp = c(1, 1, 0), las = 1)
mtext("Flowers and Fruits", side = 2, cex = 2, line = 5)
box()
for(i in 1:n_cuts_size){
  points(LATR_fruits_dat_plot$mean_density[LATR_fruits_dat_plot$size_bin == i],
         LATR_fruits_dat_plot$mean_fruits[LATR_fruits_dat_plot$size_bin == i], pch = 16, col = PlotCol[i],
         cex = (LATR_fruits_dat_plot$bin_n[LATR_fruits_dat_plot$size_bin == i]/max(LATR_fruits_dat_plot$bin_n))*3)
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
     xlim = c(0, 200), ylim = c(0, 1), cex.lab = 2, axes = FALSE, ann = FALSE)
axis(1, at = seq(0, 200, length.out = 5), cex.axis = 1.5, mgp = c(1, 1, 0))
axis(2, at = seq(0, 1, length.out = 6), cex.axis = 1.5, mgp = c(1, 1, 0), las = 1)
mtext("Weighted density", side = 1, cex = 2, line = 3.5)
mtext("Pr(Flowering)", side = 2, cex = 2, line = 5)
box()
for(i in 1:n_cuts_size){
  points(LATR_flow_dat_plot$mean_density[LATR_flow_dat_plot$size_bin == i],
         LATR_flow_dat_plot$mean_flower[LATR_flow_dat_plot$size_bin == i], pch = 16, col = PlotCol[i],
         cex = (LATR_flow_dat_plot$bin_n[LATR_flow_dat_plot$size_bin == i]/max(LATR_flow_dat_plot$bin_n))*3)
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
jpeg(filename = "Figure 4.jpeg", width = 1000, height = 1200, units = "px")

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
     xlim = c(0, 200), ylim = c(0, 1), cex.lab = 2, axes = FALSE, ann = FALSE)
for(i in 1:n_cuts_size){
  points(LATR_surv_nat_plot$mean_density[LATR_surv_nat_plot$size_bin == i],
         LATR_surv_nat_plot$mean_surv[LATR_surv_nat_plot$size_bin == i], pch = 16, col = PlotCol[i], cex = 1.5)
  lines(LATR_surv_nat_pred$weighted.dens[LATR_surv_nat_pred$size_bin == i],
        invlogit(LATR_surv_nat_pred$pred[LATR_surv_nat_pred$size_bin == i]), col = PlotCol[i], lwd = 3)}
points(LATR_surv_exp_plot$mean_density, LATR_surv_exp_plot$mean_surv, ylim = c(0, 1), pch = 2)
lines(LATR_surv_exp_pred$weighted.dens, invlogit(LATR_surv_exp_pred$pred), lty = 2)
axis(1, at = seq(0, 200, length.out = 5), labels = FALSE, mgp = c(1, 1, 0))
axis(2, at = seq(0, 1, length.out = 6), cex.axis = 1.5, mgp = c(1, 1, 0), las = 1)
mtext("Pr(Survival)", side = 2, cex = 2, line = 5)
box()
popViewport()

# Setup to plot data with points scaled on number of observations
pushViewport(viewport(layout = gly, layout.pos.row = 540:1175, layout.pos.col = 25:1000))
par(fig = gridFIG())
par(new = TRUE)
par(mar = c(5, 6, 4, 2))

# Plot natural census and transplant survival data with points scaled on number of observations
plot(LATR_surv_nat_plot$mean_density, LATR_surv_nat_plot$mean_surv, type = "n",
     xlim = c(0, 200), ylim = c(0, 1), cex.lab = 2, axes = FALSE, ann = FALSE)
for(i in 1:n_cuts_size){
  points(LATR_surv_nat_plot$mean_density[LATR_surv_nat_plot$size_bin == i],
         LATR_surv_nat_plot$mean_surv[LATR_surv_nat_plot$size_bin == i], pch = 16, col = PlotCol[i],
  cex = (LATR_surv_nat_plot$bin_n[LATR_surv_nat_plot$size_bin == i]/max(LATR_surv_nat_plot$bin_n))*5)
  lines(LATR_surv_nat_pred$weighted.dens[LATR_surv_nat_pred$size_bin == i],
        invlogit(LATR_surv_nat_pred$pred[LATR_surv_nat_pred$size_bin == i]), col = PlotCol[i], lwd = 3)}
points(LATR_surv_exp_plot$mean_density, LATR_surv_exp_plot$mean_surv, ylim = c(0, 1), pch = 2,
       cex = (LATR_surv_exp_plot$bin_n/max(LATR_surv_nat_plot$bin_n))*5)
lines(LATR_surv_exp_pred$weighted.dens, invlogit(LATR_surv_exp_pred$pred), lty = 2)
axis(1, at = seq(0, 200, length.out = 5), cex.axis = 1.5, mgp = c(1, 1, 0))
axis(2, at = seq(0, 1, length.out = 6), cex.axis = 1.5, mgp = c(1, 1, 0), las = 1)
mtext("Weighted density", side = 1, cex = 2, line = 3.5)
mtext("Pr(Survival)", side = 2, cex = 2, line = 5)
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





# ----- ALL CODE BELOW WILL SOON BE DEPRECATED ------------------------------------------------------------

# ----- Flowering as a function of volume and density (1 of 4) --------------------------------------------

# Remove entries for which there is no recorded volume
# Then standardise density and introduce unique transect identifier
CData.s <- CData %>% 
  drop_na(volume_t) %>% 
  mutate(d.stand = (weighted.dens - mean(weighted.dens, na.rm = TRUE)) / sd(weighted.dens, na.rm = TRUE),
         unique.transect = interaction(transect, site))

# Calculate range of data
v.range <- range(CData.s$volume_t)
d.range <- range(CData.s$d.stand)

# Calculate predictions for 4 "bins" of volume categories
v.cut <- CData.s %>% 
  mutate(v.bin = cut_number(volume_t, n = 4)) %>% 
  group_by(v.bin) %>% 
  summarise(v.mean = mean(volume_t))
d.and.v <- as.tibble(expand.grid(v.cut$v.mean, seq(d.range[1], d.range[2], 0.5)))
names(d.and.v) <- c("v.mean", "x.d") 
F.pred <- d.and.v %>% 
  mutate(F.mean = invlogit(Mod.F.top.cf[1] + 
                           Mod.F.top.cf[2] * v.mean + 
                           Mod.F.top.cf[3] * x.d +
                           Mod.F.top.cf[4] * v.mean * x.d +
                           Mod.F.top.cf[5] * (x.d^2)+
                           Mod.F.top.cf[6] * v.mean * (x.d^2)))

# Get binned means of the data with respect to volume and density
F.mean.df <- CData.s %>% 
  mutate(v.bin = cut_number(volume_t, n = 4),
         d.bin = cut_number(d.stand, n = 6)) %>% 
  group_by(d.bin, v.bin) %>% 
  summarise(F.mean = mean(did.flower),
            d.mean = mean(d.stand),
            v.mean = mean(volume_t))

# Generate visualisation; flowering probability as function of density, with 4 volume bins
F.pred %>% 
  ggplot() +
  geom_line(aes(x = x.d, y = F.mean, colour = as.factor(v.mean)), size = 3) +
  scale_colour_manual(values = c("darkorchid2", "dodgerblue1", "navy", "firebrick2",
                                 "darkorchid2", "dodgerblue1", "navy", "firebrick2")) +
  geom_point(data = F.mean.df, aes(x = d.mean, y = F.mean, 
                                   colour = as.factor(v.bin)), size = 7) +
  labs(x = NULL, y = "Flowering Probability") +
  ylim(0, 1) +
  xlim(-1.5, 3) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        plot.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 38, margin = margin(t = 4, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 38, margin = margin(t = 0, r = 4, b = 0, l = 0)),
        axis.title.x = element_text(size = 43, margin = margin(t = 16, r = 0, b = 0, l = 6)),
        axis.title.y = element_text(size = 43, margin = margin(t = 80, r = 16, b = 0, l = 0)),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(8, "points")) -> Mod.F.top.dPlot

# Do same process as previous 3 steps, but switch density and volume
# This produces a graph of flowering probability as a function of volume, with 4 density bins
d.cut <- CData.s %>% 
  mutate(d.bin = cut_number(d.stand, n = 4)) %>% 
  group_by(d.bin) %>% 
  summarise(d.mean = mean(d.stand))
d.and.v.2 <- as.tibble(expand.grid(d.cut$d.mean, seq(v.range[1], v.range[2], 0.5)))
names(d.and.v.2) <- c("d.mean", "x.v") 
F.pred <- d.and.v.2 %>% 
  mutate(F.mean = invlogit(Mod.F.top.cf[1] + 
                           Mod.F.top.cf[2] * x.v + 
                           Mod.F.top.cf[3] * d.mean +
                           Mod.F.top.cf[4] * x.v * d.mean +
                           Mod.F.top.cf[5] * (d.mean^2) +
                           Mod.F.top.cf[6] * x.v * (d.mean^2)))
F.mean.df <- CData.s %>% 
  mutate(d.bin = cut_number(d.stand, n = 4),
         v.bin = cut_number(volume_t, n = 6)) %>% 
  group_by(v.bin, d.bin) %>% 
  summarise(F.mean = mean(did.flower),
            d.mean = mean(d.stand),
            v.mean = mean(volume_t))
F.pred %>% 
  ggplot() +
  geom_line(aes(x = x.v, y = F.mean, colour = as.factor(d.mean)), size = 3) +
  scale_colour_manual(values = c("darkorchid2", "firebrick2", "dodgerblue1", "darkorchid2",
                                 "navy", "firebrick2", "dodgerblue1", "navy")) +
  scale_y_continuous(position = "right") +
  geom_point(data = F.mean.df, aes(x = v.mean, y = F.mean, 
                                   colour = as.factor(d.bin)), size = 7) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        plot.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 38, margin = margin(t = 4, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 38, margin = margin(t = 0, r = 4, b = 0, l = 0)),
        axis.title.x = element_text(size = 43, margin = margin(t = 16, r = 0, b = 0, l = 6)),
        axis.title.y = element_text(size = 43, margin = margin(t = 80, r = 16, b = 0, l = 0)),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(8, "points")) -> Mod.F.top.vPlot





# ----- Growth as a function of volume and density (2 of 4) -----------------------------------------------

# Restore CData.s since we deleted stuff earlier
# Remove entries for which there is no recorded volume and log growth is NA
# Then standardise density and introduce unique transect identifier
CData.s <- CData %>% 
  drop_na(volume_t, logGR) %>% 
  mutate(d.stand = (weighted.dens - mean(weighted.dens, na.rm = TRUE)) / sd(weighted.dens, na.rm = TRUE),
         unique.transect = interaction(transect, site))

# Calculate range of data
v.range <- range(CData.s$volume_t)
d.range <- range(CData.s$d.stand)

# Calculate predictions for 4 "bins" of volume categories
v.cut <- CData.s %>% 
  mutate(v.bin = cut_number(volume_t, n = 4)) %>% 
  group_by(v.bin) %>% 
  summarise(v.mean = mean(volume_t))
d.and.v <- as.tibble(expand.grid(v.cut$v.mean, seq(d.range[1], d.range[2], 0.5)))
names(d.and.v) <- c("v.mean", "x.d") 
G.pred <- d.and.v %>% 
  mutate(G.mean = Mod.G.top.cf[1] + 
                  Mod.G.top.cf[2] * v.mean + 
                  Mod.G.top.cf[3] * x.d +
                  Mod.G.top.cf[4] * v.mean * x.d +
                  Mod.G.top.cf[5] * (x.d^2)+
                  Mod.G.top.cf[6] * v.mean * (x.d^2))

# Get binned means of the data with respect to volume and density
G.mean.df <- CData.s %>% 
  mutate(v.bin = cut_number(volume_t, n = 4),
         d.bin = cut_number(d.stand, n = 6)) %>% 
  group_by(d.bin, v.bin) %>% 
  summarise(G.mean = mean(logGR),
            d.mean = mean(d.stand),
            v.mean = mean(volume_t))

# Generate visualisation; growth ratio as function of density, with 4 volume bins
G.pred %>% 
  ggplot() +
  geom_line(aes(x = x.d, y = G.mean, colour = as.factor(v.mean)), size = 3) +
  scale_colour_manual(values = c("darkorchid2", "dodgerblue1", "navy", "firebrick2",
                                 "darkorchid2", "dodgerblue1", "navy", "firebrick2")) +
  scale_y_continuous(breaks = seq(0, 1, length.out = 5), limits = c(-0.1, 1)) +
  geom_point(data = G.mean.df, aes(x = d.mean, y = G.mean, 
                                   colour = as.factor(v.bin)), size = 7) +
  labs(x = NULL, y = "Log Annual Growth Ratio") +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        plot.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 38, margin = margin(t = 4, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 38, margin = margin(t = 0, r = 4, b = 0, l = 0)),
        axis.title.x = element_text(size = 43, margin = margin(t = 16, r = 0, b = 0, l = 6)),
        axis.title.y = element_text(size = 43, margin = margin(t = 80, r = 16, b = 0, l = 0)),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(8, "points")) -> Mod.G.top.dPlot

# Do same process as previous 3 steps, but switch density and volume
# This produces a graph of annual growth ratio as a function of volume, with 4 density bins
d.cut <- CData.s %>% 
  mutate(d.bin = cut_number(d.stand, n = 4)) %>% 
  group_by(d.bin) %>% 
  summarise(d.mean = mean(d.stand))
d.and.v.2 <- as.tibble(expand.grid(d.cut$d.mean, seq(v.range[1], v.range[2], 0.5)))
names(d.and.v.2) <- c("d.mean", "x.v") 
G.pred <- d.and.v.2 %>% 
  mutate(G.mean = Mod.G.top.cf[1] + 
                  Mod.G.top.cf[2] * x.v + 
                  Mod.G.top.cf[3] * d.mean +
                  Mod.G.top.cf[4] * x.v * d.mean +
                  Mod.G.top.cf[5] * (d.mean^2) +
                  Mod.G.top.cf[6] * x.v * (d.mean^2))
G.mean.df <- CData.s %>% 
  mutate(d.bin = cut_number(d.stand, n = 4),
         v.bin = cut_number(volume_t, n = 6)) %>% 
  group_by(v.bin, d.bin) %>% 
  summarise(G.mean = mean(logGR),
            d.mean = mean(d.stand),
            v.mean = mean(volume_t))
G.pred %>% 
  ggplot() +
  geom_line(aes(x = x.v, y = G.mean, colour = as.factor(d.mean)), size = 3) +
  scale_colour_manual(values = c("darkorchid2", "firebrick2", "dodgerblue1", "darkorchid2",
                                 "navy", "firebrick2", "dodgerblue1", "navy")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01), position = "right",
                     breaks = seq(0, 1, length.out = 5), limits = c(-0.05, 1)) +
  geom_point(data = G.mean.df, aes(x = v.mean, y = G.mean, 
                                   colour = as.factor(d.bin)), size = 7) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        plot.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 38, margin = margin(t = 4, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 38, margin = margin(t = 0, r = 4, b = 0, l = 0)),
        axis.title.x = element_text(size = 43, margin = margin(t = 16, r = 0, b = 0, l = 6)),
        axis.title.y = element_text(size = 43, margin = margin(t = 80, r = 16, b = 0, l = 0)),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(8, "points")) -> Mod.G.top.vPlot





# ----- Reproduction as a function of volume and density (3 of 4) -----------------------------------------

# Restore CData.s since we deleted stuff earlier
# Remove entries for which there is no recorded volume or that did not flower
# Then standardise density and introduce unique transect identifier
CData.s <- CData %>%
  drop_na(volume_t) %>%
  filter(did.flower == 1) %>% 
  mutate(d.stand = (weighted.dens - mean(weighted.dens, na.rm = TRUE)) / sd(weighted.dens, na.rm = TRUE),
         unique.transect = interaction(transect, site))

# Calculate range of data
v.range <- range(CData.s$volume_t)
d.range <- range(CData.s$d.stand)

# Calculate predictions for 4 "bins" of volume categories
v.cut <- CData.s %>% 
  mutate(v.bin = cut_number(volume_t, n = 4)) %>% 
  group_by(v.bin) %>% 
  summarise(v.mean = mean(volume_t))
d.and.v <- as.tibble(expand.grid(v.cut$v.mean, seq(d.range[1], d.range[2], 0.5)))
names(d.and.v) <- c("v.mean", "x.d") 
R.pred <- d.and.v %>% 
  mutate(R.mean = Mod.R.top.cf[1] + 
                  Mod.R.top.cf[2] * v.mean + 
                  Mod.R.top.cf[3] * x.d +
                  Mod.R.top.cf[4] * v.mean * x.d +
                  Mod.R.top.cf[5] * (x.d^2) +
                  Mod.R.top.cf[6] * v.mean * (x.d^2))

# Get binned means of the data with respect to volume and density
R.mean.df <- CData.s %>% 
  mutate(v.bin = cut_number(volume_t, n = 4),
         d.bin = cut_number(d.stand, n = 6)) %>% 
  group_by(d.bin, v.bin) %>% 
  summarise(R.mean = mean(logTR1),
            d.mean = mean(d.stand),
            v.mean = mean(volume_t))

# Generate visualisation; reproductive count as function of density, with 4 volume bins
R.pred %>% 
  ggplot() +
  geom_line(aes(x = x.d, y = R.mean, colour = as.factor(v.mean)), size = 3) +
  scale_colour_manual(values = c("darkorchid2", "dodgerblue1", "navy", "firebrick2",
                                 "firebrick2", "darkorchid2", "dodgerblue1", "navy")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  geom_point(data = R.mean.df, aes(x = d.mean, y = R.mean, 
                                   colour = as.factor(v.bin)), size = 7) +
  labs(x = "Standardised Weighted Density", y = "Log Reproductive Structure Count") +
  xlim(-1.5, 3) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        plot.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 38, margin = margin(t = 4, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 38, margin = margin(t = 0, r = 4, b = 0, l = 0)),
        axis.title.x = element_text(size = 43, margin = margin(t = 16, r = 0, b = 0, l = 6)),
        axis.title.y = element_text(size = 43, margin = margin(t = 80, r = 16, b = 0, l = 0)),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(8, "points")) -> Mod.R.top.dPlot

# Do same process as previous 3 steps, but switch density and volume
# This produces a graph of reproductive count as a function of volume, with 4 density bins
d.cut <- CData.s %>% 
  mutate(d.bin = cut_number(d.stand, n = 4)) %>% 
  group_by(d.bin) %>% 
  summarise(d.mean = mean(d.stand))
d.and.v.2 <- as.tibble(expand.grid(d.cut$d.mean, seq(v.range[1], v.range[2], 0.5)))
names(d.and.v.2) <- c("d.mean", "x.v") 
R.pred <- d.and.v.2 %>% 
  mutate(R.mean = Mod.R.top.cf[1] + 
                  Mod.R.top.cf[2] * x.v + 
                  Mod.R.top.cf[3] * d.mean +
                  Mod.R.top.cf[4] * x.v * d.mean +
                  Mod.R.top.cf[5] * (d.mean^2) +
                  Mod.R.top.cf[6] * x.v * (d.mean^2))
R.mean.df <- CData.s %>% 
  mutate(d.bin = cut_number(d.stand, n = 4),
         v.bin = cut_number(volume_t, n = 6)) %>% 
  group_by(v.bin, d.bin) %>% 
  summarise(R.mean = mean(logTR1),
            d.mean = mean(d.stand),
            v.mean = mean(volume_t))
R.pred %>% 
  ggplot() +
  geom_line(aes(x = x.v, y = R.mean, colour = as.factor(d.mean)), size = 3) +
  scale_colour_manual(values = c("darkorchid2", "firebrick2", "dodgerblue1", "darkorchid2",
                                 "navy", "firebrick2", "dodgerblue1", "navy")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01), position = "right",
                     breaks = seq(-3, 9, length.out = 5), limits = c(-3, 9)) +
  scale_x_continuous(breaks = c(5, 8, 12, 15)) +
  geom_point(data = R.mean.df, aes(x = v.mean, y = R.mean, 
                                   colour = as.factor(d.bin)), size = 7) +
  labs(x = "Log-Volume", y = NULL) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        plot.background = element_rect(fill = "white"),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 38, margin = margin(t = 4, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 38, margin = margin(t = 0, r = 4, b = 0, l = 0)),
        axis.title.x = element_text(size = 43, margin = margin(t = 16, r = 0, b = 0, l = 6)),
        axis.title.y = element_text(size = 43, margin = margin(t = 80, r = 16, b = 0, l = 0)),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(8, "points")) -> Mod.R.top.vPlot





# ----- Combine figures from the 3 sections above (4 of 4) ------------------------------------------------

# Prepare graphics device
jpeg(filename = "Figure 3.jpeg", width = 2600, height = 2950, units = "px")

# Create blank page
grid.newpage()
plot.new()

# Set grid layout and activate it
gly <- grid.layout(2950, 2600)
pushViewport(viewport(layout = gly))

# Place graphs in a grid
print(Mod.F.top.dPlot, vp = viewport(layout.pos.row = 1:950, layout.pos.col = 1:1300))
print(Mod.F.top.vPlot, vp = viewport(layout.pos.row = 1:950, layout.pos.col = 1300:2600))
print(Mod.G.top.dPlot, vp = viewport(layout.pos.row = 1000:1950, layout.pos.col = 1:1300))
print(Mod.G.top.vPlot, vp = viewport(layout.pos.row = 1000:1950, layout.pos.col = 1300:2600))
print(Mod.R.top.dPlot, vp = viewport(layout.pos.row = 2000:2950, layout.pos.col = 1:1300))
print(Mod.R.top.vPlot, vp = viewport(layout.pos.row = 2000:2950, layout.pos.col = 1300:2600))

# Create legend
grid.rect(vp = viewport(layout.pos.row = 190:260, layout.pos.col = 1400:1470), gp = gpar(fill = "navy"))
grid.rect(vp = viewport(layout.pos.row = 290:360, layout.pos.col = 1400:1470), gp = gpar(fill = "dodgerblue1"))
grid.rect(vp = viewport(layout.pos.row = 390:460, layout.pos.col = 1400:1470), gp = gpar(fill = "darkorchid2"))
grid.rect(vp = viewport(layout.pos.row = 490:560, layout.pos.col = 1400:1470), gp = gpar(fill = "firebrick2"))
grid.text(label = c("[75, 100]", "[50, 75)", "[25, 50)", "[0, 25)"),
          x = rep(0.580, 4), y = c(0.925, 0.892, 0.858, 0.823), just = "left", gp = gpar(fontsize = 37))
grid.text(label = c("Volume percentile (left graphs)", "Density percentile (right graphs)"),
          x = rep(0.537, 2), y = c(0.970, 0.952), just = "left", gp = gpar(fontsize = 40))

# Deactivate grid layout; finalise graphics save
popViewport()
dev.off()

# Clean up mess from global environment
remove(d.and.v, d.and.v.2, d.cut, F.mean.df, F.pred, G.mean.df, G.pred, R.mean.df, R.pred,
       v.cut, d.range, v.range, Mod.F.top.dPlot, Mod.F.top.vPlot, Mod.G.top.dPlot, 
       Mod.G.top.vPlot, Mod.R.top.dPlot, Mod.R.top.vPlot, gly)






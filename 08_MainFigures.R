# ----- Plot wavespeeds and growth rates ------------------------------------------------------------------

# Range of s values for plotting wavespeeds
s.values <- c(0.0001, 0.0005, 0.001, 0.005, seq(0.01, 2, length.out = 196))

# Find wavespeed minima
c.min <- min(c.values)
c.min.2 <- min(c.values.2)

# Find values of s at which wavespeed minima occur
s.c.min <- s.values[as.numeric(match(c.min, c.values))]
s.c.min.2 <- s.values[as.numeric(match(c.min.2, c.values.2))]

# Prepare graphics device
jpeg(filename = "Figure 1.jpeg", width = 2600, height = 900, units = "px")

# Create blank page
grid.newpage()
plot.new()

# Set grid layout and activate it
gly <- grid.layout(900, 2600)
pushViewport(viewport(layout = gly))

# Setup to plot wavespeeds
pushViewport(viewport(layout = gly, layout.pos.row = 1:900, layout.pos.col = 1:1300))
par(fig = gridFIG())
par(new = TRUE)
par(mar = c(7, 14, 2.1, 2.1))

# Plot wavespeed as function of shape parameter for both scenarios
plot(s.values, c.values, 
     ylim = c(0, 0.03), xlim = c(0, 1.5), type = "l", col = "yellow3", lwd = 6,
     mgp = c(5, 2, 0), cex.lab = 4, cex.axis = 3, axes = FALSE, ann = FALSE)
axis(1, at = seq(0, 1.5, length.out = 4), cex.axis = 3, mgp = c(5, 2, 0))
axis(2, at = seq(0, 0.03, length.out = 7), cex.axis = 3, mgp = c(5, 2, 0), las = 1)
box()
mtext("Wave Shape Parameter", side = 1, cex = 4, line = 5)
mtext("Wave Speed (m/yr)", side = 2, cex = 4, line = 10)
lines(s.values, c.values.2, 
      col = "green4", lwd = 6)
segments(x0 = c(-0.2, -0.2, s.c.min, s.c.min.2), x1 = c(s.c.min, s.c.min.2, s.c.min, s.c.min.2),
         y0 = c(c.min, c.min.2, -0.05, -0.05), y1 = c(c.min, c.min.2, c.min, c.min.2), lty = 2)
text(labels = c(paste0("~ ", round(c.min.2, digits = 4), " m/yr"),
                paste0("~ ", round(c.min, digits = 4), " m/yr")),
     x = rep(0.110, 2), y = c(c.min.2, c.min) + 0.0007, cex = 2.5)
popViewport()

# Setup to plot growth rates
pushViewport(viewport(layout.pos.row = 1:900, layout.pos.col = 1300:2600))
par(fig = gridFIG())
par(new = TRUE)
par(mar = c(7, 2.1, 2.1, 14))

# Plot population growth rates
plot(seq(-1.3, max(CData.s$d.stand), length.out = 100), lambda.2, 
     type = "l", lwd = 6, col = "green4", cex.lab = 4, cex.axis = 2, mgp = c(5, 2, 0),
     xlim = c(-1.15, 2.5), ylim = c(1, 1.012), axes = FALSE, ann = FALSE)
axis(1, cex.axis = 3, mgp = c(5, 2, 0))
axis(4, at = seq(1, 1.012, length.out = 7), cex.axis = 3, mgp = c(5, 2, 0), las = 1)
box()
mtext("Standardised Weighted Density", side = 1, cex = 4, line = 5)
mtext("Growth Rate", side = 4, cex = 4, line = 11)
lines(seq(-1.3, max(CData.s$d.stand), length.out = 100), lambda, 
      lwd = 6, col = "yellow3")
popViewport()

# Create legend
grid.rect(vp = viewport(layout.pos.row = 80:130, layout.pos.col = 1900:1950), gp = gpar(fill = "green4"))
grid.rect(vp = viewport(layout.pos.row = 150:200, layout.pos.col = 1900:1950), gp = gpar(fill = "yellow3"))
grid.text(label = c("Higher seedling survival", "Normal conditions"),
          x = rep(0.759, 2), y = c(0.888, 0.810), just = "left", gp = gpar(fontsize = 31))

# Deactivate grid layout; finalise graphics save
popViewport()
dev.off()

# Clean up mess from global environment
remove(s.values, c.min, c.min.2, s.c.min, s.c.min.2, gly)





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


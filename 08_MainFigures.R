# create vital rate figure (needs script 05 loaded) ------------------------------------------------

if(boot.switch == FALSE){
  
  size_breaks <- 4
  density_breaks <- 5
  LATR_cols <- wes_palette("Zissou1", size_breaks, type = "continuous")
  
  pdf("Manuscript/Figures/vital_rates.pdf",useDingbats = F,height=9,width=7)
  par(mfrow=c(3,2),mar=c(5,5,1,1))
  ## survival
  LATR_surv_dat %>% 
    filter(transplant==F) %>% 
    mutate(size_bin = as.numeric(cut(log_volume_t,breaks = size_breaks)),
           density_bin = as.numeric(cut(weighted.dens,breaks = density_breaks))) %>% 
    group_by(size_bin,density_bin) %>% 
    summarise(mean_size = mean(log_volume_t),
              mean_dens = mean(weighted.dens),
              mean_surv = mean(survival_t1),
              bin_n = n()) -> LATR_surv_plot_transF
  LATR_surv_plot_transF$size_col <- LATR_cols[LATR_surv_plot_transF$size_bin]
  
  mean_size <- LATR_surv_plot_transF %>% group_by(size_bin) %>% summarise(mean_size=mean(mean_size))
  surv_predict <- expand.grid(
    weighted.dens = seq(min(LATR_surv_dat$weighted.dens),max(LATR_surv_dat$weighted.dens),1),
    size_bin = mean_size$size_bin,
    transplant=F,
    unique.transect = LATR_flow_dat$unique.transect[1]
  )
  ## associate sizes with bin numbers for plotting
  surv_predict$log_volume_t<-mean_size$mean_size[surv_predict$size_bin]
  surv_predict$pred <- predict.gam(LATR_surv_best, newdata = surv_predict, type = "response", exclude = "s(unique.transect)")
  surv_predict$size_col <- LATR_cols[surv_predict$size_bin]
  
  plot(LATR_surv_plot_transF$mean_dens,LATR_surv_plot_transF$mean_surv,col=alpha(LATR_surv_plot_transF$size_col,0.5),pch=16,
       ylim=c(0,1),xlab="Weighted density",ylab="Survival probability",cex.lab=1.4,
       cex=LATR_surv_plot_transF$bin_n/max(LATR_surv_plot_transF$bin_n)*2+1,
       xlim=c(min(LATR_surv_dat$weighted.dens),max(LATR_surv_dat$weighted.dens)))
  points(surv_predict$weighted.dens,surv_predict$pred,col=surv_predict$size_col,pch=16,cex=0.75)
  title("A",adj=0)
  
  LATR_surv_dat %>% 
    filter(transplant==T) %>% 
    mutate(density_bin = as.numeric(cut(weighted.dens,breaks = density_breaks))) %>% 
    group_by(density_bin) %>% 
    summarise(mean_dens = mean(weighted.dens),
              mean_surv = mean(survival_t1),
              bin_n = n()) -> LATR_surv_plot_transT
  surv_trans_predict <- expand.grid(
    weighted.dens = seq(min(LATR_surv_dat$weighted.dens),max(LATR_surv_dat$weighted.dens),1),
    log_volume_t = mean(LATR_surv_dat$log_volume_t[LATR_surv_dat$transplant==T]),
    transplant=T,
    unique.transect = LATR_flow_dat$unique.transect[1]
  )
  surv_trans_predict$pred <- predict.gam(LATR_surv_best, newdata = surv_trans_predict, type = "response", exclude = "s(unique.transect)")
  
  points(LATR_surv_plot_transT$mean_dens,LATR_surv_plot_transT$mean_surv,col=alpha("black",0.5),pch=16,
         cex.lab=1.4,cex=LATR_surv_plot_transT$bin_n/max(LATR_surv_plot_transT$bin_n)*2+1)
  lines(surv_trans_predict$weighted.dens,surv_trans_predict$pred,lwd=3)
  #points(LATR_surv_dat$weighted.dens[LATR_surv_dat$transplant==T],LATR_surv_dat$survival_t1[LATR_surv_dat$transplant==T],
  #       col=alph)
  text(30,0.1,"Transplants",font=3)
  legend("right",legend=c("smallest"," "," ","largest"),lwd=2,
         col=LATR_cols,bg="white")
  
  ## growth
  ## starting with raw data visualization, binning sizes and densities
  ## this plot freaked me out because it clearly shows positive density
  ## dependence at the smallest size -- no matter how you slice up sizes
  LATR_grow %>% 
    mutate(size_bin = as.numeric(cut(log_volume_t,breaks = size_breaks)),
           density_bin = as.numeric(cut(weighted.dens,breaks = density_breaks))) %>% 
    group_by(size_bin,density_bin) %>% 
    summarise(mean_size = mean(log_volume_t),
              mean_dens = mean(weighted.dens),
              mean_sizet1 = mean(log_volume_t1),
              sd_sizet1 = sd(log_volume_t1),
              bin_n = n()) -> LATR_grow_plot
  LATR_grow_plot$size_col <- LATR_cols[LATR_grow_plot$size_bin]
  
  #plot(LATR_grow_plot$mean_dens,LATR_grow_plot$mean_sizet1,col=LATR_grow_plot$size_bin,pch=16,cex=2)
  ## however, I realized that the bins actually differ in initial size, and
  ## this can explain the apparent positive DD, here plotting the change in size to control for size differences
  #plot(LATR_grow_plot$mean_dens,(LATR_grow_plot$mean_sizet1-LATR_grow_plot$mean_size),col=LATR_grow_plot$size_bin,pch=16,cex=2)
  
  ## let's avoid binning the growth data for this read
  ## create dummy data frame for gam prediction
  mean_size <- LATR_grow %>% 
    mutate(size_bin = as.numeric(cut(log_volume_t,breaks = size_breaks))) %>% 
    group_by(size_bin) %>% 
    summarise(mean_size=mean(log_volume_t))
  grow_predict <- expand.grid(
    weighted.dens = seq(min(LATR_grow$weighted.dens),max(LATR_grow$weighted.dens),1),
    size_bin = mean_size$size_bin,
    unique.transect = LATR_grow$unique.transect[1]
  )
  grow_predict$log_volume_t<-mean_size$mean_size[grow_predict$size_bin]
  grow_predict[,c("pred_mu","pred_sigma")] <- predict.gam(LATR_grow_best, newdata = grow_predict, type = "response", exclude = "s(unique.transect)")
  LATR_grow$size_col <- LATR_cols[as.numeric(cut(LATR_grow$log_volume_t,breaks = size_breaks))]
  grow_predict$size_col <- LATR_cols[grow_predict$size_bin]
  
  ## growth mean and SD inset
  plot(LATR_grow$weighted.dens,LATR_grow$log_volume_t1,pch=16,cex=0.5,cex.lab=1.4,
       col=alpha(LATR_grow$size_col,0.5),xlab="Weighted density",ylab=expression(paste(Size[t+1]," (log ",cm^3,")")))
  points(grow_predict$weighted.dens,grow_predict$pred_mu,col=grow_predict$size_col,pch=16,cex=0.75)
  rect(100, -2.5, 220, 6, col="white",lwd=0.5)
  plotInset(100, -1, 220, 6,
            expr=plot(grow_predict$weighted.dens,1/(grow_predict$pred_sigma^2),
                      col=grow_predict$size_col,pch=16,cex=0.25,xlab="",ylab=expression(paste("SD(",Size[t+1],")")),
                      cex.axis=0.5,mgp=c(3/2, 1/2, 0)),
            mar=c(0, 3, 0, 0))
  title("B",adj=0)
  
  ## flowering
  LATR_flow_dat %>% 
    mutate(size_bin = as.numeric(cut(log_volume_t,breaks = size_breaks)),
           density_bin = as.numeric(cut(weighted.dens,breaks = density_breaks))) %>% 
    group_by(size_bin,density_bin) %>% 
    summarise(mean_size = mean(log_volume_t),
              mean_dens = mean(weighted.dens),
              mean_flow = mean(total.reproduction_t > 0),
              bin_n = n()) -> LATR_flow_plot
  LATR_flow_plot$size_col <- LATR_cols[LATR_flow_plot$size_bin]
  
  mean_size <- LATR_flow_plot %>% group_by(size_bin) %>% summarise(mean_size=mean(mean_size))
  flow_predict <- expand.grid(
    weighted.dens = seq(min(LATR_flow_dat$weighted.dens),max(LATR_flow_dat$weighted.dens),1),
    size_bin = mean_size$size_bin,
    unique.transect = LATR_flow_dat$unique.transect[1]
  )
  ## associate sizes with bin numbers for plotting
  flow_predict$log_volume_t<-mean_size$mean_size[flow_predict$size_bin]
  flow_predict$pred <- predict.gam(LATR_flower_best, newdata = flow_predict, type = "response", exclude = "s(unique.transect)")
  flow_predict$size_col <- LATR_cols[flow_predict$size_bin]
  
  plot(LATR_flow_plot$mean_dens,LATR_flow_plot$mean_flow,col=alpha(LATR_flow_plot$size_col,0.5),pch=16,
       ylim=c(0,1),xlab="Weighted density",ylab="Flowering probability",
       cex=LATR_flow_plot$bin_n/max(LATR_flow_plot$bin_n)*2+1,cex.lab=1.4,
       xlim=c(min(LATR_flow_dat$weighted.dens),max(LATR_flow_dat$weighted.dens)))
  points(flow_predict$weighted.dens,flow_predict$pred,col=flow_predict$size_col,pch=16,cex=0.75)
  title("C",adj=0)
  
  ## fruits
  mean_size <- LATR_fruits_dat %>% 
    mutate(size_bin=as.numeric(cut(log_volume_t,breaks = size_breaks))) %>% 
    group_by(size_bin) %>% 
    summarise(mean_size=mean(log_volume_t))
  fruits_predict <- expand.grid(
    weighted.dens = seq(min(LATR_fruits_dat$weighted.dens),max(LATR_fruits_dat$weighted.dens),1),
    size_bin = mean_size$size_bin,
    unique.transect = LATR_fruits_dat$unique.transect[1]
  )
  ## associate sizes with bin numbers for plotting
  LATR_fruits_dat$size_col <- LATR_cols[as.numeric(cut(LATR_fruits_dat$log_volume_t,breaks = size_breaks))]
  ## note that these size breaks are *different* breaks for fruits than for the other vital rate (bc this is the flowering subset)
  fruits_predict$log_volume_t<-mean_size$mean_size[fruits_predict$size_bin]
  fruits_predict$pred <- predict.gam(LATR_fruits_best, newdata = fruits_predict, type = "response", exclude = "s(unique.transect)")
  fruits_predict$size_col <- LATR_cols[fruits_predict$size_bin]
  
  plot(LATR_fruits_dat$weighted.dens,LATR_fruits_dat$total.reproduction_t,pch=16,cex=0.5,ylim=c(0,4000),cex.lab=1.4,
       col=alpha(LATR_fruits_dat$size_col,0.5),xlab="Weighted density",ylab="Flowers and fruits")
  points(fruits_predict$weighted.dens,fruits_predict$pred,col=fruits_predict$size_col,pch=16,cex=0.75)
  title("D",adj=0)
  
  ## recruitment
  recruit_predict <- expand.grid(
    weighted.dens = seq(min(LATR_recruitment$weighted.dens),max(LATR_recruitment$weighted.dens),1),
    unique.transect = LATR_fruits_dat$unique.transect[1]
  )
  recruit_predict$pred <- predict.gam(LATR_recruit_best, newdata = recruit_predict, type = "response", exclude = "s(unique.transect)")
  
  plot(LATR_recruitment$weighted.dens,LATR_recruitment$recruits/LATR_recruitment$total_seeds,cex.lab=1.4,
       col=alpha("black",0.5),pch=1,xlab="Weighted density",ylab="Per-seed recruitment probability",ylim=c(0,0.001))
  lines(recruit_predict$weighted.dens,recruit_predict$pred,lwd=3)
  title("E",adj=0)
  
  plot(LATR_recruit_size$weighted.dens,LATR_recruit_size$log_volume,col=alpha("black",0.5),cex.lab=1.4,
       pch=16,xlab="Weighted density",ylab=expression(paste("Recruit size (log ",cm^3,")")))
  lines(LATR_recruit_size$weighted.dens,LATR_recruit_size$pred[,1],lwd=3)
  title("F",adj=0)
  
  dev.off()
  
}

# ----- Plot bootstrapped wavespeeds ----------------------------------------------------------------------

# Prepare graphics device
jpeg(filename = "Figure 1.jpeg", width = 750, height = 500, units = "px")

# Create blank page
grid.newpage()
plot.new()

# Set grid layout and activate it
gly <- grid.layout(500, 750)
pushViewport(viewport(layout = gly))

# Setup to plot data
pushViewport(viewport(layout = gly, layout.pos.row = 1:500, layout.pos.col = 1:750))
par(fig = gridFIG())
par(new = TRUE)
par(mar = c(5, 8, 4, 2))

# Plot PDF of bootstrapped wavespeeds
plot(density(boot.cv1, from = 0, to = 0.25), lwd = 3, xlim = c(0, 0.25), ylim = c(0, 25), 
     axes = FALSE, ann = FALSE, zero.line = FALSE)
axis(1, at = seq(0, 0.25, length.out = 6), cex.axis = 1.5, mgp = c(1, 1, 0))
axis(2, at = seq(0, 25, length.out = 6), cex.axis = 1.5, mgp = c(1, 1, 0), las = 1)
mtext("Wave Speed (m/yr)", side = 1, cex = 2, line = 3.5)
mtext("Probability Density", side = 2, cex = 2, line = 4)
box()
popViewport()

# Deactivate grid layout; finalise graphics save
popViewport()
dev.off()





# ----- Plot growth rate as function of density -----------------------------------------------------------

# get bootstrapped results
boot.lambda <- read.csv("BootLambda.csv")
# compute mean result
mean.lambda <- LambdaD()

plot(boot.lambda$density,boot.lambda[,3],type="n",xlab="Weighted density",
     ylab=expression(paste(lambda)),cex.lab=1.2,ylim=c(1,1.06))
for(i in 3:dim(boot.lambda)[2]){
  lines(boot.lambda$density,boot.lambda[,i],col=alpha("black",0.15))
}
lines(boot.lambda$density,mean.lambda,lwd=3,col="red")

boot.lambda.stats <- data.frame(boot.lambda.df[, 1],
                               apply(X = boot.lambda.df[, -1], MARGIN = 1, FUN = mean),
                               apply(X = boot.lambda.df[, -1], MARGIN = 1, FUN = quantile, probs = 0.025),
                               apply(X = boot.lambda.df[, -1], MARGIN = 1, FUN = quantile, probs = 0.975)) 

# Prepare graphics device
jpeg(filename = "Figure 2.jpeg", width = 750, height = 500, units = "px")

# Create blank page
grid.newpage()
plot.new()

# Set grid layout and activate it
gly <- grid.layout(500, 750)
pushViewport(viewport(layout = gly))

# Setup to plot data
pushViewport(viewport(layout = gly, layout.pos.row = 1:500, layout.pos.col = 1:750))
par(fig = gridFIG())
par(new = TRUE)
par(mar = c(5, 8, 4, 2))

# Note: choose only ONE of the following four plots before deactivating grid and finalising graphics save

# [1] Plot lambda as function of density, with individual curves
plot(x = boot.lambda.stats[, 1], y = boot.lambda.stats[, 2], type = "l", lwd = 3,
     ylim = c(1, 1.06), xlim = c(0, 200), 
     axes = FALSE, ann = FALSE)
for(i in 2:length(boot.lambda)){
  lines(x = boot.lambda[[1]], y = boot.lambda[[i]], col = "grey")}
lines(x = boot.lambda.stats[, 1], y = boot.lambda.stats[, 2], lwd = 3)
axis(1, at = seq(0, 200, length.out = 5), cex.axis = 1.5, mgp = c(1, 1, 0))
axis(2, at = seq(1, 1.06, length.out = 7), cex.axis = 1.5, mgp = c(1, 1, 0), las = 1)
mtext("Weighted Density", side = 1, cex = 2, line = 3.5)
mtext("Population Growth Rate", side = 2, cex = 2, line = 5)
box()

# [2] Plot lambda as function of density, with individual curves (w/ transparency)
plot(x = boot.lambda.stats[, 1], y = boot.lambda.stats[, 2], type = "l", lwd = 3,
     ylim = c(1, 1.06), xlim = c(0, 200), 
     axes = FALSE, ann = FALSE)
for(i in 2:length(boot.lambda)){
  lines(x = boot.lambda[[1]], y = boot.lambda[[i]], col = rgb(red = 0, blue = 0, green = 0, alpha = 0.2))}
lines(x = boot.lambda.stats[, 1], y = boot.lambda.stats[, 2], lwd = 3, col = "red")
axis(1, at = seq(0, 200, length.out = 5), cex.axis = 1.5, mgp = c(1, 1, 0))
axis(2, at = seq(1, 1.06, length.out = 7), cex.axis = 1.5, mgp = c(1, 1, 0), las = 1)
mtext("Weighted Density", side = 1, cex = 2, line = 3.5)
mtext("Population Growth Rate", side = 2, cex = 2, line = 5)
box()

# [3] Plot lambda as function of density, with individual curves
# Curves coloured based on starting distance from mean
lmeans <- apply(X = boot.lambda.df[, -1], MARGIN = 1, FUN = mean)
colours <- colorRampPalette(c("gray15", "gray90"))
boot.lambda.df2 <- data.frame(t(boot.lambda.df))[-1, ]
data.frame(cbind(boot.lambda.df2, 1:(ncol(boot.lambda.df) - 1),
                 as.numeric((boot.lambda.df[1, -1] - lmeans[1])^2))) %>% 
  arrange(.[, 27]) %>% 
  cbind(colours(ncol(boot.lambda.df) - 1)) %>% 
  arrange(.[, 26]) -> plotdata
plot(x = boot.lambda.stats[, 1], y = boot.lambda.stats[, 2], type = "l", lwd = 3,
     ylim = c(1, 1.06), xlim = c(0, 200), 
     axes = FALSE, ann = FALSE)
for(i in 2:length(boot.lambda)){
  lines(x = boot.lambda.df[[1]], y = plotdata[i, 1:25], col = as.character(plotdata[i, 28]))}
lines(x = boot.lambda.stats[, 1], y = boot.lambda.stats[, 2], lwd = 3)
axis(1, at = seq(0, 200, length.out = 5), cex.axis = 1.5, mgp = c(1, 1, 0))
axis(2, at = seq(1, 1.06, length.out = 7), cex.axis = 1.5, mgp = c(1, 1, 0), las = 1)
mtext("Weighted Density", side = 1, cex = 2, line = 3.5)
mtext("Population Growth Rate", side = 2, cex = 2, line = 5)
box()

# [4] Plot lambda as function of density, with 95% bootstrap interval
plot(x = boot.lambda.stats[, 1], y = boot.lambda.stats[, 2], type = "l", lwd = 3,
     ylim = c(1, 1.06), xlim = c(0, 200), 
     axes = FALSE, ann = FALSE)
polygon(x = c(boot.lambda.stats[, 1], rev(boot.lambda.stats[, 1])), 
        y = c(boot.lambda.stats[, 3], rev(boot.lambda.stats[, 4])),
        col = alpha("black", alpha = 0.2), border = NA)
axis(1, at = seq(0, 200, length.out = 5), cex.axis = 1.5, mgp = c(1, 1, 0))
axis(2, at = seq(1, 1.06, length.out = 7), cex.axis = 1.5, mgp = c(1, 1, 0), las = 1)
mtext("Weighted Density", side = 1, cex = 2, line = 3.5)
mtext("Population Growth Rate", side = 2, cex = 2, line = 5)
box()

# Deactivate grid layout; finalise graphics save
popViewport(2)
dev.off()

# Clean up mess from global environment
remove(gly, boot.lambda.df)





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
plot(density(na.omit(WALD.b.h(H = 0.6, elas = 0, seed = 8675309, reps = 10000, heights = 50)),
             from = 0, to = 10, bw = 0.05), lwd = 4, xlim = c(0, 2),
     ylim = c(0, 5), xlab = "Distance (m)", ylab = "Probability Density", axes = FALSE, ann = FALSE)
axis(1, at = seq(0, 2, length.out = 5), cex.axis = 2, mgp = c(5, 2, 0))
axis(2, at = seq(0, 5, length.out = 6), cex.axis = 2, mgp = c(5, 2, 0), las = 1)
box()
mtext("Distance (m)", side = 1, cex = 3, line = 5)
mtext("Probability Density", side = 2, cex = 3, line = 4)
for(i in 1:4){
  colours <- c("black", "blue", "forestgreen", "red", "orange")
  lines(density(na.omit(WALD.b.h(H = 0.6 + 0.2*i, elas = 0, seed = 8675309, reps = 10000, heights = 50)),
                from = 0, to = 10, bw = 0.05), 
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
mtext("Weighted Density", side = 1, cex = 2.5, line = 3.5)
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
mtext("Weighted Density", side = 1, cex = 2.5, line = 3.5)
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
mtext("Weighted Density", side = 1, cex = 2, line = 3.5)
mtext("Pr(Recruitment)", side = 2, cex = 2, line = 5)
lines(LATR_recruitment_line$weighted.dens, invlogit(LATR_recruitment_line$pred), col = "red")
box()
popViewport()

# Deactivate grid layout; finalise graphics save
popViewport()
dev.off()


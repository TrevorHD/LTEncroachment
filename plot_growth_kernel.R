source("C:/Users/tm9/Dropbox/github/IPM_size_transitions/Diagnostics.R")
## need to run the growth model in 05_CDataAnalysis.R

## compare to original (gaussian) gam parameter estimates
plot(LATR_beta,out$estimate[1:66])
abline(0,1)
## compare the expected value of the two models
plot(LATR_Xp[,1:(grow_sd_index-1)]%*%LATR_beta[1:(grow_sd_index-1)],
     LATR_Xp[,1:(grow_sd_index-1)]%*%out$estimate[1:(grow_sd_index-1)])
abline(0,1)
## I'm satisfied that the sgt can recover the same expected value as the gaussian gam()

# compare simulated and real data -----------------------------------------
n_sim <- 500
LATR_sim_NO<-matrix(NA,nrow=nrow(LATR_grow),ncol=n_sim)
grow_pred<-predict(LATR_grow_best,type="response",data=LATR_grow)
for(i in 1:n_sim){
  print(i)
  #LATR_sim_SGT[,i] <- rsgt(n = nrow(LATR_grow), 
  #                         mu = LATR_Xp[,1:(grow_sd_index-1)]%*%out$estimate[1:(grow_sd_index-1)], 
#                           sigma = exp(LATR_Xp[,grow_sd_index:gam_coef_length]%*%out$estimate[grow_sd_index:gam_coef_length]),
#                           lambda=-invlogit(out$estimate[(gam_coef_length+1)]+out$estimate[(gam_coef_length+2)]*LATR_grow$log_volume_t),
#                           p=exp(out$estimate[(gam_coef_length+3)]),
#                           q=exp(out$estimate[(gam_coef_length+4)]),
#                           mean.cent=T,
#                           var.adj=T)
  LATR_sim_NO[,i] <- rnorm(n = nrow(LATR_grow),mean=grow_pred[,1],sd=1/grow_pred[,2])  
}

n_bins = 6
alpha_scale = 0.7
LATR_moments <- LATR_grow %>% 
  arrange(log_volume_t) %>% 
  mutate(size_bin = cut_interval(log_volume_t,n=n_bins)) %>% 
  group_by(size_bin) %>% 
  summarise(mean_t1 = mean(log_volume_t1),
            sd_t1 = sd(log_volume_t1),
            skew_t1 = NPskewness(log_volume_t1),
            kurt_t1 = NPkurtosis(log_volume_t1),
            bin_mean = mean(log_volume_t),
            bin_n = n()) 

par(mfrow=c(2,2),mar=c(4,4,2,1),cex.axis=1.3,cex.lab=1.3,mgp=c(2,1,0),bty="l"); 
sim_bin_means=sim_moment_means=sim_moment_means_norm = matrix(NA,n_bins,n_sim); 
for(i in 1:n_sim){
  sim_moments <- bind_cols(LATR_grow,data.frame(sim=LATR_sim_NO[,i])) %>% 
    arrange(log_volume_t) %>% 
    mutate(size_bin = cut_interval(log_volume_t,n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_t1 = mean(sim),
              bin_mean = mean(log_volume_t))
  sim_bin_means[,i]=sim_moments$bin_mean; 
  sim_moment_means[,i]=sim_moments$mean_t1; 		  
}
matplot(LATR_moments$bin_mean, sim_moment_means,col=alpha("gray",0.5),pch=16,xlab="Mean size t0",ylab="Mean(Size t1)",cex=1,
        xlim=c(min(LATR_moments$bin_mean),max(LATR_moments$bin_mean)+0.4))
points(LATR_moments$bin_mean+0.2, LATR_moments$mean_t1,pch=16,lwd=2,col=alpha("red",alpha_scale),cex=1.6)
points(LATR_moments$bin_mean, apply(sim_moment_means,1,median),pch=1,lwd=2,col=alpha("black",alpha_scale),cex=1.6)
legend("topleft",legend=c("NO","Data"),
       col=c("gray","red"),pch=16,bty="n",cex=1.2,pt.lwd=2,pt.cex = 1.2) 
add_panel_label("a")

for(i in 1:n_sim){
  sim_moments <- bind_cols(LATR_grow,data.frame(sim=LATR_sim_NO[,i])) %>% 
    arrange(log_volume_t) %>% 
    mutate(size_bin = cut_interval(log_volume_t,n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_t1 = sd(sim),
              bin_mean = mean(log_volume_t))
  sim_bin_means[,i]=sim_moments$bin_mean; 
  sim_moment_means[,i]=sim_moments$mean_t1;	  
}
matplot(LATR_moments$bin_mean, sim_moment_means,col=alpha("gray",0.5),pch=16,xlab="Mean size t0",ylab="SD(Size t1)",cex=1,
        xlim=c(min(LATR_moments$bin_mean),max(LATR_moments$bin_mean)+0.4)) 
points(LATR_moments$bin_mean+0.2, LATR_moments$sd_t1,pch=16,lwd=2,col=alpha("red",alpha_scale),cex=1.6)
points(LATR_moments$bin_mean, apply(sim_moment_means,1,median),pch=1,lwd=2,col=alpha("black",alpha_scale),cex=1.6)
add_panel_label("b")

for(i in 1:n_sim){
  sim_moments <- bind_cols(LATR_grow,data.frame(sim=LATR_sim_NO[,i])) %>% 
    arrange(log_volume_t) %>% 
    mutate(size_bin = cut_interval(log_volume_t,n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_t1 = NPskewness(sim),
              bin_mean = mean(log_volume_t))
  sim_bin_means[,i]=sim_moments$bin_mean; 
  sim_moment_means[,i]=sim_moments$mean_t1;   
}
matplot(LATR_moments$bin_mean, sim_moment_means,col=alpha("gray",0.5),pch=16,xlab="Mean size t0",ylab="Skew(Size t1)",cex=1,
        xlim=c(min(LATR_moments$bin_mean),max(LATR_moments$bin_mean)+0.4))
points(LATR_moments$bin_mean+0.2, LATR_moments$skew_t1,pch=16,lwd=2,col=alpha("red",alpha_scale),cex=1.6)
points(LATR_moments$bin_mean, apply(sim_moment_means,1,median),pch=1,lwd=2,col=alpha("black",alpha_scale),cex=1.6)
add_panel_label("c")

for(i in 1:n_sim){
  sim_moments <- bind_cols(LATR_grow,data.frame(sim=LATR_sim_NO[,i])) %>% 
    arrange(log_volume_t) %>% 
    mutate(size_bin = cut_interval(log_volume_t,n=n_bins)) %>% 
    group_by(size_bin) %>% 
    summarise(mean_t1 = NPkurtosis(sim),
              bin_mean = mean(log_volume_t))
  sim_bin_means[,i]=sim_moments$bin_mean; 
  sim_moment_means[,i]=sim_moments$mean_t1;  
}
matplot(LATR_moments$bin_mean, sim_moment_means,col=alpha("gray",0.5),pch=16,xlab="Mean size t0",ylab="Kurtosis(Size t1)",cex=1,
        xlim=c(min(LATR_moments$bin_mean),max(LATR_moments$bin_mean)+0.4))
points(LATR_moments$bin_mean+0.2, LATR_moments$kurt_t1,pch=16,lwd=2,col=alpha("red",alpha_scale),cex=1.6)
points(LATR_moments$bin_mean, apply(sim_moment_means,1,median),pch=1,lwd=2,col=alpha("black",alpha_scale),cex=1.6)
add_panel_label("d")


# plot growth kernel ------------------------------------------------------
## checking that the fitted kernel gives a visually satisfying fit to the data
size_dim <- 200
size_dum <- seq(min(LATR_grow$log_volume_t1),max(LATR_grow$log_volume_t1),length.out = size_dim)
## problem here is that growth has a continuous covariate. to visualize this I will
## discretize 4 quantiles of the weighted density distribution
hist(LATR_grow$weighted.dens)
quantile(LATR_grow$weighted.dens,c(0.25,0.5,0.75))
bins.quantiles(LATR_grow$weighted.dens,target.bins = 4,max.breaks = 4)

growth_kernel <- matrix(NA,size_dim,size_dim)
for(i in 1:size_dim){
  mu_size <- fixef(CYIM_lmer_best)[1] + fixef(CYIM_lmer_best)[2] * size_dum[i]
  CYIM_lmer_best_kernel[,i] <- dnorm(size_dum,
                                     mean = mu_size,
                                     sd = exp(pars[[best_model]][1] + pars[[best_model]][2]*mu_size))
}

levelplot(CYIM_lmer_best_kernel,row.values = size_dum, column.values = size_dum,cuts=30,
          col.regions=rainbow(30),xlab="log Size t",ylab="log Size t+1",main="Gaussian, non-constant variance",
          panel = function(...) {
            panel.levelplot(...)
            grid.points(log(CYIM$vol_t), log(CYIM$vol_t1), pch = ".",gp = gpar(cex=3,col=alpha("black",0.5)))
          })
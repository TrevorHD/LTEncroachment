
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
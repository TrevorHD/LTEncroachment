##### Fix important parameters for SIPM -------------------------------------------------------------------

# Give dimension of the size vector (# bins of the approximating matrix)
TM.matdim <- 200

# Eviction extensions for upper and lower size limits
TM.lower.extension <- -8
TM.upper.extension <- 2

##### IPM functions ---------------------------------------------------------------------------------------

# Growth from size x to y at density d, using best GAM -- GAUSSIAN
TM.growth <- function(x, y, d, elas=0){
  xb = pmin(pmax(x, LATR_size_bounds$min_size), LATR_size_bounds$max_size)
  lpmat <- predict.gam(LATR_grow_best,
                       newdata = data.frame(weighted.dens = d, log_volume_t = xb, unique.transect = "1.FPS"),
                       type = "lpmatrix",
                       exclude = "s(unique.transect)")
  # Linear predictor for mean and log sigma 
  grow_mu <- lpmat[, 1:(grow_sd_index-1)] %*% coef(LATR_grow_best)[1:(grow_sd_index-1)]
  if(elas=="growth"){grow_mu=grow_mu*(1+pert)}
  grow_sigma <- exp(lpmat[, grow_sd_index:length(coef(LATR_grow_best))] %*% coef(LATR_grow_best)[grow_sd_index:length(coef(LATR_grow_best))])
  return(dnorm(y, mean = grow_mu, sd = grow_sigma))}

# Growth from size x to y at density d, using best GAM -- SGT!!
#TM.growth <- function(x, y, d){
#  xb = pmin(pmax(x, LATR_size_bounds$min_size), LATR_size_bounds$max_size)
#  lpmat <- predict.gam(LATR_grow_best,
#                       newdata = data.frame(weighted.dens = d, log_volume_t = xb, unique.transect = "1.FPS"),
#                       type = "lpmatrix",
#                       exclude = "s(unique.transect)")
#  # Linear predictor for mean, sigma, and lambda
#  grow_mu <- lpmat[, 1:(grow_sd_index-1)] %*% coef_grow_best[1:(grow_sd_index-1)]
#  grow_sigma <- exp(lpmat[, grow_sd_index:gam_coef_length] %*% coef_grow_best[grow_sd_index:gam_coef_length])
#  grow_lambda <- -invlogit(coef_grow_best[(gam_coef_length+1)]+coef_grow_best[(gam_coef_length+2)]*xb)
#  return(dsgt(x = y, 
#              mu=grow_mu,
#              sigma=grow_sigma,
#              lambda=grow_lambda,
#              p=exp(coef_grow_best[(gam_coef_length+3)]),
#              q=exp(coef_grow_best[(gam_coef_length+4)]),
#              mean.cent=T,
#              var.adj=T))}


# Survival of size x at density d using best GAM
# For nnaturally occuring plants (transplant = FALSE)
TM.survival <- function(x, d, elas=0){
  xb = pmin(pmax(x, LATR_size_bounds$min_size), LATR_size_bounds$max_size)
  lpmat <- predict.gam(LATR_surv_best,
                       newdata = data.frame(weighted.dens = d, log_volume_t = xb, transplant = FALSE,
                                            unique.transect = "1.FPS"),
                       type = "lpmatrix",
                       exclude = "s(unique.transect)")
  pred <- lpmat %*% coef(LATR_surv_best)
  if(elas=="survival"){pred=pred*(1+pert)}
  return(invlogit(pred))
  }

# Combined growth and survival at density d
TM.growsurv <- function(x, y, d, elas=0){
  TM.survival(x, d, elas) * TM.growth(x, y, d, elas)}

# Flowering at size x and density d using best GAM
TM.flower <- function(x, d, elas=0){
  xb = pmin(pmax(x, LATR_size_bounds$min_size), LATR_size_bounds$max_size)
  lpmat <- predict.gam(LATR_flower_best,
                       newdata = data.frame(weighted.dens = d, log_volume_t = xb, unique.transect = "1.FPS"),
                       type = "lpmatrix",
                       exclude = "s(unique.transect)")
  pred <- lpmat %*% coef(LATR_flower_best)
  if(elas=="flower"){pred=pred*(1+pert)}
  return(invlogit(pred))}

# Seed production (fruits * seeds/fruit) at size x and density d using best GAM
# Note: we assume 6 seeds per fruit, and best GAM is actually not density dependent
TM.seeds <- function(x, d, seeds.per.fruit = 5, elas=0){
  xb = pmin(pmax(x, LATR_size_bounds$min_size), LATR_size_bounds$max_size)
  lpmat <- predict.gam(LATR_fruits_best,
                       newdata = data.frame(weighted.dens = d, log_volume_t = xb, unique.transect = "1.FPS"),
                       type = "lpmatrix",
                       exclude = "s(unique.transect)")
  pred <- lpmat %*% coef(LATR_fruits_best)
  if(elas=="fertility"){pred=pred*(1+pert)}
  return(exp(pred*delta)*seeds.per.fruit)}

# Seed-to-Seedling recruitment probability at density d
TM.recruitment <- function(d,elas=0){
  lpmat <- predict.gam(LATR_recruit_best,
                       newdata = data.frame(weighted.dens = d, unique.transect = "1.FPS"),
                       type = "lpmatrix",
                       exclude = "s(unique.transect)")
  pred <- lpmat%*% coef(LATR_recruit_best)
  if(elas=="recruitment"){pred=pred*(1+pert)}
  return(invlogit(pred[[1]]*delta))}

# Recruit size distribution at size y
TM.recruitsize <- function(y,d,elas=0){
  lpmat <- predict.gam(LATR_recruitsize_best,
                       newdata = data.frame(weighted.dens = d, unique.transect = "1.FPS"),
                       type = "lpmatrix",
                       exclude = "s(unique.transect)")
  recruitsize_mu <- lpmat[, 1:(recruit_size_sd_index-1)] %*% coef(LATR_recruitsize_best)[1:(recruit_size_sd_index-1)]
  if(elas=="recruitsize"){recruitsize_mu=recruitsize_mu*(1+pert)}
  recruitsize_sigma <- exp(lpmat[, recruit_size_sd_index:recruit_size_coef_length] %*% coef(LATR_recruitsize_best)[recruit_size_sd_index:recruit_size_coef_length])
  return(dnorm(x = y, mean = recruitsize_mu*delta, sd = recruitsize_sigma))
  }

# Combined flowering, fertility, and recruitment
TM.fertrecruit <- function(x, y, d, elas=0){
  TM.flower(x, d, elas) * TM.seeds(x, d, elas) * TM.recruitment(d, elas) * TM.recruitsize(y,d,elas)}

# Put it all together; projection matrix is a function of weighted density (dens)
# We need a large lower extension because growth variance (gaussian) is greater for smaller plants
TransMatrix <- function(dens, ext.lower = TM.lower.extension, ext.upper = TM.upper.extension,
                        min.size = LATR_size_bounds$min_size, max.size = LATR_size_bounds$max_size,
                        mat.size = TM.matdim, elas=0){
  
  # Matrix size and size extensions (upper and lower integration limits)
  n <- mat.size
  L <- min.size + ext.lower
  U <- max.size + ext.upper
  
  # Bin size for n bins
  h <- (U - L)/n
  
  # Lower boundaries of bins 
  b <- L + c(0:n)*h
  
  # Bin midpoints
  y <- 0.5*(b[1:n] + b[2:(n + 1)])
  
  # Growth/Survival matrix
  Pmat <- t(outer(y, y, TM.growsurv, d = dens, elas = elas)) * h 
  
  # Fertility/Recruiment matrix
  Fmat <- t(outer(y, y, TM.fertrecruit, d = dens, elas = elas)) * h 
  
  # Put it all together
  IPMmat <- Pmat + Fmat
  
  #and transition matrix
  return(list(IPMmat = IPMmat, Fmat = Fmat, Pmat = Pmat, meshpts = y))}

##### Find minimum wave speed -----------------------------------------------------------------------------

# Function to calculate the minimum wavespeed across a range of s
Wavespeed <- function(n = TM.matdim, elas=0){
  
  # Fit equation to convert volume to height for dispersal kernel use
  LATR_full %>%
    select(max.ht_t, volume_t) %>% 
    drop_na(max.ht_t, volume_t) %>% 
    rename("h" = max.ht_t, "v" = volume_t) %>% 
    arrange(v) %>% 
    nlsLM(h ~ A*v^(1/3),
          start = list(A = 0), data = .) %>% 
    coef() %>% 
    as.numeric() -> A
  
  #To see data, store data frame as test and then use plot(test$v, test$h)
  #To see fit line, store the model as fit and then use lines(sort(test$v), fitted(fit), col = "red")
  
  # Function converting volume to height
  vol.to.height <- function(v){
    h <- A*v^(1/3)
    return(h)}
  
  # Vector of n heights across which dispersal kernel will be evaluated
  z.list <- sapply(exp(TM$meshpts), vol.to.height)/100
  
  # List of simulated dispersal distances for each height
  r.list <- as.list(sapply(z.list[z.list >= 0.15], elas=elas, WALD.f.e.h))
  
  # Define modified bessel function for product of s and dispersal distance
  bessel <- function(r, t){
    return(besselI(t*r, 0))}
  
  # Function to evaluate MGF for each height; returns MGFs at each height for a specific value of s
  MGF.s <- function(s){
    
    # No dispersal when z < 0.15; use dirac delta function, with MGF of 1
    mgf.values.a <- rep(1, length(z.list[z.list < 0.15]))
    
    # For all other heights, evaluate bessel at each dispersal distance; find mean across all distances
    sapply(r.list[1:length(r.list)], bessel, t = s) %>% 
      sapply(., mean) -> mgf.values.b
    
    # Return MGF values for all heights
    mgf.values <- c(mgf.values.a, mgf.values.b)
    return(mgf.values)}
  
  # Set up range of s values over which to calculate wavespeeds
  s.seq <- c(0.0001, 0.0005, 0.001, 0.005, seq(0.01, 1, length.out = 46))
  
  # Apply MGF for each value of s
  mgf.over.s <- mapply(MGF.s, s = s.seq)
  
  # Define function to calculate wavespeed for each value of s
  ws.calc <- function(m, s){
    H <- TM$Fmat %*% diag(as.vector(m)) + TM$Pmat
    rho <- Re(eigen(H)$values[1])
    (1/s)*log(rho)}
  
  # Create empty vector on which to add wave speeds
  vec <- c()
  
  # Calculate wavespeed for each s and add it to the vector
  for(i in 1:length(s.seq)){
    val <- ws.calc(m = mgf.over.s[, i], s = s.seq[i])
    vec <- append(vec, val)}
  
  # Return vector of wavespeeds
  return(vec)}





##### Find lambda as function of density ------------------------------------------------------------------

# Function to calculate lambda from density
LambdaD <- function(d.only = FALSE){
  d.values <- seq(min(LATR_full$weighted.dens, na.rm = TRUE), 
                  max(LATR_full$weighted.dens, na.rm = TRUE), length.out = 25)
  if(d.only == TRUE){
    return(d.values)}
  if(d.only == FALSE){
    boot.lambda <- c()
    for(d in 1:length(d.values)){
      boot.lambda[d] <- lambda(TransMatrix(dens = d.values[d], mat.size = 200)$IPMmat)}
    return(boot.lambda)}}


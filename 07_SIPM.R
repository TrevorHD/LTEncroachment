##### IPM functions ---------------------------------------------------------------------------------------

# Growth from size x to y at density d, using best GAM
TM.growth <- function(x, y, d){
  xb = pmin(pmax(x, LATR_size_bounds$min_size), LATR_size_bounds$max_size)
  lpmat <- predict.gam(LATR_grow_best,
                       newdata = data.frame(weighted.dens = d, log_volume_t = xb, unique.transect = "1.FPS"),
                       type = "lpmatrix",
                       exclude = "s(unique.transect)")
  # Linear predictor for mean and log sigma 
  # Need to update so these indices are not hard-coded but for now they work
  grow_mu <- lpmat[, 1:19] %*% coef(LATR_grow_best)[1:19]
  grow_sigma <- exp(lpmat[, 32:50] %*% coef(LATR_grow_best)[32:50])
  return(dnorm(y, mean = grow_mu, sd = grow_sigma))}

# Survival of size x at density d using best GAM
# For nnaturally occuring plants (transplant = FALSE)
TM.survival <- function(x, d){
  xb = pmin(pmax(x, LATR_size_bounds$min_size), LATR_size_bounds$max_size)
  lpmat <- predict.gam(LATR_surv_best,
                       newdata = data.frame(weighted.dens = d, log_volume_t = xb, transplant = FALSE,
                                            unique.transect = "1.FPS"),
                       type = "lpmatrix",
                       exclude = "s(unique.transect)")
  pred <- lpmat[, 1:21] %*% coef(LATR_surv_best)[1:21]
  return(invlogit(pred))}

# Combined growth and survival at density d
TM.growsurv <- function(x, y, d){
  TM.survival(x, d) * TM.growth(x, y, d)}

# Flowering at size x and density d using best GAM
TM.flower <- function(x, d){
  xb = pmin(pmax(x, LATR_size_bounds$min_size), LATR_size_bounds$max_size)
  lpmat <- predict.gam(LATR_flower_best,
                       newdata = data.frame(weighted.dens = d, log_volume_t = xb, unique.transect = "1.FPS"),
                       type = "lpmatrix",
                       exclude = "s(unique.transect)")
  pred <- lpmat[, 1:20] %*% coef(LATR_flower_best)[1:20]
  return(invlogit(pred))}

# Seed production (fruits * seeds/fruit) at size x and density d using best GAM
# Note: we assume 6 seeds per fruit, and best GAM is actually not density dependent
TM.seeds <- function(x, d, seeds.per.fruit = 6){
  xb = pmin(pmax(x, LATR_size_bounds$min_size), LATR_size_bounds$max_size)
  lpmat <- predict.gam(LATR_fruits_best,
                       newdata = data.frame(weighted.dens = d, log_volume_t = xb, unique.transect = "1.FPS"),
                       type = "lpmatrix",
                       exclude = "s(unique.transect)")
  pred <- lpmat[, 1:19] %*% coef(LATR_fruits_best)[1:19]
  return(exp(pred)*seeds.per.fruit)}

# Seed-to-Seedling recruitment probability at density d
TM.recruitment <- function(d){
  lpmat <- predict.gam(LATR_recruit_best,
                       newdata = data.frame(weighted.dens = d, unique.transect = "1.FPS"),
                       type = "lpmatrix",
                       exclude = "s(unique.transect)")
  pred <- lpmat[, 1] %*% coef(LATR_recruit_best)[1]
  return(invlogit(pred[[1]]))}

# Recruit size distribution at size y
TM.recruitsize <- function(y){
  dnorm(x = y, mean = boot.LATR_recruit_size$recruit_mean, sd = boot.LATR_recruit_size$recruit_sd)}

# Combined flowering, fertility, and recruitment
TM.fertrecruit <- function(x, y, d){
  TM.flower(x, d) * TM.seeds(x, d) * TM.recruitment(d) * TM.recruitsize(y)}

# Put it all together; projection matrix is a function of weighted density (dens)
# We need a large lower extension because growth variance (gaussian) is greater for smaller plants
TransMatrix <- function(lower.extension = -8, upper.extension = 2,
                        min.size = LATR_size_bounds$min_size, max.size = LATR_size_bounds$max_size,
                        mat.size = 200, dens){
  
  # Matrix size and size extensions (upper and lower integration limits)
  n <- mat.size
  L <- min.size + lower.extension
  U <- max.size + upper.extension
  
  # Bin size for n bins
  h <- (U - L)/n
  
  # Lower boundaries of bins 
  b <- L + c(0:n)*h
  
  # Bin midpoints
  y <- 0.5*(b[1:n] + b[2:(n + 1)])
  
  # Growth/Survival matrix
  Pmat <- t(outer(y, y, TM.growsurv, d = dens)) * h 
  
  # Fertility/Recruiment matrix
  Fmat <- t(outer(y, y, TM.fertrecruit, d = dens)) * h 
  
  # Put it all together
  IPMmat <- Pmat + Fmat
  
  #and transition matrix
  return(list(IPMmat = IPMmat, Fmat = Fmat, Pmat = Pmat, meshpts = y))}

# Construct transition matrix
TM <- TransMatrix(mat.size = 100, dens = -1.3)





##### Find minimum wave speed -----------------------------------------------------------------------------

# Function to calculate the minimum wavespeed across a range of s
Wavespeed <- function(n){
  
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
  z.list <- sapply(seq(exp(LATR_size_bounds$min_size - 8), exp(LATR_size_bounds$max_size + 2), length.out = n), 
                   vol.to.height)/100
  
  # List of simulated dispersal distances for each height
  r.list <- as.list(sapply(z.list[z.list >= 0.15], WALD.f.e.h))
  
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
  s.seq <- c(0.0001, 0.0005, 0.001, 0.005, seq(0.01, 2, length.out = 96))
  
  # Apply MGF for each value of s
  mgf.over.s <- mapply(MGF.s, s = s.seq)
  
  # Define function to calculate wavespeed for each value of s
  ws.calc <- function(m, s){
    H <- as.vector(m)*TM$Fmat + TM$Pmat
    rho <- Re(eigen(H)$values[1])
    (1/s)*log(rho)}
  
  # Create empty vector on which to add wave speeds
  vec <- c()
  
  # Calculate wavespeed for each s and add it to the vector
  for(i in 1:100){
    val <- ws.calc(m = mgf.over.s[, i], s = s.seq[i])
    vec <- append(vec, val)}
  
  # Return vector of wavespeeds
  return(vec)}


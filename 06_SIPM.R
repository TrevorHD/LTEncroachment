##### Collect all necessary parameters --------------------------------------------------------------------

# Create list of coefficients for demography models
Params <- c()

  # Flowering probability
  Params[1] <- Mod.F.top.cf[1]          # Intercept
  Params[2] <- Mod.F.top.cf[2]          # Volume coefficient
  Params[3] <- Mod.F.top.cf[3]          # Density coefficient
  Params[4] <- Mod.F.top.cf[4]          # Volume and density interaction coefficient
  Params[5] <- Mod.F.top.cf[5]          # Density quadratic coefficient
  Params[6] <- Mod.F.top.cf[6]          # Volume and quadratic density interaction coefficient
  
  # Annual growth
  Params[7] <- Mod.G.top.cf[1]          # Intercept
  Params[8] <- Mod.G.top.cf[2]          # Volume coefficient
  Params[9] <- Mod.G.top.cf[3]          # Density coefficient
  Params[10] <- Mod.G.top.cf[4]         # Volume and density interaction coefficient
  Params[11] <- Mod.G.top.cf[5]         # Density quadratic coefficient
  Params[12] <- Mod.G.top.cf[6]         # Volume and quadratic density interaction coefficient
  
  # Number of reproductive structures
  Params[13] <- Mod.R.top.cf[1]         # Intercept
  Params[14] <- Mod.R.top.cf[2]         # Volume coefficient
  Params[15] <- Mod.R.top.cf[3]         # Density coefficient
  Params[16] <- Mod.R.top.cf[4]         # Volume and density interaction coefficient
  Params[17] <- Mod.R.top.cf[5]         # Density quadratic coefficient
  Params[18] <- Mod.R.top.cf[6]         # Volume and quadratic density interaction coefficient
  
  # Survival (for volume < 7.294)
  Params[19] <- Mod.S.avg.cf[1]         # Intercept
  Params[20] <- Mod.S.avg.cf[2]         # Volume coefficient
  Params[21] <- Mod.S.avg.cf[3]         # Density coefficient
  Params[22] <- Mod.S.avg.cf[4]         # Volume and density interaction coefficient
  Params[23] <- Mod.S.avg.cf[5]         # Density quadratic coefficient
  Params[24] <- Mod.S.avg.cf[6]         # Volume and quadratic density interaction coefficient
  
  # Survival (for volume >= 7.294)
  Params[25] <- 1
  
  # Number of seeds per fruit
  Params[26] <- 5
  
  # Per-seed recruitment probability
  Params[27] <- as.numeric(invlogit(fixef(Mod.P1[[1]])["(Intercept)"]))         # Intercept
  Params[28] <- 0                                                               # Density coefficient
  
  # Minimum and maximum shrub sizes
  Params[29] <- min(boot.CData.s$volume_t, na.rm = TRUE)
  Params[30] <- max(boot.CData.s$volume_t, na.rm = TRUE)
  
  # Minimum and maximum shrub density
  Params[31] <- min(boot.CData.s$d.stand)
  Params[32] <- max(boot.CData.s$d.stand)
  
  # Mean recruit size
  Params[33] <- mean(CData.Recruits$volume_t1, na.rm = TRUE)
  
  # Standard deviation of recruit size
  Params[34] <- sd(CData.Recruits$volume_t1, na.rm = TRUE)
  
  # Growth residual standard deviation
  Params[35] <- sigma(Mod.G[[7]])

# Yes, I know I could have done this much more succinctly
# However, this makes it easier to see what each parameter actually is





##### Construct IPM Kernel --------------------------------------------------------------------------------

# Construct transition matrix
TransMatrix <- function(n, d){
  
  # Growth from size x to y
  xy.Growth <- function(x, y){
    xb <- pmin(pmax(x, Params[29]), Params[30])
    return(dnorm(y, mean = xb + (Params[7] + Params[8]*xb + Params[9]*d + Params[10]*d*xb +
                                 Params[11]*(d^2) + Params[12]*(xb*(d^2))), sd = Params[35]))}
  
  # Survival of size x  
  x.Survival <- function(x){
    xb <- pmin(pmax(x, Params[29]), Params[30])
    val <- ifelse(xb < 7.294, invlogit(Params[19] + Params[20]*xb + Params[21]*d + Params[22]*d*xb +
                                       Params[23]*(d^2) + Params[24]*(xb*(d^2))), Params[25])
    return(val)}
  
  # Growth from size x to y and subsequent survival
  xy.GrowSurv <- function(x, y){
    return(x.Survival(x) * xy.Growth(x, y))}
  
  # Probability of flowering for size x
  x.Flowering <- function(x){
    xb <- pmin(pmax(x, Params[29]), Params[30])
    return(invlogit(Params[1] + Params[2]*xb + Params[3]*d + Params[4]*d*xb +
                    Params[5]*(d^2) + Params[6]*(xb*(d^2))))}
  
  # Number of reproductive structures produced by flowering size x
  x.Reproduction <- function(x){
    xb <- pmin(pmax(x, Params[29]), Params[30])
    return(exp(Params[13] + Params[14]*xb + Params[15]*d + Params[16]*d*xb +
               Params[17]*(d^2) + Params[18]*(xb*(d^2))))}
  
  # Production of y-sized recruits from x-sized adults
  xy.Recruitment <- function(x, y){
    return(x.Flowering(x) * x.Reproduction(x) * Params[26] * 
           Params[27] * dtruncnorm(y, a = Params[29], b = Params[30], 
                                   mean = Params[33], sd = Params[34]))}
  
  # Set upper and lower boundaries slightly beyond data range
  L <- Params[29]*0.9
  U <- Params[30]*1.1
  
  # Create n bins
  b <- L + c(0:n)*(U - L)/n
  
  # Create midpoints in each bin at which functions are evaluated
  y <- 0.5*(b[1:n] + b[2:(n + 1)])
  
  # Evaluate midpoints to construct matrix
  P <- t(outer(y, y, xy.GrowSurv))	
  B <- t(outer(y, y, xy.Recruitment))
  M <- P + B
  M <- (U - L)*M/n
  P <- (U - L)*P/n
  B <- (U - L)*B/n
  return(list(matrix = M, meshpts = y, Pmatrix = P,
              Bmatrix = B))}

# Evaluate the matrix for n subdivisions at the lowest density
TM <- TransMatrix(n = 200, d = -1.3)

# Use Re(eigen(TM$matrix)$values[1]) for geometric growth rate (dominant eigenvalue of TM)





##### Find minimum wave speed -----------------------------------------------------------------------------

# Function to calculate the minimum wavespeed across a range of s
Wavespeed <- function(n){
  
  # Fit equation to convert volume to height for dispersal kernel use
  CData %>%
    select(max.ht_t, volume_t) %>% 
    drop_na(max.ht_t, volume_t) %>% 
    rename("h" = max.ht_t, "v" = volume_t) %>% 
    arrange(v) %>% 
    nlsLM(h ~ A*(exp(v)^(1/3)),
          start = list(A = 0), data = .) %>% 
    coef() %>% 
    as.numeric() -> A
  
  #To see data, store data frame as test and then use plot(test$v, test$h)
  #To see fit line, store the model as fit and then use lines(sort(test$v), fitted(fit), col = "red")
  
  # Function converting volume to height
  vol.to.height <- function(v){
    h <- A*(exp(v)^(1/3))
    return(h)}
  
  # Vector of n heights across which dispersal kernel will be evaluated
  z.list <- sapply(seq(0.9*Params[29], 1.1*Params[30], length.out = n), vol.to.height)/100
  
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
  s.seq <- c(0.0001, 0.0005, 0.001, 0.005, seq(0.01, 2, length.out = 196))
  
  # Apply MGF for each value of s
  mgf.over.s <- mapply(MGF.s, s = s.seq)
  
  # Define function to calculate wavespeed for each value of s
  ws.calc <- function(m, s){
    H <- as.vector(m)*TM$Bmatrix + TM$Pmatrix
    rho <- Re(eigen(H)$values[1])
    (1/s)*log(rho)}
  
  # Create empty vector on which to add wave speeds
  vec <- c()
  
  # Calculate wavespeed for each s and add it to the vector
  for(i in 1:200){
    val <- ws.calc(m = mgf.over.s[, i], s = s.seq[i])
    vec <- append(vec, val)}
  
  # Return vector of wavespeeds
  return(vec)}


##### "Base" WALD PDF -------------------------------------------------------------------------------------

# "Base" WALD PDF using empirical means of observed wind speeds and terminal velocities
# Code adapted from Skarpaas and Shea (2007)
# Using this will assume that height H is the release point for all seeds
WALD.b <- function(H, elas="none",n,seed){
  
  # Add option for height perturbation analysis
  if(elas == "dispersal"){H <- H*(1 + pert)}
  
  # Initialise physical constants
  K <- 0.4      # von Karman constant
  C0 <- 3.125   # Kolmogorov constant
  
  # Initialise other fixed quantities
  Aw <- 1.3     # Ratio of sigmaw to ustar
  h <- 0.15     # Grass cover height
  d <- 0.7*h    # Zero-plane displacement
  z0 <- 0.1*h   # Roughness length
  zm <- 3       # Wind speed measurement height
  
  # Let n be the number of simulation replications
  # Let H be the seed release height
  
  # Find mean wind speed
  Um <- mean(boot.ws.raw)
  
  # Find mean terminal velocity
  f <- mean(boot.tv.raw)
  
  # Calculate ustar, the friction velocity
  ustar <- K*Um*(log((zm - d)/z0))^(-1)

  # Set up integrand for wind speed between vegetation surface and drop height H
  integrand <- function(z){
    (1/K)*(log((z - d)/z0))}
  
  # Integrate to obtain U
  U <- (ustar/H)*integrate(integrand, lower = d + z0, upper = H)$value

  # Calculate instability parameter
  sigma <- 2*(Aw^2)*sqrt((K*(H - d)*ustar)/(C0*U))

  # Calculate scale parameter lambda
  lambda <- (H/sigma)^2
  
  # Calculate location parameter nu
  nu <- H*U/f
  
  set.seed(seed)
  return(rinvGauss(n, nu = nu, lambda = lambda))}



WALD.b.tom <- function(H, elas="none"){

  # Initialise physical constants
  K <- 0.4      # von Karman constant
  C0 <- 3.125   # Kolmogorov constant
  
  # Initialise other fixed quantities
  Aw <- 1.3     # Ratio of sigmaw to ustar
  h <- 0.15     # Grass cover height
  d <- 0.7*h    # Zero-plane displacement
  z0 <- 0.1*h   # Roughness length
  zm <- 3       # Wind speed measurement height
  
  # Let n be the number of simulation replications
  # Let H be the seed release height
  
  # Find mean wind speed
  Um <- mean(boot.ws.raw)
  
  # Find mean terminal velocity
  f <- mean(boot.tv.raw)
  
  # Calculate ustar, the friction velocity
  ustar <- K*Um*(log((zm - d)/z0))^(-1)
  
  # Set up integrand for wind speed between vegetation surface and drop height H
  integrand <- function(z){
    (1/K)*(log((z - d)/z0))}
  
  # Integrate to obtain U
  U <- (ustar/H)*integrate(integrand, lower = d + z0, upper = H)$value
  
  # Calculate instability parameter
  sigma <- 2*(Aw^2)*sqrt((K*(H - d)*ustar)/(C0*U))
  
  # Calculate scale parameter lambda
  lambda <- (H/sigma)^2
  if(elas == "dispersal.scale"){lambda <- lambda*(1 + pert)}
  
  # Calculate location parameter nu
  nu <- H*U/f
  if(elas == "dispersal.location"){nu <- nu*(1 + pert)}
  
  return(list(lambda=lambda,nu=nu))
}

##### "Full" WALD PDF with empirical distributions --------------------------------------------------------
# "Full" WALD PDF, including distributions of wind speeds and terminal velocities
# Code adapted from Skarpaas and Shea (2007)
WALD.f.e <- function(n, H, elas, seed = NULL){

  # Initialise physical constants
  K <- 0.4      # von Karman constant
  C0 <- 3.125   # Kolmogorov constant
  
  # Initialise other fixed quantities
  Aw <- 1.3     # Ratio of sigmaw to ustar
  h <- 0.15     # Grass cover height
  d <- 0.7*h    # Zero-plane displacement
  z0 <- 0.1*h   # Roughness length
  zm <- 3       # Wind speed measurement height
  
  # Let n be the number of simulation replications
  # Let H be the seed release height
  
  # Simulate wind speeds from empirical distribution of wind speeds
  set.seed(seed)
  Um <- rnorm(n, sample(boot.ws.raw, size = n, replace = TRUE), boot.ws.PDF$bw)
  
  # Simulate terminal velocities from empirical distribution of terminal velocities
  set.seed(seed)
  f <- rnorm(n, sample(boot.tv.raw, size = n, replace = TRUE), boot.tv.PDF$bw)
  
  # Calculate ustar, the friction velocity
  ustar <- K*Um*(log((zm - d)/z0))^(-1)
  
  # Set up integrand for wind speed between vegetation surface and drop height H
  integrand <- function(z){
    (1/K)*(log((z - d)/z0))}
  
  # Integrate to obtain U
  U <- (ustar/H)*integrate(integrand, lower = d + z0, upper = H)$value
  
  # Calculate instability parameter
  sigma <- 2*(Aw^2)*sqrt((K*(H - d)*ustar)/(C0*U))
  
  # Calculate scale parameter lambda
  lambda <- (H/sigma)^2
  if(elas == "dispersal.scale"){lambda <- lambda*(1 + pert)}
  
  # Calculate location parameter nu
  nu <- H*U/f
  if(elas == "dispersal.location"){nu <- nu*(1 + pert)}
  
  # Generate inverse Gaussian distribution
  set.seed(seed)
  return(rinvGauss(n, nu = nu, lambda = lambda))
  }





##### "Base" WALD PDF release across entire height --------------------------------------------------------

WALD.b.h <- function(H, elas, seed = NULL, reps, heights){
  
  # Use 10000 replicates for each height #TM: could we do less?
  n <- reps
  
  # Create "continuous" sequence of release heights
  h.range <- seq(0.15, H, length.out = heights) ##TM: could we make this chunkier
  
  # Simulate seed release events for each height; returns n*length.out dispersal events
  return(na.omit(as.vector(sapply(h.range, WALD.b, n = n, elas = elas, seed = seed))))}


##### "Full" empirical WALD PDF release across entire height ----------------------------------------------

WALD.f.e.h <- function(H, elas, seed = NULL, reps, heights){
  
  # Use 10000 replicates for each height #TM: could we do less?
  n <- reps
  
  # Create "continuous" sequence of release heights
  h.range <- seq(0.15, H, length.out = heights) ##TM: could we make this chunkier
  
  # Simulate seed release events for each height; returns n*length.out dispersal events
  return(na.omit(as.vector(sapply(h.range, WALD.f.e, n = n, elas = elas, seed = seed))))}



##### Tom's version of full WALD
WALD.f.e.h.tom <- function(n, H, elas, h=0.15, seed = NULL){
  
  # Add option for height perturbation analysis
  if(elas == "dispersal"){H <- H*(1 + pert)}
  
  # Initialise physical constants
  K <- 0.4      # von Karman constant
  C0 <- 3.125   # Kolmogorov constant
  
  # Initialise other fixed quantities
  Aw <- 1.3     # Ratio of sigmaw to ustar
  h <- 0.15     # Grass cover height
  d <- 0.7*h    # Zero-plane displacement
  z0 <- 0.1*h   # Roughness length
  zm <- 3       # Wind speed measurement height
  
  # Simulate wind speeds from empirical distribution of wind speeds
  set.seed(seed)
  Um <- sample(boot.ws.raw, size = n, replace = TRUE)
  
  # Simulate terminal velocities from empirical distribution of terminal velocities
  set.seed(seed)
  f <- sample(boot.tv.raw, size = n, replace = TRUE)
  
  # Calculate ustar, the friction velocity
  ustar <- K*Um*(log((zm - d)/z0))^(-1)
  
  # Set up integrand for wind speed between vegetation surface and drop height H
  integrand <- function(z){
    (1/K)*(log((z - d)/z0))}
  
  # Integrate to obtain U
  U <- (ustar/H)*integrate(integrand, lower = d + z0, upper = H)$value
  
  # Calculate instability parameter
  sigma <- 2*(Aw^2)*sqrt((K*(H - d)*ustar)/(C0*U))
  
  # Calculate scale parameter lambda
  lambda <- (H/sigma)^2
  if(elas == "dispersal.scale"){lambda <- lambda*(1 + pert)}
  
  # Calculate location parameter nu
  nu <- H*U/f
  if(elas == "dispersal.location"){nu <- nu*(1 + pert)}
  
  # Generate inverse Gaussian distribution
  set.seed(seed)
  return(rinvGauss(n, nu = nu, lambda = lambda))
}

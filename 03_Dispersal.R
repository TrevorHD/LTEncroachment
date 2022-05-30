##### "Base" WALD PDF -------------------------------------------------------------------------------------

# "Base" WALD PDF, including distributions of wind speeds and terminal velocities
# Code adapted from Skarpaas and Shea (2007)
WALD.b <- function(n, H){
  
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
  
  # Generate inverse Gaussian distribution
  return(rinvGauss(n, nu = nu*delta, lambda = lambda))}





##### "Full" WALD PDF with empirical distributions --------------------------------------------------------

# "Full" WALD PDF, including distributions of wind speeds and terminal velocities
# Code adapted from Skarpaas and Shea (2007)
WALD.f.e <- function(n, H){
  
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
  Um <- rnorm(n, sample(boot.ws.raw, size = n, replace = TRUE), boot.ws.PDF$bw)
  
  # Simulate terminal velocities from empirical distribution of terminal velocities
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
  
  # Calculate location parameter nu
  nu <- H*U/f

  # Generate inverse Gaussian distribution
  return(rinvGauss(n, nu = nu, lambda = lambda))}





##### "Full" empirical WALD PDF release across entire height ----------------------------------------------

WALD.f.e.h <- function(H){
  
  # Use 10000 replicates for each height
  n <- 10000
  
  # Create "continuous" sequence of release heights
  h.range <- seq(0.15, H, length.out = 50)
  
  # Simulate seed release events for each height
  return(na.omit(as.vector(sapply(h.range, WALD.f.e, n = n))))}


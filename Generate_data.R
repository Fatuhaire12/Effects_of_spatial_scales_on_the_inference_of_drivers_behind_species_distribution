# Simulation parameters
Numsim <- 1001  # Number of simulations

# Landscape parameters
maxstep <- 7
np <- 2^(maxstep - 1)  # Size of the fractal landscape
xl <- 10  # Length of domain in x and y directions
dx <- 2 * xl / np  # Grid spacing
x <- seq(-xl, xl - dx, dx)
y <- x
[X, Y] <- meshgrid(x, y)

# Time parameters
ngens <- 200  # Number of generations
dt <- 1.0  # Time step length

# Dispersal kernel
df <- 1.0  # Diffusion coefficient
hker <- dx^2 * exp(-(X^2 + Y^2) / (2 * df * dt)) / (2 * pi * dt * df)
Fhker <- fft2(hker)

g <- 1
while (g < Numsim) {
  # Fractal k
  h <- runif(1)
  k0 <- fractal(h, maxstep)  # Fractal one
  k <- transformation(k0, 0.25, 1)  # Case one
  k <- k[1:(end - 1), 1:(end - 1)]
  
  # Fractal r
  h <- runif(1)
  r0 <- fractal(h, maxstep)  # Fractal two
  r <- transformation(r0, -0.5, 1)
  r <- r[1:(end - 1), 1:(end - 1)]
  
  # Initialize population
  p <- matrix(0, nrow = np, ncol = np)
  origin <- floor(np / 2)
  p[origin - 3:origin + 3, origin - 3:origin + 3] <- 1
  
  for (j in 1:ngens) {
    hn <- p * exp(r * (1 - p / k))  # Ricker growth model
    fhn <- fft2(hn)  # FFT of the species
    
    # Shift to center the probability functions
    p <- Re(fftshift(ifft2(Fhker * fhn)))
    maxpop <- max(p, na.rm = TRUE)
    
    if (is.nan(maxpop)) {
      break
    }
  }
  
  if (j == ngens) {
    # Write CSV datasets here
    fname <- sprintf('valuesofR%d.csv', g)
    write.csv(r, file = fname)
    
    fname <- sprintf('valuesofK%d.csv', g)
    write.csv(k, file = fname)
    
    fname <- sprintf('valuesofN%d.csv', g)
    write.csv(p, file = fname)
    
    g <- g + 1
  }
}

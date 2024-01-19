# Load required libraries
library("gtools")
library("hier.part")
library("ape")
library("fractaldim")
library("lattice")
library("reshape2")
library("MESS")

# Function to calculate scaling effect
myfun <- function(popden, valuesofR, valuesofK, tmpdist) {
  values <- c()
  
  # Hierarchical partitioning
  independent <- data.frame(as.vector(valuesofR), as.vector(valuesofK))
  hierarchy <- hier.part(as.vector(popden), independent, family = "gaussian", gof = "Rsqu", barplot = FALSE)
  values <- c(values, hierarchy$gfs[c(4, 2, 3)])
  
  # Morans' I
  values <- c(values, c(Moran.I(as.vector(valuesofR), tmpdist)$observed))
  values <- c(values, c(Moran.I(as.vector(valuesofK), tmpdist)$observed))
  values <- c(values, c(Moran.I(as.vector(popden), tmpdist)$observed))
  
  # Fractal dimension
  values <- c(values, fractaldim::fd.estim.boxcount(valuesofR)$fd)
  values <- c(values, fractaldim::fd.estim.boxcount(valuesofK)$fd)
  values <- c(values, fractaldim::fd.estim.boxcount(popden)$fd)
  
  return(values)
}

# Function to scale data
scaledata <- function(mm, ss) {
  dd <- dim(mm)
  dx <- dd[1]
  dy <- dd[2]
  dscale = floor(dx / ss)
  tmpmat <- matrix(0, dscale, dscale)
  for (i in 0:(dscale - 1)) {
    for (j in 0:(dscale - 1)) {
      tmpmat[i + 1, j + 1] = mean(mm[(i * ss + 1):((i + 1) * ss), (j * ss + 1):((j + 1) * ss)])
    }
  }
  return(tmpmat)
}

# Function to get scaled data
getdataScales <- function(i) {
  scaleddata <- c()
  valuesofR <- t(read.csv(paste0("valuesofR", i, ".csv"), header = FALSE, sep = ","))
  valuesofK <- t(read.csv(paste0("valuesofK", i, ".csv"), header = FALSE, sep = ","))
  popden <- t(read.csv(paste0("valuesofN", i, ".csv"), header = FALSE, sep = ","))
  
  for (ss in c(2, 4, 8, 16)) {
    mergedR = scaledata(valuesofR, ss)
    mergedK = scaledata(valuesofK, ss)
    mergedP = scaledata(popden, ss)
    np = floor(64 / ss)
    tmpx = ss * rep(1:np, np)
    dim(tmpx) <- c(np, np)
    tmpdist = 1 / as.matrix(dist(cbind(as.vector(tmpx), as.vector(t(tmpx)))))
    
    diag(tmpdist) <- 0
    scaleddata = rbind(scaleddata, myfun(mergedP, mergedR, mergedK, tmpdist))
  }
  
  np = 64
  tmpx = rep(1:np, np)
  dim(tmpx) <- c(np, np)
  tmpdist = 1 / as.matrix(dist(cbind(as.vector(tmpx), as.vector(t(tmpx)))))
  
  diag(tmpdist) <- 0
  scaleddata = rbind(scaleddata, myfun(popden, valuesofR, valuesofK, tmpdist))
  
  return(scaleddata)
}

# Initialize result variables
results.rawdata <- c()
results.scale2 <- c()
results.scale4 <- c()
results.scale8 <- c()
results.scale16 <- c()

numsim <- 100

# Loop through simulations
for (i in 1:numsim) {
  results <- getdataScales(i)
  results.rawdata <- rbind(results.rawdata, results[1, ])
  results.scale2 <- rbind(results.scale2, results[2, ])
  results.scale4 <- rbind(results.scale4, results[3, ])
  results.scale8 <- rbind(results.scale8, results[4, ])
  results.scale16 <- rbind(results.scale16, results[5, ])
}

# Name the columns
colnames(results.rawdata) <- c('R2RK', 'R2R', 'R2K', 'MoransR', 'MoransK', 'MoransN', 'fdR', 'fdK', 'fdN')
colnames(results.scale2) <- colnames(results.scale4) <- colnames(results.scale8) <- colnames(results.scale16) <- c('R2RK', 'R2R', 'R2K', 'MoransR', 'MoransK', 'MoransN', 'fdR', 'fdK', 'fdN')

# Data frames for contributions
aR <- data.frame(results.rawdata[, 1], results.scale2[, 1], results.scale4[, 1], results.scale8[, 1], results.scale16[, 1])
aK <- data.frame(results.rawdata[, 2], results.scale2[, 2], results.scale4[, 2], results.scale8[, 2], results.scale16[, 2])
aN <- data.frame(results.rawdata[, 3], results.scale2[, 3], results.scale4[, 3], results.scale8[, 3], results.scale16[, 3])

# Plot Moran's I
aR1 <- aR
names(aR1) <- as.character(c("All data", 2, 4, 8, 16))
head(aR1, 2)
head(melt(aR1))
tail(melt(aR1))

aK1 <- aK
names(aK1) <- as.character(c("All data", 2, 4, 8, 16))

aN1 <- aN
names(aN1) <- as.character(c("All data", 2, 4, 8, 16))

a <- melt(aR1)
b <- melt(aK1)
c <- melt(aN1)
dim(a)

# Label Parameters
PARAMETERS <- rep(c("R2RK", "RR", "RK"), each = 600)
head(PARAMETERS)

# Combine data frames
e <- rbind(a, c, b)
dim(e)

# Create a data frame
all.data <- data.frame(e, PARAMETERS)
head(all.data)

# Plot
ggplot(data = all.data, aes(x = variable, y = value)) +
  geom_boxplot(aes(fill = PARAMETERS), width = 1.0, main = "Morans") +
  theme_bw() + theme(axis.text = element_text(size = 25, face = "bold"),
                     axis.title = element_text(size = 25, face = "bold")) +
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

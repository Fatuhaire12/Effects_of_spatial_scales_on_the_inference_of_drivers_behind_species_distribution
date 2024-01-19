#
# Install and load required packages
library("gtools")
library("hier.part")
library("abind")
library("ape")
library("fractaldim")
library("lattice")
library("reshape2")
library("reshape")
library("ggplot2")

# Function to calculate HP, Moran's I, and FD of sampled data
myfun <- function(popden, valuesofR, valuesofK, tmpdist) {
  values <- c()
  
  # HP
  independent <- data.frame(valuesofR, valuesofK)
  hierarchy <- hier.part(popden, independent, family="gaussian", gof="Rsqu", barplot=FALSE)
  values <- c(values, hierarchy$gfs[c(4, 2, 3)])
  
  # Moran's I
  values <- c(values, c(Moran.I(valuesofR, tmpdist)$observed))
  values <- c(values, c(Moran.I(valuesofK, tmpdist)$observed))
  values <- c(values, c(Moran.I(popden, tmpdist)$observed))
  
  # FD
  values <- c(values, fractaldim::fd.estim.boxcount(valuesofR)$fd)
  values <- c(values, fractaldim::fd.estim.boxcount(valuesofK)$fd)
  values <- c(values, fractaldim::fd.estim.boxcount(popden)$fd)
  
  return(values)
}

# Function to get sampled data
getdataSampling <- function(i) {
  sampleddata <- c()
  valuesofR <- t(read.csv(paste0("valuesofR", i, ".csv"), header=FALSE, sep=","))
  valuesofK <- t(read.csv(paste0("valuesofK", i, ".csv"), header=FALSE, sep=","))
  popden <- t(read.csv(paste0("valuesofN", i, ".csv"), header=FALSE, sep=","))
  
  for (nsample in c(5, 10, 50, 100, 500, 1000, 4000)) {
    locs <- sample(0:4095, nsample, replace=FALSE)
    locx <- 1 + (locs %% 64)
    locy <- 1 + (locs %% 64)
    
    sampledR <- matrix(nrow=nsample, ncol=1)
    sampledK <- matrix(nrow=nsample, ncol=1)
    sampledN <- matrix(nrow=nsample, ncol=1)
    
    for (t in 1:nsample) {
      sampledR[t] = valuesofR[locx[t], locy[t]]
      sampledK[t] = valuesofK[locx[t], locy[t]]
      sampledN[t] = popden[locx[t], locy[t]]
    }
    
    tmpdist = 1/as.matrix(dist(cbind(locx, locy)))
    diag(tmpdist) <- 0
    tmpdist[is.infinite(tmpdist)] <- 0
    
    sampleddata = abind(sampleddata, myfun(as.vector(sampledN), as.vector(sampledR), as.vector(sampledK), tmpdist))
  }
  
  return(sampleddata)
}

# Initiate an empty list for the sampled data
results.sample5 = c()
results.sample10 = c()
results.sample50 = c()
results.sample100 = c()
results.sample500 = c()
results.sample1000 = c()
results.sample4000 = c()

numsim = 100

for (i in 1:numsim) {
  results = getdataSampling(i)
  results.sample5 = rbind(results.sample5, results[1,])
  results.sample10 = rbind(results.sample10, results[2,])
  results.sample50 = rbind(results.sample50, results[3,])
  results.sample100 = rbind(results.sample100, results[4,])
  results.sample500 = rbind(results.sample500, results[5,])
  results.sample1000 = rbind(results.sample1000, results[6,])
  results.sample4000 = rbind(results.sample4000, results[7,])
}

# Name the columns in the results of the sampled data
colnames(results.sample5) <- c('R2_RK', 'R2_R', 'R2_K', 'MoransR', 'MoransK', 'MoransN', 'fdR', 'fdK', 'fdN')
colnames(results.sample10) <- colnames(results.sample5)
colnames(results.sample50) <- colnames(results.sample5)
colnames(results.sample100) <- colnames(results.sample5)
colnames(results.sample500) <- colnames(results.sample5)
colnames(results.sample1000) <- colnames(results.sample5)
colnames(results.sample4000) <- colnames(results.sample5)

# Plot the contribution of each variable
aK <- data.frame(results.sample5[, 1], results.sample10[, 1], results.sample50[, 1], 
                 results.sample100[, 1], results.sample1000[, 1], results.sample4000[, 1])
aN <- data.frame(results.sample5[, 2], results.sample10[, 2], results.sample50[, 2], 
                 results.sample100[, 2], results.sample1000[, 2], results.sample4000[, 2])
aR <- data.frame(results.sample5[, 3], results.sample10[, 3], results.sample50[, 3], 
                 results.sample100[, 3], results.sample1000[, 3], results.sample4000[, 3])

# Plot Moran's I
aR <- data.frame(results.sample5[, 4], results.sample10[, 4], results.sample50[, 4], 
                 results.sample100[, 4], results.sample1000[, 4], results.sample4000[, 4])
aK <- data.frame(results.sample5[, 5], results.sample10[, 5], results.sample50[, 5], 
                 results.sample100[, 5], results.sample1000[, 5], results.sample4000[, 5])
aN <- data.frame(results.sample5[, 6], results.sample10[, 6], results.sample50[, 6], 
                 results.sample100[, 6], results.sample1000[, 6], results.sample4000[, 6])

# Plot the Fractal dimension
aR <- data.frame(results.sample5[, 7], results.sample10[, 7], results.sample50[, 7], 
                 results.sample100[, 7], results.sample1000[, 7], results.sample4000[, 7])
aK <- data.frame(results.sample5[, 8], results.sample10[, 8], results.sample50[, 8], 
                 results.sample100[, 8], results.sample1000[, 8], results.sample4000[, 8])
aN <- data.frame(results.sample5[, 9], results.sample10[, 9], results.sample50[, 9], 
                 results.sample100[, 9], results.sample1000[, 9], results.sample4000[, 9])

# Plot Moran's I of sampled data
a <- melt(aR)
b <- c(5, 10, 50, 100, 1000)
aR1 <- aR
names(aR1) <- as.character(c(b, "All data"))
head(aR1, 2)
head(melt(aR1))
tail(melt(aR1))

aK1 <- aK
names(aK1) <- as.character(c(b, "All data"))

aN1 <- aN
names(aN1) <- as.character(c(b, "All data"))

a <- melt(aR1)
b <- melt(aK1)
c <- melt(aN1)
dim(a)

# Set parameters
PARAMETERS <- rep(c("R2RK", "RR", "RK"), each=600)
head(PARAMETERS)

# Combine data frames
e <- rbind(a, c, b)
dim(e)

# Create a data frame
all.data <- data.frame(e, PARAMETERS)
head(all.data)

# Plot
ggplot(data = all.data, aes(x=variable, y=value)) +
  geom_boxplot(aes(fill=PARAMETERS), width=1.0, main="Morans") +
  theme_bw() + theme(axis.text=element_text(size=25, face="bold"), axis.title=element_text(size=25, face="bold")) +
  theme(axis.line=element_line(colour="black", size=1, linetype="solid"),panel.border = element_rect(colour="black", fill=NA, size=1)) +
  geom_line() + labs(x="Sample Size", y="Contribution") +
  theme(plot.title = element_text(family="Trebuchet MS", color="black", face="bold", size=32, hjust=0))

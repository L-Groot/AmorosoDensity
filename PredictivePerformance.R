##########################################
### Investigate predictive performance ###
##########################################

#------------------
# Source Functions
#------------------

# -> Source functions
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoDensity/refs/",
              "heads/main/estimate_amoroso_np.R"))

library(AmoRosoDistrib)


#-----------------------
# Predictive Performance
#-----------------------

# Generate data from Amoroso with certain parameters
# dat <- rgg4()
dat <- rgg4(1000,a=-1,l=2,c=1,mu=5)
#hist(dat, xlim = c(0,10), breaks = 20)
#range(dat)
#plot(density(dat))

# Then fit NP methods, R density() and Amoroso to that data
# estimate_amoroso_np()

# --> Amoroso should win (have highest likelihood)
amo <- estimate_amoroso(dat, criterion = "ML")
res <- estimate_amoroso_np(dat, hist = TRUE, amorosocrit = "ML")


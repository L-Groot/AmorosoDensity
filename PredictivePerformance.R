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

### Generate data from Amoroso

# dat <- rgg4()
dat <- rgg4(1000,a=-1,l=2,c=1,mu=5)
#hist(dat, xlim = c(0,10), breaks = 20)
#range(dat)
#plot(density(dat))


### Split data in training and test set



### Fit Amoroso and NP methods
res <- estimate_amoroso_np(dat, hist = TRUE, amorosocrit = "ML")


### Make continuous functions from NP fits

par(mfrow = c(2,4))

# -> Bernstein 1
f_bern1 <- splinefun(res$bern1$x,res$bern1$y, method = "monoH.FC")
plot(f_bern1, xlim = c(min(res$bern1$x),max(res$bern1$x)))
points(res$bern1$x,res$bern1$y, col = "magenta", cex = 0.2)

# -> Bernstein 2
f_bern2 <- splinefun(res$bern2$x,res$bern2$y, method = "monoH.FC")
plot(f_bern2, xlim = c(min(res$bern2$x),max(res$bern2$x)))
#plot(f_bern2, xlim = c(min(res$bern2$x),6))
points(res$bern2$x,res$bern2$y, col = "springgreen3", cex = 0.2)

# adj. KDE (unimodal)
f_scKDE_uni <- splinefun(res$scKDE_uni$x,res$scKDE_uni$y, method = "monoH.FC")
plot(f_scKDE_uni, xlim = c(min(res$scKDE_uni$x),max(res$scKDE_uni$x)))
points(res$scKDE_uni$x,res$scKDE_uni$y, col = "cornflowerblue", cex = 0.2)

# adj. KDE (2 inflections)
f_scKDE_2inf <- splinefun(res$scKDE_2inf$x,res$scKDE_2inf$y, method = "monoH.FC")
plot(f_scKDE_2inf, xlim = c(min(res$scKDE_2inf$x),max(res$scKDE_2inf$x)))
points(res$scKDE_2inf$x,res$scKDE_2inf$y, col = "darkorange2", cex = 0.2)

# adj. KDE (2 inflections +)
f_scKDE_2infplus <- splinefun(res$scKDE_2infplus$x,res$scKDE_2infplus$y, method = "monoH.FC")
plot(f_scKDE_2infplus, xlim = c(min(res$scKDE_2infplus$x),max(res$scKDE_2infplus$x)))
points(res$scKDE_2infplus$x,res$scKDE_2infplus$y, col = "purple3", cex = 0.2)

# R density
f_rdens <- splinefun(res$rdens$x,res$rdens$y, method = "monoH.FC")
plot(f_rdens, xlim = c(min(res$rdens$x),max(res$rdens$x)))
points(res$rdens$x,res$rdens$y, col = "red", cex = 0.2)

### Also plot Amoroso
plot(dgg4(res$amo$x,res$amo$pars[1],res$amo$pars[2],res$amo$pars[3],res$amo$pars[4]),
     xlim = c(min(res$amo$x),max(res$amo$x)), ylab = "Amoroso")
points(res$amo$x,res$amo$y, col = "gold3", cex = 0.2)



#-------------------------------------------------------------------------------

# Load functions and packages
# -> for estimating Amoroso
source("https://raw.githubusercontent.com/L-Groot/AmorosoDensity/refs/heads/main/estimate_amoroso.R")
# -> for estimating Bernstein
source("https://raw.githubusercontent.com/L-Groot/AmorosoDensity/refs/heads/main/estimate_bernstein.R")
# -> for estimating adjusted KDE
if (!requireNamespace("scdensity", quietly = TRUE)) {
  install.packages("scdensity")}
library(scdensity)
# -> for Amoroso
# -> for estimating adjusted KDE
if (!requireNamespace("AmoRosoDistrib", quietly = TRUE)) {
  install.packages("AmoRosoDistrib")}
library(AmoRosoDistrib)

#-------------------------------------------------------------------------------
#dat <- palmerpenguins::penguins$bill_length_mm
#dat <- rnrom(100)
#-------------------------------------------------------------------------------

estimate_amoroso_np <- function(dat, main = NULL, minimal = TRUE, plot = TRUE) {
  
  # Remove NA
  num_nas_to_remove <- length(dat) - length(na.omit(dat))
  dat <- as.vector(na.omit(dat))
  
  # Print how many NAs were removed
  if(num_nas_to_remove > 0) {
    message(c("WARNING: ", as.character(num_nas_to_remove), " NAs were removed from the data.\n"))
    cat("--------------------------------------------------------------------\n")
  }
  
  # Get n
  n <- length(dat)
  
  ## Estimate density with different methods ##
  # Amoroso
  amo <- estimate_amoroso(dat, plot=0, criterion="maxL")
  # Bernstein
  bern1 <- estimate_bernstein(dat, bound_type = "sd")
  bern2 <- estimate_bernstein(dat, bound_type = "Carv")
  # Adjusted KDE
  scKDE_2infplus <- scdensity(dat, constraint = "twoInflections+")
  scKDE_2inf <- scdensity(dat, constraint = "twoInflections")
  scKDE_uni <- scdensity(dat, constraint = "unimodal")
  
  # Extract Amoroso parameters
  amo_xx <- amo$x
  amo_ML <- amo$max_L_model # get ML Amoroso
  amo_ML_pars <- as.vector(unlist(amo_ML[,3:6]))
  amo_ML_name <- paste0(amo_ML$method, " (", amo_ML$space, ")")
  
  
  # Get Amoroso density values
  
  dAmoroso(amo_xx,19.2342455,0.3355688,5.0602110,31.9464149)
  
  
  amo_yy <- dAmoroso(amo_xx, amo_ML_pars[1], amo_ML_pars[2], amo_ML_pars[3], amo_ML_pars[4])
  amo_yy <- dgg4(amo_xx, amo_ML_pars[1], amo_ML_pars[2], amo_ML_pars[3], amo_ML_pars[4])
  # Replace NAs by zeroes
  amo_yy[is.na(amo_yy)] <- 0
  # Remove NA
  # non_na_indices <- !is.na(amo_yy)
  # amo_xx <- amo_xx[non_na_indices]
  # amo_yy <- amo_yy[non_na_indices]
  
  # Define the range of x values
  #x_range <- seq(-3, 3, length.out = 1000)
  buffer <- (max(dat)-min(dat))/10
  x_range <- seq(min(dat)-buffer,max(dat)+buffer,length.out = 1000)
  
  # Define max of y values
  y_max <- max(density(dat)$y,amo_yy,bern1$y,bern2$y,
               scKDE_2infplus$y,scKDE_2inf$y,scKDE_uni$y) + 0.05
  
  if (plot == TRUE) {
    if (minimal == FALSE) {
      
      # Initialize 2x3 grid
      par(mfrow=c(2,3), oma = c(0, 0, 5, 0), cex.axis = 0.9, font.lab = 2, font.axis = 1)
      
      # Plot R density() KDE
      plot(x_range, dnorm(x_range), ylim = c(0.0,y_max), type = "l",
           lwd = 1, lty = 2, col = "grey30",
           main = "R density() KDE",
           axes = F, xlab = "x", ylab = "Density")
      axis(1)
      axis(2, las=2)
      lines(density(dat), col = 'mediumorchid2', lwd = 2)
      rug(dat, col = "blue", lwd = 1)
      
      # Plot Bernstein
      plot(x_range, dnorm(x_range), ylim = c(0.0,y_max), type = "l",
           lwd = 1, lty = 2, col = "grey30",
           main = "Bernstein Polynomial",
           axes = F, xlab = "x", ylab = "Density")
      axis(1)
      axis(2, las=2)
      lines(bern1$x, bern1$y, col = "mediumorchid2", lwd = 2)
      lines(bern2$x, bern2$y, col = "chartreuse4", lwd = 2)
      rug(dat, col = "blue", lwd = 1)
      legend("topright", legend = c("bound.type = 'sd'","bound.type = 'Carv'"),
             col = c("mediumorchid2","chartreuse4"), lty = 1, lwd = 2, cex = 0.8, bty = "n")
      
      # Plot Amoroso (lowest BIC)
      plot(x_range, dnorm(x_range), ylim = c(0.0,y_max), type = "l",
           lwd = 1, lty = 2, col = "grey30",
           main = amoroso_name,
           axes = F, xlab = "x", ylab = "Density")
      axis(1)
      axis(2, las=2)
      lines(amo_xx, amo_yy, col = "mediumorchid2", lwd = 2)
      rug(dat, col = "blue", lwd = 1)
      
      # Plot adj. KDE, constraint = "unimodal"
      plot(x_range, dnorm(x_range), ylim = c(0.0,y_max), type = "l",
           lwd = 1, lty = 2, col = "grey20",
           main = "Adj. KDE — 'unimodal'",
           axes = F, xlab = "x", ylab = "Density")
      axis(1)
      axis(2, las=2)
      lines(scKDE_uni$x, scKDE_uni$y, col = "mediumorchid2", lwd = 2)
      rug(dat, col = "blue", lwd = 1)
      
      # Plot adj. KDE, constraint = "twoInflections"
      plot(x_range, dnorm(x_range), ylim = c(0.0,y_max), type = "l",
           lwd = 1, lty = 2, col = "grey30",
           main = "Adj. KDE — 'twoInflections'",
           axes = F, xlab = "x", ylab = "Density")
      axis(1)
      axis(2, las=2)
      lines(scKDE_2inf$x, scKDE_2inf$y, col = "mediumorchid2", lwd = 2)
      rug(dat, col = "blue", lwd = 1)
      
      # Plot adj. KDE, constraint = "twoInflections+"
      plot(x_range, dnorm(x_range), ylim = c(0.0,y_max), type = "l",
           lwd = 1, lty = 2, col = "grey30",
           main = "Adj. KDE — 'twoInflections+'",
           axes = F, xlab = "x", ylab = "Density")
      axis(1)
      axis(2, las=2)
      lines(scKDE_2infplus$x, scKDE_2infplus$y, col = "mediumorchid2", lwd = 2)
      rug(dat, col = "blue", lwd = 1)
      
      # Add an overall title
      if (is.null(main)) {
        big_title <- paste0("rnorm(", as.character(n), ")")
      } else {
        big_title <- main
      }
      
      #mtext(big_title, outer = TRUE, cex = 1.5)
      mtext(big_title, outer = TRUE, cex = 1.5, line = 2, font = 2)
    }
    
    else {
      
      # Initialize 2x3 grid
      par(mfrow=c(2,3), oma = c(1, 6, 1, 1), mar = c(5,5,5,5), cex.axis = 1.4,
          font.lab = 2, font.axis = 1, family = "Times New Roman")
      
      # Plot R density() KDE
      plot(x_range, dnorm(x_range), xlim = c(-4, 4), ylim = c(0.0,y_max), type = "l",
           lwd = 1, lty = 2, col = "grey30",
           main = "",
           axes = F, xlab="", ylab ="")
      mtext("R density()", side=3, font=2, cex=1.5, line=1)
      axis(1, at = c(-4, 0, 4), labels = c(-4, 0, 4))
      lines(density(dat), col = 'blueviolet', lwd = 2)
      rug(dat, col = "dodgerblue3", lwd = 1)
      
      # Plot Bernstein
      plot(x_range, dnorm(x_range), xlim = c(-4, 4), ylim = c(0.0,y_max), type = "l",
           lwd = 1, lty = 2, col = "grey30",
           main = "",
           axes = F, xlab="", ylab ="")
      mtext("Bernstein", side=3, font=2, cex=1.5, line=1)
      axis(1, at = c(-4, 0, 4), labels = c(-4, 0, 4))
      lines(bern1$x, bern1$y, col = "blueviolet", lwd = 2)
      lines(bern2$x, bern2$y, col = "deeppink2", lwd = 2)
      rug(dat, col = "dodgerblue3", lwd = 1)
      # legend("topright", legend = c("bound.type = 'sd'","bound.type = 'Carv'"),
      #        col = c("blueviolet","deeppink2"), lty = 1, lwd = 2, cex = 0.8, bty = "n")
      
      # Plot Amoroso (lowest BIC)
      plot(x_range, dnorm(x_range), xlim = c(-4, 4), ylim = c(0.0,y_max), type = "l",
           lwd = 1, lty = 2, col = "grey30",
           main = "",
           axes = F, xlab="", ylab ="")
      mtext(amoroso_name, side=3, font=2, cex=1.5, line=1)
      axis(1, at = c(-4, 0, 4), labels = c(-4, 0, 4))
      lines(amo_xx, amo_yy, col = "blueviolet", lwd = 2)
      rug(dat, col = "dodgerblue3", lwd = 1)
      
      # Plot adj. KDE, constraint = "unimodal"
      plot(x_range, dnorm(x_range), xlim = c(-4, 4), ylim = c(0.0,y_max), type = "l",
           lwd = 1, lty = 2, col = "grey30",
           main = "",
           axes = F, xlab="", ylab ="")
      mtext("Adj. KDE ('unimodal')", side=3, font=2, cex=1.5, line=1)
      axis(1, at = c(-4, 0, 4), labels = c(-4, 0, 4))
      lines(scKDE_uni$x, scKDE_uni$y, col = "blueviolet", lwd = 2)
      rug(dat, col = "dodgerblue3", lwd = 1)
      
      # Plot adj. KDE, constraint = "twoInflections"
      plot(x_range, dnorm(x_range), xlim = c(-4, 4), ylim = c(0.0,y_max), type = "l",
           lwd = 1, lty = 2, col = "grey30",
           main = "",
           axes = F, xlab="", ylab ="")
      mtext("Adj. KDE ('twoInflections')", side=3, font=2, cex=1.5, line=1)
      axis(1, at = c(-4, 0, 4), labels = c(-4, 0, 4))
      lines(scKDE_2inf$x, scKDE_2inf$y, col = "blueviolet", lwd = 2)
      rug(dat, col = "dodgerblue3", lwd = 1)
      
      # Plot adj. KDE, constraint = "twoInflections+"
      plot(x_range, dnorm(x_range), xlim = c(-4, 4), ylim = c(0.0,y_max), type = "l",
           lwd = 1, lty = 2, col = "grey30",
           main = "",
           axes = F, xlab="", ylab ="")
      mtext("Adj. KDE ('twoInflections+')", side=3, font=2, cex=1.5, line=1)
      axis(1, at = c(-4, 0, 4), labels = c(-4, 0, 4))
      lines(scKDE_2infplus$x, scKDE_2infplus$y, col = "blueviolet", lwd = 2)
      rug(dat, col = "dodgerblue3", lwd = 1)
      
      # Add an overall title
      if (is.null(main)) {
        big_title <- paste0("rnorm(", as.character(n), ")")
      } else {
        big_title <- main
      }
      mtext(text=big_title,side=2,cex=3,line=1.5,outer=TRUE, font=2)
    }
    
  }
  
}

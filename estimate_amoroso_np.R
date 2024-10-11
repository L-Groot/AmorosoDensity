## GOAL: For a given data vector, plot 6 panels next to each other:
# - one per method: R KDE, Bernstein, lowest BIC Amoroso,
# adj KDE (constraint = "unimodal"), adj KDE (constraint = "twoInflections")
# adj KDE (constraint = "twoInflections+")
# - Each plot should show the true distribution and the density estimator
# - The main title above all titles should be the sample size

#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

if (!requireNamespace("scdensity", quietly = TRUE)) {
  install.packages("scdensity")}
library(scdensity)

# Load functions
source("https://raw.githubusercontent.com/L-Groot/AmorosoDensity/refs/heads/main/estimate_amoroso.R")
#source("https://raw.githubusercontent.com/L-Groot/AmorosoDensity/refs/heads/main/estimate_bernstein.R")

#-------------------------------------------------------------------------------

estimate_amoroso_np <- function(dat, main = NULL, minimal = TRUE, plot = TRUE) {
  
  # Get n
  n <- length(dat)
  
  ## Estimate density with different methods ##
  # Amoroso
  amo <- estimate_amoroso(dat, plot=0)
  # Bernstein
  bern1 <- estimate_bernstein(dat, bound_type = "sd")
  bern2 <- estimate_bernstein(dat, bound_type = "Carv")
  # Adjusted KDE
  scKDE_2infplus <- scdensity(dat, constraint = "twoInflections+")
  scKDE_2inf <- scdensity(dat, constraint = "twoInflections")
  scKDE_uni <- scdensity(dat, constraint = "unimodal")
  
  # Extract parameters to plot lowest-BIC Amoroso
  amo_xx <- amo$x
  amo <- amo$win_model
  amo_pars <- as.vector(unlist(amo[,3:6]))
  amoroso_name <- paste0(amo$method, " (", amo$space, ")")
  # Get Amoroso density values
  amo_yy <- dgg4(amo_xx, amo_pars[1], amo_pars[2], amo_pars[3], amo_pars[4])
  # Replace NAs by zeroes
  amo_yy[is.na(amo_yy)] <- 0
  # Remove NA
  # non_na_indices <- !is.na(amo_yy)
  # amo_xx <- amo_xx[non_na_indices]
  # amo_yy <- amo_yy[non_na_indices]
  
  # Define the range of x values
  x_range <- seq(-3, 3, length.out = 1000)
  
  # Define max of y values
  y_max <- max(density(dat)$y,amo_yy,bern1$y,bern2$y,
               scKDE_2infplus$y,scKDE_2inf$y,scKDE_uni$y) + 0.05
  
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
    
    #mtext(big_title, outer = TRUE, cex = 1.5)
    #mtext(big_title, side = 2, outer = TRUE, cex = 1.5, line = 2, font = 2, adj = 2)
    #limits <- par("usr")
    #text(x = limits[1], y = grconvertY(0.5, from = "ndc"),
    #     labels = "color of line", xpd = NA, srt = 270)
    mtext(text=big_title,side=2,cex=3,line=1.5,outer=TRUE, font=2)
  }
}

set.seed(68)
dat <- rnorm(20)
compare_amoroso_np(dat, main = "N = 20")

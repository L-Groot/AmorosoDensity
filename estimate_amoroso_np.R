#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

### Load functions and packages ###

# -> for estimating Amoroso
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoDensity/refs/",
              "heads/main/estimate_amoroso.R"))

# -> for estimating Bernstein
source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoDensity/refs/",
              "heads/main/estimate_bernstein.R"))

# -> for estimating adjusted KDE
if (!requireNamespace("scdensity", quietly = TRUE)) {
  install.packages("scdensity")}
library(scdensity)

# -> for estimating Amoroso
if (!requireNamespace("AmoRosoDistrib", quietly = TRUE)) {
  install.packages("AmoRosoDistrib")}
library(AmoRosoDistrib)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

estimate_amoroso_np <- function(dat, main = NULL, plot = TRUE, minimal = FALSE,
                                generatedbynormal = FALSE, withmean = 0, withsd = 1,
                                amorosocrit = "ML", xticks = NULL) {
  
  
  #############################
  ### 1. ESTIMATE DENSITIES ###
  #############################
  
  # Remove NA
  num_nas_to_remove <- length(dat) - length(na.omit(dat))
  dat <- as.vector(na.omit(dat))
  
  # Print how many NAs were removed
  if(num_nas_to_remove > 0) {
    message(c("WARNING: ", as.character(num_nas_to_remove),
              " NAs were removed from the data.\n"))
    cat("--------------------------------------------------------------------\n")
  }
  
  # Get n
  n <- length(dat)
  
  
  #### Amoroso ####
  
  amo <- estimate_amoroso(dat, plot=0, criterion="maxL")
  # Extract Amoroso parameters
  amo_xx <- amo$x
  if(amorosocrit == "ML") {
    amo_win <- amo$max_L_model
  } else if (amorosocrit == "BIC") {
    amo_win <- amo$min_BIC_model
  } else {
    stop("'amorosocrit' must be either 'ML' or 'BIC'.")
  }
  amo_pars <- as.vector(unlist(amo_win[,3:6]))
  amo_name <- paste0(amo_win$method, " (", amo_win$space, ")")
  # Get Amoroso density values
  amo_yy <- dAmoroso(amo_xx, amo_pars[1], amo_pars[2], amo_pars[3], amo_pars[4])
  #amo_yy <- dgg4(amo_xx, amo_ML_pars[1], amo_ML_pars[2], amo_ML_pars[3], amo_ML_pars[4])
  # Replace NAs by zeroes
  amo_yy[is.na(amo_yy)] <- 0
  # Put Amoroso density values in a list
  amo <- list(x = amo_xx, y = amo_yy)
  
  
  #### Bernstein ####
  
  bern1 <- estimate_bernstein(dat, bound_type = "sd")
  bern2 <- estimate_bernstein(dat, bound_type = "Carv")
  
  
  #### Adjusted KDE ####
  
  scKDE_2infplus <- scdensity(dat, constraint = "twoInflections+")
  scKDE_2inf <- scdensity(dat, constraint = "twoInflections")
  scKDE_uni <- scdensity(dat, constraint = "unimodal")
  
  
  ##### R density ####
  
  rdens <- density(dat)
  
  
  
  # Make list with all models
  modlist <- list(rdens = rdens, amo = amo, bern1 = bern1, bern2 = bern2,
                  scKDE_uni = scKDE_uni,scKDE_2inf = scKDE_2inf,
                  scKDE_2infplus = scKDE_2infplus)
  
  # Define x range for plots
  xmin <- min(sapply(modlist, function(mod) min(mod$x)))
  xmax <- max(sapply(modlist, function(mod) max(mod$x)))
  x_range <- seq(xmin,xmax,length.out = 1000)
  
  # Define y max for plots (highest density value of all methods)
  ymax <- max(density(dat)$y,amo_yy,bern1$y,bern2$y,
              scKDE_2infplus$y,scKDE_2inf$y,scKDE_uni$y) #+ 0.05
  buffer <- 0.15*ymax
  ymax <- ymax + buffer
  
  
  
  #####################
  ### 2. MAKE PLOTS ###
  #####################
  
  if (plot == TRUE) {
    
    #---------------------------------------------------------------------------
    # Non-minimal style (i.e., not for proposal)
    #---------------------------------------------------------------------------
    
    if (minimal == FALSE) {
      
      # Make main titles
      titles <- c("R density() KDE", "Amoroso",
                  "Bernstein Polynomials", "Bernstein Polynomials",
                  "Adj. KDE ('unimodal')","Adj. KDE ('twoInflections')",
                  "Adj. KDE ('twoInflections+')")
      
      # Initialize 2x3 plotting grid
      par(mfrow=c(2,3), oma = c(0, 0, 5, 0), cex.axis = 0.9, font.lab = 2,
          font.axis = 1, family = "Times New Roman")
      
      #.......................
      # Plot density estimates
      #.......................
      for (i in 1:length(modlist)) {
        
        if (names(modlist)[i] != "bern2") {
          
          # Empty plot
          plot(NA, xlim = c(xmin, xmax), ylim = c(0.0, ymax), xlab = "x",
               ylab = "Density", main = titles[i], axes = FALSE)
          ifelse(is.null(xticks),
                 axis(1),
                 axis(1, at = xticks, labels = xticks))
          axis(2, las = 2)
          rug(dat, col = "blue", lwd = 1)
          
          # Add density estimate
          lines(modlist[[i]]$x, modlist[[i]]$y, col = 'mediumorchid2', lwd = 2)
          
          # Optional: add data-generating normal
          if (generatedbynormal == TRUE) {
            lines(x_range, dnorm(x_range, mean = withmean, sd = withsd), 
                  type = "l", lwd = 1, lty = 2, col = "grey30")
          }
          
        } else {
          # Add 2nd Bernstein fit
          lines(modlist[[i]]$x, modlist[[i]]$y, col = 'chartreuse4', lwd = 2)
          # Add legend
          legend("topright", legend = c("bound.type = 'sd'","bound.type = 'Carv'"),
                 col = c("mediumorchid2","chartreuse4"), lty = 1, lwd = 2, cex = 0.8,
                 bty = "n")
        }
      }
      
      #...............
      # Add big title
      #...............
      if (is.null(main)) {
        big_title <- "Nonparametric vs Amoroso fits"
      } else {
        big_title <- main
      }
      mtext(big_title, outer = TRUE, cex = 1.5, line = 2, font = 2)
      
      
      #---------------------------------------------------------------------------
      # Minimal style (for proposal)
      #---------------------------------------------------------------------------
      
    } else {
      
      # Make main titles
      titles <- c("R density()", "Amoroso", "Bernstein", "Bernstein",
                  "Adj. KDE ('unimodal')","Adj. KDE ('twoInflections')",
                  "Adj. KDE ('twoInflections+')")
      
      # Initialize 2x3 plotting grid
      par(mfrow=c(2,3), oma = c(1, 6, 1, 1), mar = c(5,5,5,5), cex.axis = 1.4,
          font.lab = 2, font.axis = 1, family = "Times New Roman")
      
      #.......................
      # Plot density estimates
      #.......................
      for (i in 1:length(modlist)) {
        
        if (names(modlist)[i] != "bern2") {
          
          # Empty plot
          plot(NA, xlim = c(xmin, xmax), ylim = c(0.0,ymax), type = "l",
               lwd = 1, lty = 2, main = "",
               axes = F, xlab="", ylab ="")
          ifelse(is.null(xticks),
                 axis(1),
                 axis(1, at = xticks, labels = xticks))
          rug(dat, col = "dodgerblue3", lwd = 1)
          mtext(titles[i], side=3, font=2, cex=1.5, line=1)
          
          # Add density estimate
          lines(modlist[[i]]$x, modlist[[i]]$y, col = 'blueviolet', lwd = 2)
          
          # Optional: add data-generating normal
          if (generatedbynormal == TRUE) {
            lines(x_range, dnorm(x_range, mean = withmean, sd = withsd), 
                  type = "l", lwd = 1, lty = 2, col = "grey30")
          }
          
        } else {
          
          lines(modlist[[i]]$x, modlist[[i]]$y, col = "deeppink2", lwd = 2)
          # Add legend
          # legend("topright", legend = c("bound.type = 'sd'","bound.type = 'Carv'"),
          #        col = c("blueviolet","deeppink2"), lty = 1, lwd = 2, cex = 0.8,
          # bty = "n")
        }
      }
      
      #........................
      # Add title at left side
      #........................
      if (is.null(main)) {
        big_title <- paste0("rnorm(", as.character(n), ")")
      } else {
        big_title <- main
      }
      mtext(text=big_title, side=2, cex=3, line=1.5, outer=TRUE, font=2)
    }
    
  } else {
    print("no plot (since plot = FALSE")
  }
  
  
  ############################
  ### 3. RETURN MODEL LIST ###
  ############################
  
  invisible(modlist)
  
}


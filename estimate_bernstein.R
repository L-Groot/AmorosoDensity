# Source basis function for estimating Bernstein from Turnbull & Ghosh (2014)
source("https://raw.githubusercontent.com/L-Groot/AmorosoDensity/refs/heads/main/umdz_Rfunction.R")

#-------------------------------------------------------------------------------
# Make estimate_bernstein() function
#-------------------------------------------------------------------------------

estimate_bernstein <- function(dat, breaks = 20, plot = FALSE, n = 512,
                              main = "Bernstein Polynomial Fit",
                              bound_type = "sd",
                              print_busy = FALSE) {
  
  # Extract vector name
  vecname <- paste0("'",deparse(substitute(dat)),"'")
    
  # Print busy statement
  if (print_busy == TRUE) {
    cat("\nBusy estimating Bernstein...\nUsing the provided vector",
        vecname,"'\n",
        "\n--------------------------------------------------------------------\n\n") 
  }
  
  # Remove any missing values
  num_nas_to_remove <- length(dat) - length(na.omit(dat))
  dat <- na.omit(dat)
  
  # Print how many NAs were removed
  if(num_nas_to_remove > 0) {
    message(c("WARNING: ", as.character(num_nas_to_remove), " NAs were removed from the data.\n"))
    cat("--------------------------------------------------------------------\n")
  }
  
  # Calculate min and max value of density() fit
  min_x <- min(density(dat)$x)
  max_x <- max(density(dat)$x)
  
  # Estimate Bernstein
  res <- umd(dat, bound.type = bound_type)
  x_vals <- seq(min_x, max_x, length.out = n)
  y_vals <- res$dumd(x_vals)
  
  if (plot == TRUE) {
    # Plot histogram
    hist(dat, breaks = breaks, freq = FALSE, main = main,
         xlab = "Value", ylab = "Density", border = F)
    # Add estimated density line
    lines(x_vals, y_vals, col = "deeppink2", lwd = 4)
    # Add box around plot
    box()
  }
  
  # Return object
  return(invisible(
    list(x = x_vals,
         y = y_vals
    )))
}


# data <- rnorm(250)
# par(mfrow=c(1,2))
# estimate_bernstein(data, plot = TRUE, bound_type = "sd")
# estimate_bernstein(data, plot = TRUE, bound_type = "Carv")

#dat <- palmerpenguins::penguins$bill_depth_mm
#bern_res <- estimate_bernstein(dat, plot = TRUE)



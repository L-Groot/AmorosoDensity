setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("umdz_Rfunction.R")

estimate_bernstein <- function(dat, breaks = 20, plot = FALSE,
                              main = "Bernstein Polynomial Fit",
                              bound_type = "sd") {
  
  # Calculate min and max value of density() fit
  min_x <- min(density(dat)$x)
  max_x <- max(density(dat)$x)
  
  # Estimate Bernstein
  res <- umd(dat, bound.type = bound_type)
  x_vals <- seq(min_x, max_x, length.out = 1000)
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
  return(list(x = x_vals, y = y_vals))
}


# data <- rnorm(250)
# par(mfrow=c(1,2))
# estimate_bernstein(data, plot = TRUE, bound_type = "sd")
# estimate_bernstein(data, plot = TRUE, bound_type = "Carv")





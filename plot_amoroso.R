# Load packages
ifelse(!require(stats),
       install.packages("stats"),library(stats))
ifelse(!require(tidyverse),
       install.packages("tidyverse"),library(tidyverse))


# Define the dAmoroso function
dAmoroso <- function(x, a, l, c, mu) {
  c1 <- 1/(gamma(l))
  c2 <- abs(c/a)
  c3 <- ((x-mu)/a)^(l*c-1)
  c4 <- exp(-((x-mu)/a)^c)
  return(c1*c2*c3*c4)
}

# Make function that plots Amoroso
plot_amoroso <- function(x_vals, a, l, c, mu,
                         col = "blue", minimal = TRUE,
                         ymax = 1.0,
                         main = "", cex.main = 1.9) {
  
  # Alpha can't be 0
  if(a==0) {
    stop("alpha cannot be 0")
    
    # For alpha > 0
  } else if (a > 0){
    density_values <- c()
    for (x in x_vals) {
      if(x < mu) { # assign density of 0 to all x smaller than mu
        density_values <- append(density_values,0)
      } else {
        dval <- dAmoroso(x, a, l, c, mu)
        density_values <- append(density_values, dval)
      }
    }
    
    # For alpha < 0
  } else if (a < 0) {
    density_values <- c()
    for (x in x_vals) {
      if(x > mu) { # assign density of 0 to all x greater than mu
        density_values <- append(density_values,0)
      } else {
        dval <- dAmoroso(x, a, l, c, mu)
        density_values <- append(density_values, dval)
      }
    }
  }
  
  
  # Plot the density values against the x_vals
  if (minimal == FALSE) {
    
    plot(x_vals, density_values, type = "l", col = col,
         main = main,
         xlab = "", ylab = "",
         lwd = 2,
         ylim = c(0,ymax),
         axes = FALSE)
    
    title(paste0("Amoroso","(",
                 "α=", a, ", ",
                 "\u2113=", l, ", ",
                 "c=", c, ", ",
                 "μ=", mu,
                 ")"
                 
    ), font = 3, cex.main = cex.main)
    
    axis(1, cex.axis = 0.8)
    axis(2, cex.axis = 0.8)
    par(las = 0)
    mtext("x", side = 1, line = 2.2, cex = 0.8, font = 1)
    mtext("Density", side = 2, line = 2.6, cex = 0.8, font = 1)
    
  } else { # minimal = TRUE
    
    plot(x_vals, density_values, type = "l", col = col,
         main = main,
         xlab = "", ylab = "",
         lwd = 2,
         ylim = c(0,ymax),
         axes = FALSE)
    
    title(paste0("α = ", a, ", ",
                 "\u2113 = ", l, ", ",
                 "c = ", c, ", ",
                 "μ = ", mu),
          font = 3, cex.main = cex.main, line = 0)
    
    axis(1, cex.axis = 1.4)
    par(las = 0)
    #mtext("x", side = 1, line = 2.2, cex = 1.4, font = 2)
  }
  
}

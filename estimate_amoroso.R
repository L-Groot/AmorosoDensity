# Load packages
if (!requireNamespace("AmoRosoDistrib", quietly = TRUE)) {
  install.packages("AmoRosoDistrib")}
library(AmoRosoDistrib)

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")}
library(dplyr)

#-------------------
# Define Amoroso PDF
#-------------------

dAmoroso <- function(x, a, lambda, c, mu) {
  c1 <- 1/(gamma(lambda))
  c2 <- abs(c/a)
  c3 <- ((x-mu)/a)^(lambda*c-1)
  c4 <- exp(-((x-mu)/a)^c)
  return(c1*c2*c3*c4)
}

#-----------------------------------------------------
# Define Amoroso negative log likelihood (LL) function
#-----------------------------------------------------

get_negLL_amo <- function(data, params) {
  a <- params[1]
  lambda <- params[2]
  c <- params[3]
  mu <- params[4]
  
  negLL <- sum(-log(dAmoroso(data, a, lambda, c, mu)))
  
  return(negLL)
}

#----------------------------
# Define Amoroso BIC function
#----------------------------

# --> calculates BIC for a certain Amoroso distribution fit to a certain data
get_BIC_amoroso <- function(data, params) {
  N <- length(data)
  P <- length(params)
  LL <- -get_negLL_amo(data, params)
  BIC <- LL + log(N) * (P + 1)
  return(BIC)
}

#-------------------------------------------------------------------
# Define short info texts about plots for estimate_amoroso_MLE_MDE()
#-------------------------------------------------------------------

access_info <- "\nHOW TO ACCESS FULL OUTPUT:\n
  To access all of the output of the estimate_amoroso() function, you have to
  assign the results to an object like this:\n
  res <- estimate_amoroso(my_vector)\n
  Then you can access the following:\n
  - res$all_models: tibble that contains all 18 estimated Amorosos (9 methods
  x 2 parameter spaces (- and +)\n
  - res$min_BIC_models: tibble that contains, for each estimation method, the Amoroso
  fit with the lower BIC (i.e., either in -ve or +ve parameter space)\n
  - res$max_L_models: tibble that contains, for each estimation method, the Amoroso
  fit with the higher likelihood (i.e., either in -ve or +ve parameter space)\n
  - res$min_BIC_model: single-row tibble that contains the Amoroso fit with the lowest
  BIC overall\n
  - res$max_L_model: single-row tibble that contains the Amoroso fit with the highest
  likelihood overall\n
--------------------------------------------------------------------\n"

one_plot_info <- "\nABOUT THE PLOT:\n
  The plot contains the histogram of the variable (grey bars),
  the nonparametric R Kernel density estimator fit (dark grey line) and the
  Amoroso fits from the initial parameter estimation (as described by Combes
  et al. (2022). Each coloured line corresponds to the Amoroso fit of one method.
  Since each method has two sets of parameter estimates (+ve and -ve space), the
  plot shows only the one that fits the data better according to the selected
  criterion (default 'BIC' for lower BIC; 'maxL' for higher likelihood). For
  some methods this may be the parameter set in negative space and for other
  methods this may be the parameter set in positive space.\n"
grid_plots_info <- "\nABOUT THE PLOTS:\n
  Each of the 9 plots shows the two Amoroso fits for each method: the fit from
  the parameter estimates in positive space (green line) and the fit from the
  parameter estimates in negative space (red line). The dark-grey line is the
  fit of the nonparametric R Kernel density.\n"
best_plot_info <- "\nABOUT THE PLOT:\n
  The plot shows the best Amoroso overall, selected either by minimizing BIC
  (criterion == 'BIC) or maximizing likelihood (criterion == 'maxL'). More
  specifically, after fitting the Amoroso with each method described in
  Combes et al. (2022) in both parameter spaces (+ve and -ve), the likelihood/BIC
  was calculated for each of the 18 resulting Amoroso fits (9 methods x
  2 parameter spaces). Then, of all the 18 fits, the one with the lowest BIC/
  highest likelihood is selected as the 'best' fit.\n"


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

estimate_amoroso <- function(vec=NULL, dataframe=NULL, variable=NULL,
                             criterion = "BIC",
                             plot = 1,
                             breaks = 20,
                             varname = NULL,
                             include.init = FALSE) {
  
  # Provide either a datavector for 'vec' OR
  # provide the dataframe to 'dataframe' and the column name of the variable as
  # a string to 'variable.
  
  # criterion: criterion for model selection
  # ("LL" for maximum log likelihood or "BIC" for minimum BIC)
  # criterion is irrelevant when plot == 2
  
  # plot = 0: no plot
  # plot = 1: plot Amorosos of all estimation methods (but only the lower BIC one)
  # plot = 2: 3x3 grid with all methods in separate plot (both parameter space fits)
  # plot = 3: plot with only Amoroso with lowest BIC overall
  
  # breaks = n; n is the number of breaks in the histogram of the data
  
  # varname = "custom name of the variable for plot title"
  
  # include.init = TRUE/FALSE: compatible with plot = 1, you can chose to either
  # add the 2 Amorosos from the initial parameter estimates (TRUE) to the plot
  # or not (FALSE)
  
  
  #####################
  ## Define data (x) ##
  #####################
  
  if (!is.null(vec) && is.vector(vec)) {
    # If vector is provided directly, use it as y
    x <- na.omit(vec)
    cat("\n\nBusy estimating Amoroso...\nUsing the provided vector",
        paste0("'",deparse(substitute(vec)),"'\n"),
        "\n--------------------------------------------------------------------\n\n")
    # Check if any NAs were removed
    num_nas_removed <- length(vec) - length(na.omit(vec))
    if(num_nas_removed > 0) {
      message(c("WARNING: ", as.character(num_nas_removed), " NAs were removed from the data.\n"))
      cat("--------------------------------------------------------------------\n")
    }
  } else if (!is.null(dataframe) && !is.null(variable) && is.data.frame(dataframe) && is.character(variable)) {
    # If dataframe and variable are provided, extract variable from df and assign to y
    x <- na.omit(eval(substitute(dataframe[[variable]])))
    cat("Using the", variable, "variable from the", deparse(substitute(dataframe)), "dataframe. \n\n")
    # Check if any NAs were removed
    num_nas_removed <- length(eval(substitute(dataframe[[variable]]))) - length(na.omit(eval(substitute(dataframe[[variable]]))))
    if(num_nas_removed > 0) {
      message(num_nas_removed, " NAs were removed from the data.\n\n")
    }
  } else {
    stop("Invalid arguments. \n. Provide either \n (1) a vector object directly to the
         `vec` argument or \n (2) a dataframe object to the 'dataframe argument
         AND the name of the variable's column name as a string to the 'variable'
         argument. \n Otherwise check whether the objects, columns etc. you wish to use
         indeed exist and that spelling is correct. \n")
  }
  
  length_x <- length(x) # nr observations
  
  
  ######################################################################
  ## Estimate parameters for MLE and MDE methods in +ve and -ve space ##
  ######################################################################
  
  # Find initializing parameters
  # -> For a > 0
  init.a.pos = init.theta(data = x, -20, 20, length = 1000, a.pos = TRUE)
  init.pos = init.a.pos[1:4]
  # -> For a < 0
  init.a.neg = init.theta(data = x, -20, 20, length = 1000, a.pos = FALSE)
  init.neg = init.a.neg[1:4]
  # MLE
  mleresp <- fit.mle(x, init.pos)
  mleresn <- fit.mle(x, init.neg) 
  # Kullback-Leibler divergence (based on PDF)
  mkleresp <- fit.mkle(x, init.pos)
  mkleresn <- fit.mkle(x, init.neg) 
  # Jensen-Shanon divergence (based on PDF)
  mjseresp <- fit.mjse(x, init.pos) 
  mjseresn <- fit.mjse(x, init.neg) 
  # Hellinger distance (based on PDF)
  mheresp <- fit.mhe(x, init.pos) 
  mheresn <- fit.mhe(x, init.neg) 
  # Wasserstein distance (based on PDF)
  mWassderesp <- fit.mwe(x, init.pos, d = 4)   
  mWassderesn <- fit.mwe(x, init.neg, d = 4)  
  # Squared distance (based on PDF)
  msqeresp <- fit.msqe(x, init.pos) 
  msqeresn <- fit.msqe(x, init.neg) 
  # Hellinger distance (based on CDF)
  mhecresp <- fit.mhdfe(x, init.pos) 
  mhecresn <- fit.mhdfe(x, init.neg) 
  # Wasserstein distance (based on CDF)
  mwecresp <- fit.mwdfe(x, init.pos, d = 4)    
  mwecresn <- fit.mwdfe(x, init.neg, d = 4)   
  # Squared distance (based on CDF)
  msqdferesp <- fit.msqdfe(x, init.pos) 
  msqdferesn <- fit.msqdfe(x, init.neg)
  # PUT ALL IN A LIST
  fit_list <- list(mleresp, mleresn, # MLE
                   mkleresp, mkleresn, # Kullback-Leibler (PDF)
                   mjseresp, mjseresn, # Jensen-Shanon (PDF)
                   mheresp, mheresn, # Hellinger (PDF)
                   mWassderesp, mWassderesn, # Wasserstein (PDF)
                   msqeresp, msqeresn, # Squared (PDF)
                   mhecresp, mhecresn, # Hellinger (CDF)
                   mwecresp, mwecresn, # Wasserstein (CDF)
                   msqdferesp, msqdferesn # Squared (CDF)
  )
  
  ##############################################
  ## Make dataframe with all estimated models ##
  ##############################################
  
  names_vec <- rep(c("MLE", "Kullback-Leibler PDF", "Jensen-Shanon PDF",
                     "Hellinger PDF", "Wasserstein PDF", "Squared PDF",
                     "Hellinger CDF", "Wasserstein CDF", "Squared CDF"), each = 2)
  
  names_ID <- rep(c("MLE", "KL-PDF", "JS-PDF", "Hell-PDF", "WASS-PDF", "SQU-PDF",
                     "HELL-CDF", "WASS-CDF", "SQU-CDF"), each = 2)
                     
  
  all_models_df <- data.frame(method = names_vec, method_ID = names_ID)
  
  par_space_vec <- rep(c("+", "-"), times = length(names_vec)/2)
  a_vec <- sapply(fit_list, function(fit) fit$par[1])
  l_vec <- sapply(fit_list, function(fit) fit$par[2])
  c_vec <- sapply(fit_list, function(fit) fit$par[3])
  mu_vec <- sapply(fit_list, function(fit) fit$par[4])
  BIC_vec <- sapply(fit_list, function(fit) get_BIC_amoroso(x, fit$par))
  negLL_vec <- sapply(fit_list, function(fit) get_negLL_amo(x, c(fit$par)))
  
  all_models_tib <- all_models_df %>% mutate(
    space = par_space_vec,
    a = a_vec,
    l = l_vec,
    c = c_vec,
    mu = mu_vec,
    negLL = negLL_vec,
    BIC = BIC_vec) %>%
    as_tibble()
  
  
  ###################################################
  ## Make dataframe of lowest BIC model PER METHOD ##
  ###################################################
  
  min_BIC_models_tib <- all_models_tib %>%
    group_by(method) %>%
    filter(BIC == min(BIC)) %>%
    ungroup()
  
  ##################################################
  ## Make dataframe of higher LL model PER METHOD ##
  ##################################################
  
  max_L_models_tib <- all_models_tib %>%
    group_by(method) %>%
    filter(negLL == min(negLL)) %>%
    ungroup()
  
  ###############################################
  ## Make dataframe of lowest BIC model OF ALL ##
  ###############################################
  
  min_BIC_model_tib <- min_BIC_models_tib %>%
    filter(BIC == min(BIC))
  
  #########################################
  ## Make dataframe of highest LL OF ALL ##
  #########################################
  
  max_L_model_tib <- max_L_models_tib %>%
    filter(negLL == min(negLL))
  
  
  ###########
  ## PLOTS ##
  ###########
  
  # Define x
  xx <- seq(min(density(x)$x), max(density(x)$x), length = 10000)
  
  # Calculate KDE
  dens <- density(x)
  
  #---------#
  # NO PLOT #
  #---------#
  
  if (plot == 0) {
    cat("no plot",
        "\n--------------------------------------------------------------------\n\n")
    
    #---------------------------------------#
    # One plot with best Amoroso per method #
    #---------------------------------------#
    
  } else if (plot == 1) {
    
    # Choose model selection criterion
    if(criterion == "BIC") {
      win_models_tib <- min_BIC_models_tib
    } else if (criterion == "maxL") {
      win_models_tib <- max_L_models_tib
    } else {
      print("ERROR 1")
    }
    # Print plot info in console
    cat(one_plot_info,
        "\n--------------------------------------------------------------------\n\n")
    # Parameter settings
    par(cex.main = 1.4, cex.axis = 1, cex.lab = 1.2, bty = "n", font.lab = 2)
    # Initialize title
    if (!is.null(varname)) {
      title <- paste0("Amoroso fits to '", varname, "'")
    } else if (is.null(variable) && !is.null(vec)) {
      title <- paste0("Amoroso fits to '", deparse(substitute(vec)), "'")
    } else {
      title <- paste0("Amoroso fits to '", variable, "'")
    }
    # Create histogram
    hist(x, prob = T, main = title, breaks = breaks, col = "grey95",
         border = "grey85", axes = FALSE)
    axis(1)
    axis(2, las = 1)
    # Add kernel density estimate
    lines(dens$x, dens$y, col = "grey70", lty = 1, lwd = 2)
    # Either don't include the initial estimates...
    if (include.init == FALSE) {
      legend_names <- c("R density()")
      legend_colors <- c("grey70")
      lty_vec <- c(1, rep(1, nrow(win_models_tib)))
    } else {
      # ... or include them
      points(xx, dgg4(xx, init.pos[1], init.pos[2], init.pos[3], init.pos[4]), type = "l", lwd = 2, col = "black", lty = 1)
      points(xx, dgg4(xx, init.neg[1], init.neg[2], init.neg[3], init.neg[4]), type = "l", lwd = 2, col = "black", lty = 2)
      legend_names <- c("Initial Est.(+)", "Initial Est.(-)", "R density()")
      legend_colors <- c("black", "black", "grey70")
      lty_vec <- c(1, 2, 1, rep(1, nrow(win_models_tib)))
    }
    # Define colors for each model
    model_colors <- rainbow(nrow(win_models_tib))
    # Loop to add each winning model
    for (i in 1:nrow(win_models_tib)) {
      # Create name for legend
      name <- paste0(win_models_tib$method[i], " (", win_models_tib$space[i], ")")
      # Get model parameters
      pars <- c(win_models_tib$a[i],win_models_tib$l[i],win_models_tib$c[i],win_models_tib$mu[i])
      # Plot the winning model with a unique color
      points(xx, dgg4(xx, pars[1], pars[2], pars[3], pars[4]), 
             type = "l", lwd = 2, col = model_colors[i], lty = 1)
      # Add the model name to the legend
      legend_names <- c(legend_names, name)
      legend_colors <- c(legend_colors, model_colors[i])
    }
    # Add legend
    legend("topright", legend = legend_names,
           lwd = 2,
           col = legend_colors,
           lty = lty_vec,
           bty = "n",
           cex = 0.7)
    
    #--------------------------------------------#
    # 3x3 grid of plots with one plot per method #
    #--------------------------------------------#
    
  } else if (plot == 2) {
    
    # Print info about plots
    cat(grid_plots_info,
        "\n--------------------------------------------------------------------\n\n")
    # Make 3x3 grid
    par(mfrow = c(3, 3), cex.axis = 0.8)
    # Loop through methods
    for (i in seq(1, nrow(all_models_tib), by = 2)) {
      # Subset two rows of the method
      meth_tib <- all_models_tib[i:(i+1), ]
      # Extract parameters for model in +ve parameter space
      pos_par <- meth_tib %>%
        slice(1) %>%
        select(a, l, c, mu) %>% 
        unlist(use.names = FALSE)
      # Extract parameters for model in -ve parameter space
      neg_par <- meth_tib %>%
        slice(2) %>%
        select(a, l, c, mu) %>%
        unlist(use.names = FALSE)
      # Make plot title
      title <- meth_tib$method[[1]]
      # Plot histogram
      hist(x, prob = T, main = title, breaks = breaks, col = "grey95",
           border = "grey85", axes = FALSE)
      axis(1)
      axis(2, las = 1)
      # Add KDE
      lines(dens$x, dens$y, col = "grey70", lty = 1, lwd = 1.5)
      # Add a>0 model
      points(xx, dgg4(xx, pos_par[1], pos_par[2], pos_par[3], pos_par[4]),
             type = "l", lwd = 1.5, col = "sienna3", lty = 1)
      # Add a<0 model
      points(xx, dgg4(xx, neg_par[1], neg_par[2], neg_par[3], neg_par[4]),
             type = "l", lwd = 1.5, col = "steelblue3", lty = 1)
      # Add legend
      legend("topright", legend = c("R density()", "α > 0", "α < 0"),
             col = c("grey70", "sienna3", "steelblue3"), lty = c(1, 1, 1),
             lwd = 1.5, bty = "n", cex = 0.8)
    }
    
    #---------------------------#
    # Plot only winning Amoroso #
    #---------------------------#
    
  } else if (plot == 3) {
    
    # Choose model selection criterion
    if(criterion == "BIC") {
      win_model_tib <- min_BIC_model_tib
    } else if (criterion == "maxL") {
      win_model_tib <- max_L_model_tib
    } else {
      print("ERROR 2")
    }
    # Print plot info in console
    cat(best_plot_info,
        "\n--------------------------------------------------------------------\n\n")
    # Set plotting parameters
    par(cex.main = 1.4, cex.axis = 1, cex.lab = 1.2, bty = "n", font.lab = 2)
    # Extract name of lowest BIC Amoroso
    method_name <- paste0(win_model_tib[["method"]]," (",
                          win_model_tib[["space"]], ")")
    # Make main title
    if (!is.null(varname)) {
      title <- paste0(method_name, " makes the lowest-BIC Amoroso fit for '",
                      varname, "'")
    } else if (is.null(variable) && !is.null(vec)) {
      title <- paste0(method_name," makes the lowest-BIC Amoroso fit for '",
                      deparse(substitute(vec)), "'")
    } else {
      title <- paste0(method_name," makes the lowest-BIC Amoroso fit for '",
                      variable, "'")
    }
    # Make histogram
    hist(x, prob = T, main = title, breaks = breaks, col = "grey95",
         border = "grey85", axes = FALSE, bty = "n")
    axis(1)
    axis(2, las = 1)
    # Add KDE line
    lines(dens$x, dens$y, col = "grey70", lty = 1, lwd = 2)
    # Extract parameters of lowest BIC Amoroso
    pars <- win_model_tib%>%
      slice(1) %>%
      select(a, l, c, mu) %>% 
      unlist(use.names = FALSE)
    # Add Amoroso line
    points(xx, dgg4(xx, pars[1], pars[2], pars[3], pars[4]), 
           type = "l", lwd = 4, col = "mediumpurple3", lty = 1)
    # Add legend
    legend_names <- c("R density()", method_name)
    legend_colors <- c("grey70","mediumpurple3")
    legend("topright", legend = legend_names,
           lwd = c(2,4),
           col = legend_colors,
           lty = c(1, 1),
           bty = "n",
           cex = 1.2)
    
    #-----------------------------------------------#
    # Error message: inadmissible plotting argument #
    #-----------------------------------------------#
    
  } else {
    stop("plot must be one of 0, (no plot),1 (all methods in one plot), 2 (one plot
         per method) or 3 (one plot with lowest BIC Amoroso")
  }
  
  
  ####################
  ## Return objects ##
  ####################
  
  # Print best model info
  cat("MAXIMUM LIKELIHOOD AMOROSO:\n\n",
      paste0(max_L_model_tib$method," (",max_L_model_tib$space,")"),"\n")
  cat("-LL:", max_L_model_tib$negLL, "\n")
  cat(paste0(" α = ", round(max_L_model_tib$a,digits=2), ", ",
             "\u2113 = ", round(max_L_model_tib$l,digits=2), ", ",
             "c = ", round(max_L_model_tib$c,digits=2), ", ",
             "μ = ", round(max_L_model_tib$mu,digits=2)),"\n\n")
  
  cat("LOWEST BIC AMOROSO:\n\n",
      paste0(min_BIC_model_tib$method," (",min_BIC_model_tib$space,")"),"\n")
  cat(" BIC:", min_BIC_model_tib$BIC, "\n")
  cat(paste0(" α = ", round(min_BIC_model_tib$a,digits=2), ", ",
             "\u2113 = ", round(min_BIC_model_tib$l,digits=2), ", ",
             "c = ", round(min_BIC_model_tib$c,digits=2), ", ",
             "μ = ", round(min_BIC_model_tib$mu,digits=2)),"\n")
  cat("\n--------------------------------------------------------------------\n")
  
  # Print info on how to access full output
  cat(access_info, "\n")
  
  # Return all win models and final best model
  return(invisible(list(
    all_models = all_models_tib,
    min_BIC_models = min_BIC_models_tib,
    max_L_models = max_L_models_tib,
    min_BIC_model = min_BIC_model_tib,
    max_L_model = max_L_model_tib, 
    x = xx
  )))
  
}


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#-----------------------------------------
# Test the function
#-----------------------------------------

if (!requireNamespace("palmerpenguins", quietly = TRUE)) {
  install.packages("palmerpenguins")}
library(palmerpenguins)

dat <- palmerpenguins::penguins$flipper_length_mm
#dat <- palmerpenguins::penguins$bill_length_mm
#dat <- palmerpenguins::penguins$bill_depth_mm

#res <- estimate_amoroso(dat, plot = 1, criterion = "BIC")
res <- estimate_amoroso(dat, plot = 1, criterion = "maxL")
 
# res <- estimate_amoroso(dat, plot = 2)

# res <- estimate_amoroso(dat, plot = 3, criterion = "BIC")
# res <- estimate_amoroso(dat, plot = 3, criterion = "maxL")

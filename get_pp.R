#-------------------------------------------------------------------------------
### Functions and packages

source(paste0("https://raw.githubusercontent.com/L-Groot/AmorosoDensity/refs/",
              "heads/main/estimate_amoroso_np.R"))

if (!requireNamespace("AmoRosoDistrib", quietly = TRUE)) {
  install.packages("AmoRosoDistrib")}
library(AmoRosoDistrib)

if (!requireNamespace("caret", quietly = TRUE)) {
  install.packages("caret")}
library(caret)

if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")}
library(tidyverse)

options(scipen = 999)

#-------------------------------------------------------------------------------
### Custom functions

# Mean Squared Error (MSE) function
calculate_mse <- function(pred, true) {
  mean((pred - true)^2)
}

# Mean Absolute Error (MAE) function
calculate_mae <- function(pred, true) {
  mean(abs(pred - true))
}

# Root Mean Square Error (RMSE) function
calculate_rmse <- function(pred, true) {
  sqrt(mean((pred - true)^2))
}


#-------------------------------------------------------------------------------
#res <- get_pp_amo(amo_pars = c(4,1,0.9,0), method = "k-fold", sequence = TRUE)
#-------------------------------------------------------------------------------

### Function that calculates predictive performance of all methods when the
# data-generating distribution is an Amoroso

get_pp_amo <- function(amo_pars = c(-1,2,1,5),
                   n = 1000, # size of dataset to be simulated
                   method = "k-fold", # or "split-half"
                   k = 5, # nr of folds (if method = 'k-fold')
                   prop_train = 0.8, # % data in training set (if method = 'split-half')
                   #plot = FALSE, # plot the fits or not
                   sequence = FALSE, # whether seed should be set for Amoroso data generation
                   seed = 125 # seed for random fold creation
          ) 
{
  
  
  #############################
  ### GENERATE AMOROSO DATA ###
  #############################
  
  true_pars <- amo_pars
  dat <- rgg4(n,
              a = true_pars[1], l = true_pars[2],
              c = true_pars[3], mu = true_pars[4],
              sequence = sequence)
  
  
  ###############
  ### K- FOLD ###
  ###############
  
  if (method == "k-fold") {
    
    #------------------------
    # Split data into k folds
    #------------------------
    set.seed(seed)
    folds <- createFolds(dat, k = k, list = TRUE, returnTrain = TRUE)
    
    # Initialize a tibble with the method names
    error_tib <- tibble(
      Method = c("bern1", "bern2",
                 "scKDE_uni", "scKDE_2inf", "scKDE_2infplus",
                 "rdens", "amo")
    )
    
    #-----------------------
    # Loop through each fold
    #-----------------------
    for (fold in 1:k) {
      
      cat("fold =", fold, "\n")
      
      #------------------------------------------
      # (1) Get training and test data for fold i
      #------------------------------------------
      train_indices <- folds[[fold]]
      train <- dat[train_indices]
      test <- dat[-train_indices]
      
      #------------------------------------------------
      # (2) Fit Amoroso and NP methods on training data
      #------------------------------------------------
      res <- estimate_amoroso_np(train, hist = TRUE, amorosocrit = "ML")
      
      #cat("fold =", i, "(after estimation)\n")
      
      #-------------------------------------------
      # (3) Make continuous functions from NP fits
      #-------------------------------------------
      
      # Make a list to store the functions of the NP fits
      npfun_list <- list(rdens = NULL,
                          bern1 = NULL,
                          bern2 = NULL,
                          scKDE_uni = NULL,
                          scKDE_2inf = NULL,
                          scKDE_2infplus = NULL
                          )
      
      # Loop through list and add the functions
      for (i in 1:length(npfun_list)) { # the first 6 fits in res are the np fits
        if (length(res[[i]]) > 1) { # if the fit is valid, estimate function
          npfun_list[[i]] <- splinefun(res[[i]]$x,res[[1]]$y, method = "monoH.FC")
        } else { # if the fit is not valid, put NA
          npfun_list[[i]] <- NA
        }
      }
      
      #cat("fold =", i, "(after making continuous functions)\n")
      
      #--------------------------------------------------
      # (4) Generate predictions from NP and Amoroso fits
      #--------------------------------------------------
      
      # Make a list to store predictions in
      pred_list <- list(rdens = NULL,
                        bern1 = NULL,
                        bern2 = NULL,
                        scKDE_uni = NULL,
                        scKDE_2inf = NULL,
                        scKDE_2infplus = NULL,
                        amo = NULL
                        )
      
      # Loop through list and generate predictions from each fit
      for (i in 1:(length(pred_list)-1)) { # not for Amoroso
        if(!is.na(npfun_list)[[i]]) { # if the fit is valid
          pred_list[[i]] <- npfun_list[[i]](test) # add predictions
        } else { # if the fit is invalid
          pred_list[[i]] <- NA # add NA
        }
      }
      
      # Make Amoroso predictions
      (est_pars <- res$amo$pars)
      
      pred_y_amo <- dgg4(test, est_pars[1], est_pars[2], est_pars[3], est_pars[4])
      
      if(anyNA(pred_y_amo)) {
        # Count nr of NA values
        nr_na <- sum(is.na(pred_y_amo))
        # Replace NA values with zero
        pred_y_amo[is.na(pred_y_amo)] <- 0
        # Print warning message
        warning(paste(nr_na, "NA values were replaced with zero in pred_y_amo."))
      }
      
      # Add Amoroso predictions to pred_list
      pred_list$amo <- pred_y_amo
      
      #-------------------------------------------
      # (5) Get 'true' density values for test set
      #-------------------------------------------
      true_y <- dgg4(test, true_pars[1], true_pars[2], true_pars[3], true_pars[4])
      
      #-----------------------------
      # (6) Calculate error measures
      #-----------------------------
      mse_values <- c()
      rmse_values <- c()
      mae_values <- c()
      
      names(pred_list)
      
      for (i in 1:length(pred_list)) {
        # If fit is valid (i.e., there are predicted values)
        if (length(pred_list[[i]] > 1)) { 
          # Calculate error measures
          mse_values <- append(mse_values,
                               calculate_mse(pred_list[[i]], true_y))
          rmse_values <- append(rmse_values,
                                calculate_rmse(pred_list[[i]], true_y))
          mae_values <- append(mae_values,
                               calculate_mae(pred_list[[i]], true_y))
        # If fit is not valid
        } else {
          # Add NA
          mse_values <- append(mse_values, NA)
          rmse_values <- append(rmse_values, NA)
          mae_values <- append(mae_values, NA)
        }
      }
      
      # Add error measures of the current fold to the tibble
      error_tib[[paste0("mse_f", fold)]] <- mse_values
      error_tib[[paste0("rmse_f", fold)]] <- rmse_values
      error_tib[[paste0("mae_f", fold)]] <- mae_values
    }
    
    print(error_tib)
    
    #---------------------------------
    # Average the error over the folds
    #---------------------------------
    # Count nr of invalid folds per method
    na_count <- error_tib %>%
      #group_by(Method) %>%
      select(Method, starts_with("mse")) %>%
      mutate(na_count = rowSums(is.na(.))) %>%
      pull(na_count)
    
    # Add column with nr of invalid folds to error_tib
    error_tib <- error_tib %>% mutate(nr_na_folds = na_count)
    
    # Make tibble of average error across folds
    # -> excluding the invalid folds
    av_error_tib <- error_tib %>%
      mutate(
        mse_av = rowMeans(select(., starts_with("mse_")), na.rm = TRUE),
        rmse_av = rowMeans(select(., starts_with("rmse_")), na.rm = TRUE),
        mae_av = rowMeans(select(., starts_with("mae_")), na.rm = TRUE)
      ) %>% select(Method,mse_av,rmse_av,mae_av,nr_na_folds)
    
    
    # av_error_tib <- error_tib %>%
    #   mutate(
    #     mse_av = rowMeans(select(., starts_with("mse_")), na.rm = TRUE),
    #     rmse_av = rowMeans(select(., starts_with("rmse_")), na.rm = TRUE),
    #     mae_av = rowMeans(select(., starts_with("mae_")), na.rm = TRUE)
    #   ) %>% select(Method,starts_with(c("mse_av","rmse_av","mae_av")))
    # 
    #print("av_error_tib")
    #print(av_error_tib)
    
    
    ##################
    ### SPLIT-HALF ###
    ##################
    
  } else if (method == "split-half") {
    
    #------------------------
    # Split data in two parts
    #------------------------
    inTrain <- createDataPartition(
      y = dat,
      p = prop_train, # prop. of data in the training set
      list = FALSE
    )
    train <- dat[inTrain]
    test <- dat[-inTrain]
    
    #--------------------------------------------
    # Fit Amoroso and NP methods on training data
    #--------------------------------------------
    res <- estimate_amoroso_np(train, hist = TRUE, amorosocrit = "ML")
    
    #---------------------------------------
    # Make continuous functions from NP fits
    #---------------------------------------
    f_bern1 <- splinefun(res$bern1$x,res$bern1$y, method = "monoH.FC")
    f_bern2 <- splinefun(res$bern2$x,res$bern2$y, method = "monoH.FC")
    f_scKDE_uni <- splinefun(res$scKDE_uni$x,res$scKDE_uni$y, method = "monoH.FC")
    f_scKDE_2inf <- splinefun(res$scKDE_2inf$x,res$scKDE_2inf$y, method = "monoH.FC")
    f_scKDE_2infplus <- splinefun(res$scKDE_2infplus$x,res$scKDE_2infplus$y, method = "monoH.FC")
    f_rdens <- splinefun(res$rdens$x,res$rdens$y, method = "monoH.FC")
    
    #----------------------------------------------
    # Generate predictions from NP and Amoroso fits
    #----------------------------------------------
    pred_y_bern1 <- f_bern1(test)
    pred_y_bern2 <- f_bern2(test)
    pred_y_scKDE_uni <- f_scKDE_uni(test)
    pred_y_scKDE_2inf <- f_scKDE_2inf(test)
    pred_y_scKDE_2infplus <- f_scKDE_2infplus(test)
    pred_y_rdens <- f_rdens(test)
    
    (est_pars <- res$amo$pars)
    
    pred_y_amo <- dgg4(test, est_pars[1], est_pars[2], est_pars[3], est_pars[4])
    if(anyNA(pred_y_amo)) {
      # Count nr of NA values
      nr_na <- sum(is.na(pred_y_amo))
      # Replace NA values with zero
      pred_y_amo[is.na(pred_y_amo)] <- 0
      # Print warning message
      warning(paste(nr_na, "NA values were replaced with zero in pred_y_amo."))
    }
    
    #-------------------------
    # Calculate error measures
    #-------------------------
    # Create a list of methods and their corresponding predictions
    predictions_list <- list(
      bern1 = pred_y_bern1,
      bern2 = pred_y_bern2,
      scKDE_uni = pred_y_scKDE_uni,
      scKDE_2inf = pred_y_scKDE_2inf,
      scKDE_2infplus = pred_y_scKDE_2infplus,
      rdens = pred_y_rdens,
      Amoroso = pred_y_amo
    )
    
    # Calculate error metrics for each method
    error_on_test_tib <- map_df(names(predictions_list), function(method) {
      pred <- predictions_list[[method]]
      
      tibble(
        method = method,
        mse = calculate_mse(pred, true_y),
        rmse = calculate_rmse(pred, true_y),
        mae = calculate_mae(pred, true_y)
      )
    })
    
    # Print the final error results tibble
    #print("error_on_test")
    #print(error_on_test_tib)
    
    
  } else {
    
    stop("method must be either 'k-fold' or 'split-half'")
    
  }
  
  ################################
  ### PRINT AND RETURN RESULTS ###
  ################################
  
  if(method == "k-fold") {
    cat("Method = 'k-fold', with k =", k, "folds.\n")
    cat("Tibble contains the average MSE, RMSE and MAE of each method over the k folds.\n\n")
    print(av_error_tib)
    if (any(av_error_tib$nr_na_folds > 0)) {
      warning("For at least one method, the average error measures were
              calculated excluding one or more folds. Check the 'nr_na_folds'
              column in the 'av_tib' tibble to see how many folds were excluded
              for which method.")
    }
    return(list(av_error_tib = av_error_tib, error_tib = error_tib))
  } else {
    cat("Method = 'split-half' with", prop_train*100, "% of data in the training set.\n")
    cat("Tibble contains the MSE, RMSE and MAE of each method on the test set.\n\n")
    print(error_on_test_tib)
    return(error_on_test_tib)
  }
}
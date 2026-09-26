# ============================================================================
# STEPWISE BACKWARDS SELECTION
# ============================================================================

optimize_dfm_bic <- function(data, config, maxit = 2000) {
  
  # Initialize tracking variables
  current_data   <- data
  current_config <- config
  
  # Helper function to calculate BIC
  calc_bic <- function(fit, dataset) {
    k   <- length(fit$th_hat)       
    TT  <- nrow(dataset)                      # Sample size needed for BIC
    bic <- (k * log(TT)) - (2 * fit$loglik)   # BIC formula
    return(bic)
  }
  
  # ----------------------------------------------------------------------------
  # RUN THE BASELINE MODEL
  # ----------------------------------------------------------------------------
  cat(sprintf("Fitting baseline model with %d variables...\n", ncol(current_data)))
  
  # Ensure the global config matches our starting point
  assign("indicator_config", current_config, envir = .GlobalEnv)
  
  best_fit <- dfm(current_data, maxit = maxit)
  best_bic <- calc_bic(best_fit, current_data)
  
  cat(sprintf("Baseline BIC: %.2f\n", best_bic))
  cat("--------------------------------------------------\n")
  
  # ----------------------------------------------------------------------------
  # STEPWISE BACKWARDS SELECTION LOOP
  # ----------------------------------------------------------------------------
  repeat {
    
    # Stop if we are down to just GDP and 5 indicators
    if (ncol(current_data) <= 6) {
      cat("Minimum threshold reached (GDP + 5 indicators). Stopping optimization.\n")
      break
    }
    
    # Extract the estimated factor
    factor_est <- best_fit$factor
    
    # Calculate absolute correlation between the factor and each variable
    correlations <- numeric(ncol(current_data))
    for (i in seq_len(ncol(current_data))) {
      correlations[i] <- abs(cor(current_data[, i], 
                                 factor_est, use = "pairwise.complete.obs"))
    }
    
    # Never drop GDP by setting its correlation with the factor to infinity
    correlations[1] <- Inf 
    
    # Identify the variable with the weakest correlation to the factor
    drop_idx  <- which.min(correlations)
    drop_name <- current_config[[drop_idx]]$name
    
    cat(sprintf("Evaluating removal of weakest link: '%s' (Abs Cor: %.3f)...\n", 
                drop_name, correlations[drop_idx]))
    
    # Create temporary dataset and config without this variable
    temp_data   <- current_data[, -drop_idx, drop = FALSE]
    temp_config <- current_config[-drop_idx]
    
    # Update global config for the dfm() function to use
    assign("indicator_config", temp_config, envir = .GlobalEnv)
    
    # Fit the test model
    cat("Fitting new model...\n")
    temp_fit <- dfm(temp_data, maxit = maxit)
    temp_bic <- calc_bic(temp_fit, temp_data)
    
    cat(sprintf("New BIC: %.2f | Old BIC: %.2f\n", temp_bic, best_bic))
    
    # --------------------------------------------------------------------------
    # COMPARE BICs AND DECIDE
    # --------------------------------------------------------------------------
    if (temp_bic < best_bic) {
      # The model improved, confirm the changes and loop again
      cat(sprintf(">> Improvement found! Dropping '%s' permanently.\n", drop_name))
      cat("--------------------------------------------------\n")
      
      current_data   <- temp_data
      current_config <- temp_config
      best_fit       <- temp_fit
      best_bic       <- temp_bic
      
    } else {
      
      # The model got worse, reject the change and break the loop
      cat(sprintf(">> No improvement. The removal of '%s' hurt the model.\n", drop_name))
      cat(">> Optimization complete. Restoring best model.\n\n")
      assign("indicator_config", current_config, envir = .GlobalEnv)
      break
    }
  }
  
  # Return a list containing the optimized model and data
  list(
    fit         = best_fit,
    data        = current_data,
    config      = current_config,
    final_bic   = best_bic
  )
}
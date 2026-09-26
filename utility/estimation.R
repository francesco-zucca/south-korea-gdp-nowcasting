# ==============================================================================
# ESTIMATION FUNCTION
# ==============================================================================

# y is the dataset, maxit is the number of iterations of the likelihood process
dfm <- function(data, maxit = 2000) {
  
  # ----------------------------------------------------------------------------
  # SHIFTING THE DATA IF LEADING
  # ----------------------------------------------------------------------------
  
  # Copy raw data
  y_shifted <- data
  
  # Define the number of indicators
  n <- ncol(data)
  
  # Check if indicator has a lead and lag the data by k periods to make it 
  # coincident
  for (i in seq_len(n)) {
    cfg <- indicator_config[[i]]
    
    # Check if the indicator has a lead time defined and it is greater than 0
    if (!is.null(cfg$lead) && cfg$lead > 0) {
      k <- cfg$lead
      
      # Lag the data by 'k' periods to align it with the coincident factor
      y_shifted[, i] <- c(rep(NA, k), head(data[, i], -k))
    }
  }
  
  # ----------------------------------------------------------------------------
  # STANDARDIZING THE DATA AND CREATING THE NA INDEX
  # ----------------------------------------------------------------------------
  
  # Compute mean and standard deviation, without including the NAs
  y_means <- apply(y_shifted, 2, mean, na.rm = TRUE)
  y_sds   <- apply(y_shifted, 2, sd, na.rm = TRUE)
  
  # Standardizing the data
  y <- scale(y_shifted, center = y_means, scale = y_sds)
  
  # Create NA index matrix
  na_idx <- ifelse(is.na(y), 0, 1)
  
  # Get random draws from a N(0,1)
  random_draws <- matrix(rnorm(prod(dim(y)), mean = 0, sd = 1), 
                         nrow = nrow(y), ncol = ncol(y))
  
  # Replace NAs with random draws
  y[is.na(y)] <- random_draws[is.na(y)]
  
  # Define the number of rows in the dataset
  TT <- nrow(y)
  
  # Variance normalization
  varQ <- 1
  
  # ----------------------------------------------------------------------------
  # GETTING STARTING VALUES AND SETTING BOUNDS FOR CONSTRAINED OPTIMIZATION
  # ----------------------------------------------------------------------------
  
  # Starting values for the parameters
  th0 <- make_startval(y, n)
  
  # Compute the number of parameters
  n_params <- length(th0)
  
  # Initialize all parameters to be unbounded
  lower_bounds <- rep(-Inf, n_params)
  upper_bounds <- rep(Inf, n_params)
  
  # Constraint on the loading of production
  upper_bounds[2] <- 0.55
  
  lower_bounds[(length(th0) - n + 1):length(th0)] <- 0.05 # Variance floor
  
  
  
  # ----------------------------------------------------------------------------
  # AMOUNT OF FACTORS TO BE INCLUDED IN THE LATENT VECTOR
  # ----------------------------------------------------------------------------
  
  # Initialize the vectors
  max_factor_states  <- 0
  total_error_states <- 0
  
  # Get information about the data and compute the number of errors to be 
  # included in the h vector
  for (i in seq_len(n)) {
    cfg <- indicator_config[[i]]
    
    if (cfg$type == "MoM") {
      max_factor_states  <- max(max_factor_states, 1)
      total_error_states <- total_error_states + 2
    } else if (cfg$type == "QoQ") {
      max_factor_states  <- max(max_factor_states, 5)
      total_error_states <- total_error_states + 5
    } else if (cfg$type == "soft") {
      max_factor_states  <- max(max_factor_states, 12)
      total_error_states <- total_error_states + 2
    }
  }
  
  # Compute the number of factors to be included in the h vector
  max_factor_states <- max(max_factor_states, 2)
  
  # Compute the dimension of the h vector
  h_dim <- max_factor_states + total_error_states
  
  # ----------------------------------------------------------------------------
  # LIKELIHOOD MAXIMIZATION
  # ----------------------------------------------------------------------------
  
  # Likelihood maximization
  opt <- optim(
    par    = th0,
    fn     = function(th) 
      kalman(th, y, na_idx, n, TT, h_dim, 
             indicator_config, max_factor_states, varQ)$neg_loglik,
    method = "L-BFGS-B",
    lower = lower_bounds,
    upper = upper_bounds,
    hessian = TRUE,
    control = list(maxit = maxit, factr = 1e7)
  )
  th_hat <- opt$par
  
  # Standard errors estimates
  H      <- opt$hessian
  se_hat <- sqrt(diag(solve(H)))
  
  # Rerunning Kalman filter with optimal parameters to extract the factors
  kf_hat     <- kalman(th_hat, y, na_idx, n, TT, h_dim, 
                       indicator_config, max_factor_states, varQ)
  factor_hat <- kf_hat$filter_mat[, 1, drop = FALSE]
  
  # Auxiliary indexes
  start_lam <- 1
  end_lam   <- n
  start_phi <- end_lam + 1
  end_phi   <- end_lam + 2
  start_psi <- end_phi + 1
  end_psi   <- end_phi + 2*n
  start_sig <- end_psi + 1
  end_sig   <- end_psi + n 
  
  # Returning everything as a named list
  list(
    convergence = opt$convergence,
    loglik      = -opt$value,
    factor      = factor_hat,
    Hmat        = kf_hat$Hmat,
    Fmat        = kf_hat$Fmat,
    Pmat_final  = kf_hat$Pmat_final,
    Qmat        = kf_hat$Qmat,
    filter_mat  = kf_hat$filter_mat,  
    y_means     = y_means,            
    y_sds       = y_sds,              
    lambda      = list(est = th_hat[start_lam:end_lam], 
                       se = se_hat[start_lam:end_lam]),
    phi         = list(est = th_hat[start_phi:end_phi], 
                       se = se_hat[start_phi:end_phi]),
    psi         = list(est = th_hat[start_psi:end_psi], 
                       se = se_hat[start_psi:end_psi]),
    sd          = list(est = th_hat[start_sig:end_sig], 
                       se = se_hat[start_sig:end_sig]),
    th_hat      = th_hat,
    se_hat      = se_hat
  )
}
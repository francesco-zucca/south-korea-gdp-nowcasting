# ==============================================================================
# MATRIX FUNCTION
# ==============================================================================

matrices <- function(th, n, h_dim, indicator_config, max_factor_states, varQ) {
  
  # ----------------------------------------------------------------------------
  # H MATRIX
  # ----------------------------------------------------------------------------
  
  # Initialize the matrix
  Hmat <- matrix(0, nrow = n, ncol = h_dim)
  
  # Index to get the starting positon of the error in the h vector
  current_error_col <- max_factor_states + 1
  
  # Get the information about the indicator behavior and pull their loadings
  for (i in seq_len(n)) {
    cfg <- indicator_config[[i]]
    lambda_i <- th[i] 
    
    # Determine weights based on indicator type
    if (cfg$type == "MoM") {
      # Monthly indicator loads only on one factor with weight 1
      factor_weights <- c(1)
      error_weights  <- c(1)
      # Numer of errors in the h vector
      n_errors <- 2 
      
    } else if (cfg$type == "QoQ") {
      factor_weights <- c(1/3, 2/3, 1, 2/3, 1/3)
      error_weights  <- c(1/3, 2/3, 1, 2/3, 1/3)
      n_errors <- 5 
      
    } else if (cfg$type == "soft") {
      factor_weights <- rep(1, 12)
      error_weights  <- c(1) 
      n_errors <- 2
    }
    
    # FACTOR BLOCK
    
    # Compute the weighted loadings of each indicator
    scaled_loadings <- lambda_i * factor_weights
    # Fill each indicator row with the respecting scaled loadings
    Hmat[i, 1:length(scaled_loadings)] <- scaled_loadings 
    
    # ERROR BLOCK
    
    # Compute the index of the end of the current error column
    error_end_col <- current_error_col + length(error_weights) - 1
    # Fill the H matrix with the error weights
    Hmat[i, current_error_col:error_end_col] <- error_weights
    # Update the index of the current error column
    current_error_col <- current_error_col + n_errors
  }
  
  # ----------------------------------------------------------------------------
  # F MATRIX
  # ----------------------------------------------------------------------------
  
  # Starting the matrix
  Fmat <- matrix(0, nrow = h_dim, ncol = h_dim)
  
  # Calculating offset of parameters in the theta vector
  start_phi <- n + 1         # (position of phi in the vector of parameters)
  start_psi <- start_phi + 2 # (position of psi in the vector of parameters)
  
  # FACTOR BLOCK
  
  # Top row with AR(2) coefficients for the latent factor
  Fmat[1,1] <- th[start_phi]      # phi_1 entry
  Fmat[1,2] <- th[start_phi + 1]  # phi_2 entry 
  
  # Building the diagonal matrix with 1s on the diagonal
  for (r in 2:max_factor_states) {
    Fmat[r, r - 1] <- 1
  }
  
  # ERROR BLOCK
  
  # Reset the index of the current error column
  current_error_col <- max_factor_states + 1
  
  # Extract the indicator information
  for (i in seq_len(n)) {
    cfg <- indicator_config[[i]]
    
    # Determine the number of errors for each indicator's error
    if (cfg$type == "MoM") {
      n_errors <- 2
    } else if (cfg$type == "QoQ") {
      n_errors <- 5
    } else if (cfg$type == "soft") {
      n_errors <- 2
    }
    
    # Top row of this indicator's block: AR(2) dynamics for the error
    Fmat[current_error_col, 
         current_error_col]     <- th[start_psi + 2*(i - 1)]     # Psi_1
    Fmat[current_error_col, 
         current_error_col + 1] <- th[start_psi + 2*(i - 1) + 1] # Psi_2
    
    # Building the diagonal matrix
    for (r in (current_error_col + 1):(current_error_col + n_errors - 1)) {
      Fmat[r, r - 1] <- 1
      
    }
    
    # Update the index of the current error column
    current_error_col <- current_error_col + n_errors
  }
  
  # ----------------------------------------------------------------------------
  # R MATRIX
  # ----------------------------------------------------------------------------
  
  # R is an (n x n) matrix of zeroes
  Rmat <- matrix(0, n, n)
  
  # ----------------------------------------------------------------------------
  # Q MATRIX
  # ----------------------------------------------------------------------------
  
  # Initialize the Q matrix as a matrix of zeroes
  Qmat <- matrix(0, nrow = h_dim, ncol = h_dim)
  
  # Impute the factor error variance
  Qmat[1,1] <- varQ
  
  # Calculate the offset of parameters in the theta vector
  start_sigma <- start_psi + (2*n) 
  
  # Calculating the variance of each error
  sigma2 <- th[(start_sigma):(start_sigma + n - 1)]^2
  
  # Index to get the starting positon of the error in the error vector
  current_error_col <- max_factor_states + 1
  
  # Extract the indicator information
  for (i in seq_len(n)) {
    cfg <- indicator_config[[i]]
    
    # Determine the number of errors for each indicator's error
    if (cfg$type == "MoM") {
      n_errors <- 2
    } else if (cfg$type == "QoQ") {
      n_errors <- 5
    } else if (cfg$type == "soft") {
      n_errors <- 2
    }
    
    # Inject the variable shock in the diagonal
    Qmat[current_error_col, current_error_col] <- sigma2[i]
    
    # Update the index of the current error column
    current_error_col <- current_error_col + n_errors
  }
  
  # ----------------------------------------------------------------------------
  # RETURNING LIST OF MATRICES 
  # ----------------------------------------------------------------------------
  return(list(
    Hmat = Hmat, 
    Fmat = Fmat, 
    Rmat = Rmat, 
    Qmat = Qmat
  ))
}
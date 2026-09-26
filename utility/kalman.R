# ==============================================================================
# KALMAN FILTER FUNCTION
# ==============================================================================

kalman <- function(th, y, na_idx, n, TT, h_dim, 
                   indicator_config, max_factor_states, varQ) {
  
  # Pulling the matrices
  mats      <- matrices(th, n, h_dim, indicator_config, max_factor_states, varQ)
  Hmat_base <- mats$Hmat  
  Fmat      <- mats$Fmat
  Rmat_base <- mats$Rmat
  Qmat      <- mats$Qmat
  
  # Initializing the vectors to store the likelihood and the factor estimate
  like       <- numeric(TT)
  filter_mat <- matrix(0, TT, h_dim)
  
  # Initial state vector (h_{0|0} = 0) and variance (P_{0|0} = I)
  h00        <- matrix(0, h_dim, 1)
  Pmat00     <- diag(h_dim)
  
  # Kalman filter loop
  for (it in seq_len(TT)) {
    
    # Identify which indicators have data for this month
    obs_idx <- na_idx[it, ] == 1
    
    # Pull H and R matrices (the ones that need imputation)
    Hit <- Hmat_base
    Rit <- Rmat_base
    
    # Zero out rows in the H matrix for missing data
    Hit[!obs_idx, ] <- 0
    
    # Put a variance of 1 in the r matrix for missing data
    diag(Rit)[!obs_idx] <- 1
    
    # Step 1: prediction of ht given info at t-1
    h10      <- Fmat %*% h00
    
    # Step 2: prediction of Pt given info at t-1
    Pmat10   <- Fmat %*% Pmat00 %*% t(Fmat) + Qmat
    
    # Step 3: prediction of yt given info at t-1
    y10      <- Hit %*% h10
    
    # Step 4: estimation of the forecast error
    e10      <- matrix(y[it, ], ncol = 1) - y10
    
    # Step 5: estimation of the variance of the forecast error 
    # We also add very small noise to prevent the matrix from being singular
    # (determinant of 0 and hence not invertible)
    Smat10   <- Hit %*% Pmat10 %*% t(Hit) + Rit + diag(1e-6, n)
    
    # Inverting the S matrix
    Smat10_inv <- solve(Smat10)
    
    # Log-Likelihood calculation
    like[it] <- -0.5 * (n * log(2 * pi) + log(det(Smat10)) +
                          drop(t(e10) %*% Smat10_inv %*% e10))
    
    # Step 6: Kalman gain
    Kmat     <- Pmat10 %*% t(Hit) %*% Smat10_inv
    
    # Step 7: updating of factor and variance
    h11      <- h10 + Kmat %*% e10
    Pmat11   <- Pmat10 - Kmat %*% Hit %*% Pmat10
    
    # Iterating
    filter_mat[it, ] <- t(h11)
    h00 <- h11
    Pmat00 <- Pmat11
  }
  
  # Returning loglikelihood and the extracted matrices
  list(neg_loglik = -sum(like), 
       filter_mat = filter_mat, 
       Hmat       = Hmat_base, 
       Fmat       = Fmat,
       Pmat_final = Pmat00,
       Qmat = Qmat)
}
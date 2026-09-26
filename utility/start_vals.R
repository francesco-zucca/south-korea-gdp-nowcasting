# ==============================================================================
# STARTING VALUES FUNCTION
# ==============================================================================

# Setting starting values (0.5 for lambdas, 0.3 for AR(2) coefficients)
make_startval <- function(y, n) {
  lambda_init <- rep(0.5, n)     # lambdas
  ar_f_init   <- rep(0.3, 2)     # factor AR(2)
  ar_e_init   <- rep(0.3, 2*n)   # errors AR(2)
  sd_init     <- apply(y, 2, sd) # standard deviations
  c(lambda_init, ar_f_init, ar_e_init, sd_init)
}
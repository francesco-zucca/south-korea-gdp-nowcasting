# ==============================================================================
# LOADING LIBRARIES
# ==============================================================================

library("tidyverse")
library("ggplot2")

# ==============================================================================
# LOADING THE DATA AND CALLING THE HELPER FUNCTIONS
# ==============================================================================

# Load raw data
data_raw <- read_csv("panel.csv")
# Select the windows in which we do not have too many NAs
df_clean <- data_raw %>%
  filter(date >= as.Date("2000-02-01") & date <= as.Date("2026-06-01"))

# Define the COVID-19 crisis
covid_start <- as.Date("2020-03-01")
covid_end   <- as.Date("2020-12-01")

# Replace all indicator data with NA during this period 
df_clean <- df_clean %>%
mutate(across(-date, ~ ifelse(date >= covid_start & date <= covid_end, NA, .)))

# Drop the date from the indicator matrix
df <- df_clean %>%
  select(-date) %>%
  as.matrix()

# Declare indicator behavior
indicator_config <- list(
  list(name = "gdp", type = "QoQ", lead = 0),
  list(name = "production", type = "MoM", lead = 0),
  list(name = "imports", type = "MoM", lead = 0),
  list(name = "shipments", type = "MoM", lead = 0),
  list(name = "employment", type = "MoM", lead = 0),
  list(name = "services", type = "MoM", lead = 0),
  list(name = "construction", type = "MoM", lead = 0),
  list(name = "memory_exports", type = "MoM", lead = 2),
  list(name = "sentiment", type = "soft", lead = 1),
  list(name = "spread", type = "MoM", lead = 9)
)

# Call helper functions
source("utility/matrices.R")
source("utility/kalman.R")
source("utility/start_vals.R")
source("utility/estimation.R")
source("utility/stepwise.R")

# ==============================================================================
# RUNNING THE ESTIMATION
# ==============================================================================

fit1 <- dfm(data = df, maxit = 2000)

# ==============================================================================
# EXTRACTING THE FITTED GDP AND PERFORMANCE MEASURES
# ==============================================================================

# Identify how many columns belong to the factor
h_dim <- ncol(fit1$filter_mat)
n     <- ncol(df)

# GDP has 5 error states, the other (n-1) variables have 2 error states each
factor_cols <- h_dim - (5 + (n - 1) * 2) 

# Create a factor only H matrix by zeroing out the error columns
Hmat_factor_only <- fit1$Hmat
Hmat_factor_only[, (factor_cols + 1):h_dim] <- 0

# Extract the fitted standardized values using only the common factor
latent <- t(Hmat_factor_only %*% t(fit1$filter_mat))

# Un-standardize the GDP series
gdp_fitted <- latent[, 1] * fit1$y_sds[1] + fit1$y_means[1]
gdp_actual <- df[, 1]

# Performance evaluations (computed only on known rows of GDP)
valid_idx        <- !is.na(gdp_actual)
gdp_actual_clean <- gdp_actual[valid_idx]
gdp_fitted_clean <- gdp_fitted[valid_idx]

# Compute MSE
residuals <- gdp_actual_clean - gdp_fitted_clean
rmse      <- sqrt(mean(residuals^2))

# Compute R^2
ss_tot    <- sum((gdp_actual_clean - mean(gdp_actual_clean))^2)
ss_res    <- sum(residuals^2)
r_squared <- 1 - (ss_res / ss_tot)

# Compute pearson correlation
gdp_correlation <- cor(gdp_actual_clean, gdp_fitted_clean)

# Print the results
cat(sprintf("Factor R-squared: %.3f\n", r_squared))
cat(sprintf("Factor RMSE: %.3f\n", rmse))
cat(sprintf("Factor-Actual Correlation: %.3f\n", gdp_correlation))

# ==============================================================================
# PLOT AND IN-SAMPLE PERFORMANCE 
# ==============================================================================

# Create date vector to match the data
date_clean <- df_clean$date[valid_idx]

# Combine vectors into a tibble
plot_data <- tibble(
  Date   = date_clean,
  Actual = gdp_actual_clean,    
  Fitted = gdp_fitted_clean
) %>%
  # Pivot to long format for ggplot2
  pivot_longer(cols = c(Actual, Fitted), names_to = "Series", values_to = "Value")

# Generate the quarterly plot
ggplot(plot_data, aes(x = Date, y = Value, color = Series, linetype = Series)) +
  
  # Add lines and points
  geom_line(linewidth = 0.75) +
  
  # Setting colors and line types
  scale_color_manual(values = c("Actual" = "dodgerblue3", "Fitted" = "red3")) +
  scale_linetype_manual(values = c("Actual" = "solid", "Fitted" = "solid")) +
  
  # Adding titles and labels
  labs(
    title = "Actual vs. Fitted Quarterly GDP Growth",
    subtitle = sprintf("In-sample performance: R-squared = %.3f | RMSE = %.3f | Correlation = %.3f", 
                       as.numeric(r_squared), as.numeric(rmse), as.numeric(gdp_correlation)),
    x = "", 
    y = "Quarterly GDP growth (%)",
    color = NULL,
    linetype = NULL
  ) +
  
  # Applying a cleaner, professional theme
  theme_minimal(base_family = "sans") + 
  theme(
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle    = element_text(size = 11, hjust = 0.5, color = "gray30"),
    axis.title       = element_text(size = 10, face = "bold"),
    axis.text        = element_text(size = 9),
    panel.grid.minor = element_blank(), 
    legend.key.width = unit(2, "cm")    
  )

# ==============================================================================
# FORECASTING THE NEXT TWO QUARTERS WITH CONFIDENCE INTERVALS
# ==============================================================================

# Set monthly forecast window
h_fc   <- 6

# Pull the matrices
Fmat   <- fit1$Fmat
Qmat   <- fit1$Qmat
h_row  <- fit1$Hmat[1, ]
h_curr <- matrix(fit1$filter_mat[nrow(fit1$filter_mat), ], ncol = 1)
Pmat_final <- fit1$Pmat_final

# Initialize vector
state_fc <- matrix(0, h_fc, ncol(Fmat))
var_std  <- numeric(h_fc)

# Define rows and columns of the matrix
h_row_mat <- matrix(h_row, nrow = 1)
h_col_mat <- matrix(h_row, ncol = 1)

for (s in seq_len(h_fc)) {
  
  # Forecast the state vector
  h_curr <- Fmat %*% h_curr
  
  # Forecast the state variance
  Pmat_final <- Fmat %*% Pmat_final %*% t(Fmat) + Qmat
  
  # Store the state forecast
  state_fc[s, ] <- t(h_curr)
  
  # Calculate the variance of the GDP forecast
  var_std[s] <- drop(h_row_mat %*% Pmat_final %*% h_col_mat)   
}

# Point forecast and SE, destandardized
gdp_forecast <- as.vector(state_fc %*% h_row) * fit1$y_sds[1] + fit1$y_means[1]
gdp_se       <- sqrt(var_std) * fit1$y_sds[1]

# Keep quarter-end months
last_date    <- max(df_clean$date)
future_dates <- seq(last_date, by = "month", length.out = h_fc + 1)[-1]
qend         <- as.integer(format(future_dates, "%m")) %in% c(3, 6, 9, 12)

forecast_tbl <- tibble(
  date  = future_dates[qend],
  fcst  = gdp_forecast[qend],
  lower = gdp_forecast[qend] - 1.96 * gdp_se[qend],
  upper = gdp_forecast[qend] + 1.96 * gdp_se[qend]
)
print(forecast_tbl)

# ==============================================================================
# STEPWISE FUNCTION
# ==============================================================================

# Run the automated selection algorithm
optimal_model <- optimize_dfm_bic(data = df, 
                                  config = indicator_config, maxit = 2000)

# Extract the winning outputs to use in your plotting and forecasting code!
fit1             <- optimal_model$fit
df               <- optimal_model$data
indicator_config <- optimal_model$config

# ==============================================================================
# EXTRACT FITTED GDP AND PLOT IN-SAMPLE PERFORMANCE FOR OPTIMAL MODEL
# ==============================================================================

h_dim <- ncol(fit1$filter_mat)
n     <- ncol(df)

factor_cols <- h_dim - (5 + (n - 1) * 2) 

Hmat_factor_only <- fit1$Hmat
Hmat_factor_only[, (factor_cols + 1):h_dim] <- 0

latent <- t(Hmat_factor_only %*% t(fit1$filter_mat))

gdp_fitted <- latent[, 1] * fit1$y_sds[1] + fit1$y_means[1]
gdp_actual <- df[, 1]

valid_idx        <- !is.na(gdp_actual)
gdp_actual_clean <- gdp_actual[valid_idx]
gdp_fitted_clean <- gdp_fitted[valid_idx]

residuals <- gdp_actual_clean - gdp_fitted_clean
rmse      <- sqrt(mean(residuals^2))

ss_tot    <- sum((gdp_actual_clean - mean(gdp_actual_clean))^2)
ss_res    <- sum(residuals^2)
r_squared <- 1 - (ss_res / ss_tot)

gdp_correlation <- cor(gdp_actual_clean, gdp_fitted_clean)

date_clean <- df_clean$date[valid_idx]

plot_data <- tibble(
  Date   = date_clean,
  Actual = gdp_actual_clean,    
  Fitted = gdp_fitted_clean
) %>%
  pivot_longer(cols = c(Actual, Fitted), names_to = "Series", values_to = "Value")

opt_plot <- ggplot(plot_data, aes(x = Date, y = Value, color = Series, linetype = Series)) +
  geom_line(linewidth = 0.75) +
  scale_color_manual(values = c("Actual" = "dodgerblue3", "Fitted" = "red3")) +
  scale_linetype_manual(values = c("Actual" = "solid", "Fitted" = "solid")) +
  labs(
    title = "Restricted Model: Actual vs. Fitted Quarterly GDP Growth",
    subtitle = sprintf("In-sample performance: R-squared = %.3f | RMSE = %.3f | Correlation = %.3f", 
                       as.numeric(r_squared), as.numeric(rmse), as.numeric(gdp_correlation)),
    x = "", 
    y = "Quarterly GDP growth (%)",
    color = NULL,
    linetype = NULL
  ) +
  theme_minimal(base_family = "sans") + 
  theme(
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle    = element_text(size = 11, hjust = 0.5, color = "gray30"),
    axis.title       = element_text(size = 10, face = "bold"),
    axis.text        = element_text(size = 9),
    panel.grid.minor = element_blank(), 
    legend.key.width = unit(2, "cm")    
  )

print(opt_plot)

# ==============================================================================
# FORECAST THE NEXT TWO QUARTERS USING THE RESTRICTED MODEL
# ==============================================================================

h_fc   <- 6
Fmat   <- fit1$Fmat
Qmat   <- fit1$Qmat
h_row  <- fit1$Hmat[1, ]
h_curr <- matrix(fit1$filter_mat[nrow(fit1$filter_mat), ], ncol = 1)
Pmat_final <- fit1$Pmat_final

state_fc <- matrix(0, h_fc, ncol(Fmat))
var_std  <- numeric(h_fc)

h_row_mat <- matrix(h_row, nrow = 1)
h_col_mat <- matrix(h_row, ncol = 1)

for (s in seq_len(h_fc)) {
  h_curr <- Fmat %*% h_curr
  Pmat_final <- Fmat %*% Pmat_final %*% t(Fmat) + Qmat
  state_fc[s, ] <- t(h_curr)
  var_std[s] <- drop(h_row_mat %*% Pmat_final %*% h_col_mat)   
}

gdp_forecast <- as.vector(state_fc %*% h_row) * fit1$y_sds[1] + fit1$y_means[1]
gdp_se       <- sqrt(var_std) * fit1$y_sds[1]

last_date    <- max(df_clean$date)
future_dates <- seq(last_date, by = "month", length.out = h_fc + 1)[-1]
qend         <- as.integer(format(future_dates, "%m")) %in% c(3, 6, 9, 12)

optimal_forecast_tbl <- tibble(
  date  = future_dates[qend],
  fcst  = gdp_forecast[qend],
  lower = gdp_forecast[qend] - 1.96 * gdp_se[qend],
  upper = gdp_forecast[qend] + 1.96 * gdp_se[qend]
)

print(optimal_forecast_tbl)
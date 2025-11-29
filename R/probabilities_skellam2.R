# Fit each match based on skellam2
match_probabilities_skellam2 <- function(Z_range, mu, sigma){

  probs = dskellam2(Z_range, mu, sigma, Ind = FALSE)

  P_home <- sum(probs[Z_range > 0 ])
  P_draw <- probs[Z_range == 0]
  P_away <- sum(probs[Z_range < 0 ])

  total <- P_home + P_draw + P_away

  return(c(
    home_win = P_home / total,
    draw     = P_draw / total,
    away_win = P_away / total
  ))

}

# By using the function fit the league
league_probabilities_skellam2 <- function(z_df, model){

  # Model matrix
  mm = model$model_matrix
  # model parameters (Team's home/away abillity)
  param = model$model$par[1:(ncol(mm))]
  # sigma^2:
  sigma_hat <- tail(model$model$par, 1)

  # Fitted values (mu's for the whole season)
  skellam2_fitted <- mm %*% param
  # Poisson lambda:
  lambda_h_hat <- 0.5 * (sigma_hat + skellam2_fitted)
  lambda_a_hat <- 0.5 * (sigma_hat - skellam2_fitted)

  # Actual outcomes:
  Z <- z_df$Z
  # Proportions of fitted values:
  fitted_probs <- t(sapply(
    1:nrow(mm),
    function(i) match_probabilities_skellam2(
      min(Z):max(Z), skellam2_fitted[i], sigma_hat)
  ))

  cbind(z_df, fitted_probs, lambda_h = lambda_h_hat, lambda_a = lambda_a_hat)
}

# Fit each match based on skellam2
match_probabilities_ZIskellam2 <- function(Z_range, mu, sigma, p){
  #Z_range = -6:6 ; mu = 1.2 ;  sigma = 2 ; p = 0.04

  probs = ifelse(
    Z_range == 0, (p + ((1 - p) * dskellam2(0, mu, sigma, Ind = FALSE))),
    (1 - p) * dskellam2(Z_range, mu, sigma, Ind = FALSE)
  )

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
league_probabilities_ZIskellam2 <- function(z_df, model){

  # Model matrix
  mm = model$model_matrix
  # model parameters (Team's home/away abillity)
  param = model$model$par[1:(ncol(mm))]
  # sigma^2:
  sigma_hat <- tail(model$model$par, 2)[1]
  # ZI p
  p_hat <- plogis(tail(model$model$par, 1))

  # Fitted values (mu's for the whole season)
  ZIskellam2_fitted <- mm %*% param
  #Poisson lambda:
  lambda_h_hat <- 0.5 * (sigma_hat + ZIskellam2_fitted)
  lambda_a_hat <- 0.5 * (sigma_hat - ZIskellam2_fitted)

  # Actual outcomes:
  Z <- z_df$Z
  # Proportions of fitted values:
  fitted_probs <- t(sapply(
    1:nrow(mm),
    function(i) match_probabilities_ZIskellam2(
      min(Z):max(Z), ZIskellam2_fitted[i], sigma_hat, p_hat
    )
  ))

  cbind(z_df, fitted_probs, lambda_h = lambda_h_hat, lambda_a = lambda_a_hat)
}




match_probabilities_Z2skellam2 <- function(
    Z_range, mu1, mu2, sigma1, sigma2){

  n <- length(Z_range)

  Z1_grid <- rep(Z_range, each = n)
  Z2_grid <- rep(Z_range, times = n)

  probs <- dskellam2(Z1_grid, mu1, sigma1, Ind = FALSE) *
    dskellam2(Z2_grid, mu2, sigma2, Ind = FALSE)

  Z_fitted = Z1_grid + Z2_grid

  P_home <- sum(probs[Z_fitted > 0 ])
  P_draw <- sum(probs[Z_fitted == 0])
  P_away <- sum(probs[Z_fitted < 0 ])

  total <- P_home + P_draw + P_away

  return(c(
    home_win = P_home / total,
    draw     = P_draw / total,
    away_win = P_away / total
  ))

}

league_probabilities_Z2skellam2 <- function(z_df, model){

  # Model matrix
  mm = model$model_matrix
  # model parameters (Team's home/away abillity)
  param = model$model$par[1:(ncol(mm))]

  # sigma^2:
  sigma1_hat <- tail(model$model$par, 2)[1]
  sigma2_hat <- tail(model$model$par, 1)

  mus <- mm %*% param

  lambda1_h_hat <- 0.5 * (sigma1_hat + mus)
  lambda1_a_hat <- 0.5 * (sigma1_hat - mus)

  lambda2_h_hat <- 0.5 * (sigma2_hat + mus)
  lambda2_a_hat <- 0.5 * (sigma2_hat - mus)

  Z1_fitted <- Z2_fitted <- mus

  Z_range <- min(z_df$Z) : max(z_df$Z)
  # Proportions of fitted values:
  fitted_probs <- t(sapply(
    1:nrow(mm),
    function(i) match_probabilities_Z2skellam2(
      Z_range, mus[i], mus[i], sigma1_hat, sigma2_hat
    )
  ))

  cbind(z_df, fitted_probs, Z1_fitted = Z1_fitted, Z2_fitted = Z2_fitted,
        l1_h = lambda1_h_hat, l1_a = lambda1_a_hat,
        l2_h = lambda2_h_hat, l2_a = lambda2_a_hat)
}

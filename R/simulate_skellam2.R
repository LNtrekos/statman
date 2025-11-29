# Simulate match:
simulate_match_skellam2 <- function(Z_range, n = 1, mu, sigma){

  match_probs = dskellam2(Z_range, mu, sigma, Ind = FALSE)
  match_result = sample(Z_range, n, prob = match_probs)
  return(match_result)

}

simulate_season_skellam2 <- function(B = 10^4, z_df, model){

  # --- Precompute everything --- #
  mm <- model$model_matrix
  param <- model$model$par[1:ncol(mm)]
  sigma_hat <- tail(model$model$par, 1)

  skellam2_fitted <- mm %*% param
  lambda_h_hat <- 0.5 * (sigma_hat + skellam2_fitted)
  lambda_a_hat <- 0.5 * (sigma_hat - skellam2_fitted)

  Z_range <- min(z_df$Z):max(z_df$Z)
  nteams <- length(levels(z_df$Home))
  teams <- levels(z_df$Home)

  # --- Results array --- #
  metrics <- c("P", "R", "GF", "GA", "GD")
  sim_array <- array(
    0, dim = c(B, nteams, length(metrics)),
    dimnames = list(NULL, teams, metrics)
  )

  # --- Monte-Carlo simulation --- #
  for (b in 1:B) {

    if (b %% 1000 == 0){
      print(b)
    }
    Z_sim <- sapply(
      1:nrow(mm),
      function(i) simulate_match_skellam2(
        Z_range, 1, skellam2_fitted[i], sigma_hat
      )
    )

    HG_sim <- rpois(nrow(mm), lambda_h_hat)
    AG_sim <- rpois(nrow(mm), lambda_a_hat)

    res <- univariate_league(Z_sim, HG_sim, AG_sim, z_df$Home, z_df$Away)

    sim_array[b, , "P"]  <- res$points
    sim_array[b, , "R"]  <- res$pos
    sim_array[b, , "GF"] <- res$GF
    sim_array[b, , "GA"] <- res$GA
    sim_array[b, , "GD"] <- res$GD
  }
  return(sim_array)
}




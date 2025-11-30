# Fit each match based on skellam2
simulate_match_Z2skellam2 <- function(
    Z_range, n = 1, mu1, mu2, sigma1, sigma2){

  n_z <- length(Z_range)

  Z1_grid <- rep(Z_range, each = n_z)
  Z2_grid <- rep(Z_range, times = n_z)

  probs <- dskellam2(Z1_grid, mu1, sigma1, Ind = FALSE) *
    dskellam2(Z2_grid, mu2, sigma2, Ind = FALSE)

  ind = sample(1:length(probs), n, prob = probs)

  return(c(Z1_grid[ind], Z2_grid[ind]))
}

simulate_match_Z2skellam2_consistent <- function(
    lambda1_h, lambda1_a,
    lambda2_h, lambda2_a) {

  HG1 <- rpois(1, lambda1_h)
  AG1 <- rpois(1, lambda1_a)
  Z1  <- HG1 - AG1

  HG2 <- rpois(1, lambda2_h)
  AG2 <- rpois(1, lambda2_a)
  Z2  <- HG2 - AG2

  Z <- Z1 + Z2

  return(list(Z=Z, Z1=Z1, Z2=Z2,
              HG1=HG1, AG1=AG1,
              HG2=HG2, AG2=AG2))
}

# By using the function fit the league
simulate_season_Z2skellam2 <- function(B = 10^4, z_df, model){

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

  Z_range <- min(z_df$Z):max(z_df$Z)
  nteams <- length(levels(z_df$Home))
  teams <- levels(z_df$Home)

  # --- Results array --- #
  metrics <- c(
    "P", "R", "GF", "GA", "GD"
    #"GF1", "GF2", "GA1", "GA2", "GD1", "GD2"
  )
  sim_array <- array(
    0, dim = c(B, nteams, length(metrics)),
    dimnames = list(NULL, teams, metrics)
  )

  # --- Monte-Carlo simulation --- #
  for (b in 1:B) {
    #print(b)
    if (b %% 1000 == 0) message("Simulation: ", b)

    Z_sims <- t(sapply(
      1:nrow(mm),
      function(i) simulate_match_Z2skellam2(
        Z_range, n = 1, mus[i], mus[i], sigma1_hat, sigma2_hat
      )
    ))

    Z_sim = Z_sims[,1] + Z_sims[,2]

    HG1_sim <- rpois(nrow(mm), lambda1_h_hat)
    AG1_sim <- rpois(nrow(mm), lambda1_a_hat)
    HG2_sim <- rpois(nrow(mm), lambda2_h_hat)
    AG2_sim <- rpois(nrow(mm), lambda2_a_hat)

    res <- bivariate_league(
      Z_sim, HG1_sim, AG1_sim, HG2_sim, AG2_sim,
      z_df$Home, z_df$Away
    )

    sim_array[b, , "P"]  <- res$points
    sim_array[b, , "R"]  <- res$pos
    sim_array[b, , "GF"] <- res$GF
    sim_array[b, , "GA"] <- res$GA
    sim_array[b, , "GD"] <- res$GD
  }

  return(sim_array)
}


simulate_season_Z2skellam2_consistent <- function(B = 10^4, z_df, model){

  # Model matrix
  mm <- model$model_matrix
  param <- model$model$par[1:ncol(mm)]

  sigma1_hat <- tail(model$model$par, 2)[1]
  sigma2_hat <- tail(model$model$par, 1)

  mus <- mm %*% param

  lambda1_h_hat <- 0.5 * (sigma1_hat + mus)
  lambda1_a_hat <- 0.5 * (sigma1_hat - mus)

  lambda2_h_hat <- 0.5 * (sigma2_hat + mus)
  lambda2_a_hat <- 0.5 * (sigma2_hat - mus)

  nteams <- length(levels(z_df$Home))
  teams <- levels(z_df$Home)
  n_matches <- nrow(mm)

  metrics <- c("P", "R", "GF", "GA", "GD")
  sim_array <- array(
    0, dim = c(B, nteams, length(metrics)),
    dimnames = list(NULL, teams, metrics)
  )

  for (b in 1:B) {
    if (b %% 1000 == 0) message("Simulation: ", b)

    # Allocate match-level vectors
    Z_sim  <- numeric(n_matches)
    Z1_sim <- numeric(n_matches)
    Z2_sim <- numeric(n_matches)
    HG1_sim <- numeric(n_matches)
    AG1_sim <- numeric(n_matches)
    HG2_sim <- numeric(n_matches)
    AG2_sim <- numeric(n_matches)

    # Simulate each match
    for (i in 1:n_matches) {
      sim <- simulate_match_Z2skellam2_consistent(
        lambda1_h_hat[i], lambda1_a_hat[i],
        lambda2_h_hat[i], lambda2_a_hat[i]
      )

      Z_sim[i]  <- sim$Z
      Z1_sim[i] <- sim$Z1
      Z2_sim[i] <- sim$Z2
      HG1_sim[i] <- sim$HG1
      AG1_sim[i] <- sim$AG1
      HG2_sim[i] <- sim$HG2
      AG2_sim[i] <- sim$AG2
    }

    # Evaluate league table
    res <- bivariate_league(
      Z_sim, HG1_sim, AG1_sim, HG2_sim, AG2_sim,
      z_df$Home, z_df$Away
    )

    sim_array[b, , "P"]  <- res$points
    sim_array[b, , "R"]  <- res$pos
    sim_array[b, , "GF"] <- res$GF
    sim_array[b, , "GA"] <- res$GA
    sim_array[b, , "GD"] <- res$GD
  }

  return(sim_array)
}






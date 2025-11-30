Z2skellam2_ll <- function(theta, Z1, Z2, X) {

  n_beta <- ncol(X)
  beta <- theta[1:n_beta]
  sigma1 <- theta[n_beta + 1]
  sigma2 <- theta[n_beta + 2]

  eta <- X %*% beta

  # Identifiability: sigma_j > |eta|
  if (any(sigma1 <= abs(eta)) || any(sigma2 <= abs(eta))) {
    return(1e10)
  }
  li <- dskellam2(Z1, eta, sigma1, Ind = TRUE) +
    dskellam2(Z2, eta, sigma2, Ind = TRUE)
  ll <- sum(li)

  return(-ll)
}

# Estimate for skellam2
fit_Z2skellam2 <- function(z_df){

  # z_df must contain Z: the difference
  model_matrix <- model.matrix(
    ~ Home + Away, data = z_df
  )

  # starting values
  initial_beta <- coef(lm(Z ~ Home + Away, data = z_df)) + rnorm(1, 0, 0.1)
  initial_theta <- c(initial_beta, sigma1 = 3, sigma2 = 3)

  model <- optim(
    initial_theta, Z2skellam2_ll,
    Z1 = z_df$Z1, Z2 = z_df$Z2, X = model_matrix,
    control = list(maxit=50000), hessian=TRUE
  )
  model$par
  return(list(
    model = model,
    model_matrix = model_matrix
  ))

}

summarize_Z2skellam2 <- function(z_df = z_df, model){

  teams = levels(z_df$Home)
  nteams = length(teams)

  ## --- Safe Hessian inversion with fallback ---
  SE <- NA   # default if Hessian fails

  SE <- tryCatch({
    H <- model$model$hessian
    Vcov <- solve(H)                   # <- may fail
    sqrt(diag(Vcov))
  },
  error = function(e){
    message("Warning: Hessian inversion failed. Standard errors set to NA.")
    return(rep(NA, length(model$model$par)))
  })

  ## --- Build results ---
  results <- list(

    parameters = data.frame(
      mean = c(model$model$par[1], tail(model$model$par, 2))
      # ,sd = c(SE[1], tail(SE, 2))     # uncomment if you want
    ),

    team_abillities = data.frame(
      home = model$model$par[2:nteams],
      #home_sd = SE[2:nteams],
      away = model$model$par[(nteams + 1):(2*nteams - 1)]
      #,away_sd = SE[(nteams + 1):(2*nteams - 1)]
    ),

    value = model$model$value
  )

  rownames(results$team_abillities) <- levels(z_df$Home)[-1]

  return(results)
}


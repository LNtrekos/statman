#Skellam2 loglikelihood with covariates
skellam2_ll <- function(theta, y, X) {

  n_beta <- ncol(X)
  beta <- theta[1:n_beta]
  s <- theta[n_beta + 1]

  m <- X %*% beta

  # Enforce constraints: lam1 > 0, lam2 > 0
  if(any(m + s <= 0) || any(s - m <= 0)) return(1e10)

  ll <- sum(dskellam2(y, m, s))

  return(-ll)
}

# Estimate for skellam2
fit_skellam2 <- function(z_df){

  # z_df must contain Z: the difference
  model_matrix <- model.matrix(
    ~ Home + Away, data = z_df
  )

  # starting values
  initial_beta <- coef(lm(Z ~ Home + Away, data = z_df)) + rnorm(1, 0, 0.1)
  initial_theta <- c(initial_beta, s = 5)

  model <- optim(
    initial_theta, skellam2_ll,
    y = z_df$Z, X = model_matrix,
    control = list(maxit=50000), hessian=TRUE
  )

  return(list(
    model = model,
    model_matrix = model_matrix
  ))

}

summarize_skellam2 <- function(model, z_df){

  teams = levels(z_df$Home)
  nteams = length(teams)

  alpha_hat <- model$model$par[1]
  beta_hat <- model$model$par[2:nteams]
  gamma_hat <- model$model$par[(nteams + 1):(2*nteams - 1)]
  sigma2_hat <- model$model$par[(2*nteams)]

  skellam2_model_results <- list(
    alpha_hat = alpha_hat,

    team_abillities = data.frame(
      betas = beta_hat,
      #betas_sd = betas_sd,
      gammas = gamma_hat,
      #gamma_sd = gamma_sd,
      home_strength = beta_hat + gamma_hat
    ),

    variance = sigma2_hat,

    value = model$model$value
  )

  rownames(skellam2_model_results$team_abillities) = teams[-1]
    skellam2_model_results$team_abillities = round(
      skellam2_model_results$team_abillities, 3
  )

  return(skellam2_model_results)
}

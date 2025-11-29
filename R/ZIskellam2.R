# ZI_Skellam2 loglikelihood with covariates
ZIskellam2_ll <- function(theta, y, X){

  n_beta <- ncol(X)
  beta <- theta[1:n_beta]
  s <- theta[n_beta + 1]

  logit_p <- theta[n_beta + 2]
  p <- 1 / (1 + exp(-logit_p))

  m <- X %*% beta

  if(any(m + s <= 0) || any(s - m <= 0)) {
    return(1e10)  # Return large value for invalid parameters
  }

  f_y <- dskellam2(y, m, s, Ind = F)
  f0 <- dskellam2(0, m, s, Ind = F)

  # log-likelihood vector
  log_lik_i <- ifelse(y == 0, log(p + (1 - p) * f0), log((1 - p) * f_y))
  log_lik <- -sum(log_lik_i)

  return(log_lik)

}

# Estimate for ZI_skellam2
fit_ZIskellam2 <- function(z_df){

  # z_df must contain Z: the difference
  model_matrix <- model.matrix(
    ~ Home + Away, data = z_df
  )

  # starting values
  initial_beta <- coef(lm(Z ~ Home + Away, data = z_df)) + rnorm(1, 0, 0.1)
  initial_theta <- c(initial_beta, s = 5, logit_p = 0)

  model <- optim(
    initial_theta, ZIskellam2_ll,
    y = z_df$Z, X = model_matrix,
    control = list(maxit=50000), hessian=TRUE
  )

  return(list(
    model = model,
    model_matrix = model_matrix
  ))

}

# Data frame for readable results
summarize_ZIskellam2 <- function(model, z_df){

  teams = levels(z_df$Home)
  nteams = length(teams)

  model = model$model

  alpha_hat <- model$par[1]
  beta_hat <- model$par[2:nteams]
  gamma_hat <- model$par[(nteams + 1):(2*nteams - 1)]
  sigma2_hat <- model$par[(2*nteams)]
  logit_p_hat <- tail(model$par, 1)

  model_results <- list(
    alpha_hat = alpha_hat,

    team_abillities = data.frame(
      betas = beta_hat,
      #betas_sd = betas_sd,
      gammas = gamma_hat,
      #gamma_sd = gamma_sd,
      home_strength = beta_hat + gamma_hat
    ),

    variance = sigma2_hat,

    p = plogis(logit_p_hat),

    value = model$value
  )

  rownames(model_results$team_abillities) = teams[-1]
  model_results$team_abillities = round(model_results$team_abillities, 3)

  return(model_results)

}

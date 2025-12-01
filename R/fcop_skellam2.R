frank_pcop <- function(u, v, theta) {

  # Near 0 → independence (avoid numerical problems)
  if (abs(theta) < 1e-6) {
    return(u * v)
  }

  num <- (exp(-theta * u) - 1) * (exp(-theta * v) - 1)
  den <- exp(-theta) - 1

  -(1 / theta) * log1p(num / den)   # log1p(x) = log(1 + x), more stable
}

## ===============
## Gumbel copula
## ===============
gumbel_cop <- function(u, v, theta) {

  #if (theta < 1) stop("Gumbel parameter theta must be >= 1")

  # near independence
  if (abs(theta - 1) < 1e-6) {
    return(u * v)
  }

  a <- (-log(u))^theta
  b <- (-log(v))^theta

  exp(-(a + b)^(1/theta))
}

# Skellam2 density function:
dskellam2 <- function(x, m, s, Ind = TRUE) {
  lam1 <- 0.5 * (m + s)
  lam2 <- s - lam1
  dskellam(x, lam1, lam2, log = Ind)
}

dmass_skellam2_cop <- function(
    Z1, Z2, mu1, mu2, sigma1, sigma2, pcop, theta){

  U1  <- pskellam2(Z1,     mu1, sigma1)
  U1m <- pskellam2(Z1 - 1, mu1, sigma1)

  V1  <- pskellam2(Z2,     mu2, sigma2)
  V1m <- pskellam2(Z2 - 1, mu2, sigma2)

  prop <-
    pcop(U1,  V1,  theta) -
    pcop(U1m, V1,  theta) -
    pcop(U1,  V1m, theta) +
    pcop(U1m, V1m, theta)

  prop[prop <= 0] <- 1e-16
  prop
}


cop_skellam2_ll <- function(param, Z1, Z2, modelmatrix, pcop) {

  # param: c(beta (k), sigma1, sigma2, cpar)
  k <- ncol(modelmatrix)

  beta <- param[1:k]
  sigma1 <- param[k + 1]
  sigma2 <- param[k + 2]
  theta <- param[k + 3]

  # Linear predictor for the Skellam mean (common to both halves)
  eta <- as.vector(modelmatrix %*% beta)

  # Identifiability: sigma_j > |eta|
  if (any(sigma1 <= abs(eta)) || any(sigma2 <= abs(eta))) {
    return(1e10)
  }

  # Compute joint pmf
  loglik_i <- dmass_skellam2_cop(Z1, Z2, eta, eta, sigma1, sigma2, pcop, theta)

  if (any(!is.finite(loglik_i))) {
    return(1e10)
  }

  loglik <- -sum(log(loglik_i))

  return(loglik)
}

fit_cop_skellam2 <- function(dataframe, pcop, initial_theta = 0){

  model_matrix <- model.matrix(~ Home + Away, data = dataframe)

  ## Starting values
  initial_beta <- coef(lm(Z1 ~ Home + Away, data = dataframe))

  initial_sigma1 <- initial_sigma2 <- 3

  initial_params <- c(
    initial_beta, sigma1 = initial_sigma1, sigma2 = initial_sigma2,
    theta = initial_theta
  )

  model <- optim(
    initial_params, fn = cop_skellam2_ll, Z1 = dataframe$Z1, Z2 = dataframe$Z2,
    modelmatrix = model_matrix, pcop = pcop,
    method = "BFGS", control = list(maxit = 5000, reltol = 1e-6),
    hessian = TRUE
  )

  return(list(
    model = model,
    model_matrix = model_matrix
  ))
}

# Data frame for readable results
summarize_cop_skellam2 <- function(dataframe, model){

  teams = levels(dataframe$Home)
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
      mean = c(model$model$par[1], tail(model$model$par, 3))
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
  rownames(results$team_abillities) <- levels(dataframe$Home)[-1]

  return(results)
}

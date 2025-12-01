devtools::load_all()

data_path <- "C:/Users/lntre/Desktop/thesis_data"
modelling_df_path <- file.path(data_path, "modelling")
overview_path <- file.path(data_path, "overview")

cop_skellam2_path <- "C:/Users/lntre/Desktop/results/cop_skellam2"

raw_leagues <- c("Premier League", "Serie A", "La Liga", "Bundesliga", "Ligue 1")
leagues <- c("pl_", "sa_", "ll_", "bl_", "l1_")
seasons <- c("2020-2021","2021-2022","2022-2023","2023-2024","2024-2025")

std_prfx <- "modelling_df_"

for (j in seq_along(leagues)){

  for (i in seq_along(seasons)){
    j = 1 ; i = 3
    message("Currently in league : ", raw_leagues[j], " and Season: ", seasons[i])

    modelling_file_path <- file.path(
      modelling_df_path, paste0(std_prfx, leagues[j], seasons[i], ".csv")
    )

    z_df <- read.csv(modelling_file_path) %>%
      mutate(Home = as.factor(Home), Away = as.factor(Away)) %>%
      select(Home, Away, Z1, Z2, Z, H1, A1, H2, A2)


    t0 = Sys.time()
    # fitted the non inflated model
    cop_skellam2_model <- fit_cop_skellam2(z_df, frank_pcop, initial_theta = 0)
    fitting_time = Sys.time() - t0

    message("Fitting time:", fitting_time)

    # Check the results:
    cop_skellam2_results <- summarize_cop_skellam2(
      model = cop_skellam2_model, z_df
    )

    fitted_league_cop_skellam2 <- fit_league_cop_skellam2(
      z_df, cop_skellam2_model)

    apply(fitted_league_cop_skellam2[,10:12], 2, sum)
    c(sum(z_df$Z > 0), sum(z_df$Z == 0), sum(z_df$Z < 0 ))

    B <- 1000
    teams <- levels(z_df$Home)
    nteams <- length(teams)

    message("Simulation:")
    t0 = Sys.time()

    ## --- Model pieces from cop_skellam2_model --- ##
    mm <- cop_skellam2_model$model_matrix
    parvec <- cop_skellam2_model$model$par

    param <- parvec[1:ncol(mm)]
    sigma1_hat <- tail(parvec, 3)[1]
    sigma2_hat <- tail(parvec, 3)[2]
    theta_hat <- tail(parvec, 1)

    mus <- as.vector(mm %*% param)

    lambda1_h_hat <- 0.5 * (sigma1_hat + mus)
    lambda1_a_hat <- 0.5 * (sigma1_hat - mus)

    lambda2_h_hat <- 0.5 * (sigma2_hat + mus)
    lambda2_a_hat <- 0.5 * (sigma2_hat - mus)

    # Goal difference range
    Z_range <- min(z_df$Z):max(z_df$Z)

    # Teams & matches
    teams <- levels(z_df$Home)
    nteams <- length(teams)
    nmatches <- nrow(mm)

    home_i <- as.integer(z_df$Home)
    away_i <- as.integer(z_df$Away)

    n_z <- length(Z_range)
    Z1_grid <- rep(Z_range, each = n_z)
    Z2_grid <- rep(Z_range, times = n_z)

    grid_list <- lapply(seq_len(nmatches), function(i) {

      probs <- dmass_skellam2_fcop_cpp(
        Z1_grid, Z2_grid,
        mu1 = mus[i], mu2 = mus[i],
        sigma1 = sigma1_hat, sigma2 = sigma2_hat,
        theta  = theta_hat, pskellam2 = pskellam2
      )

      cdf <- cumsum(probs)

      list(cdf = cdf)
    })

    sim_array <- simulate_seasons_cop_skellam2(
      B = B, mus = mus,
      lambda1_h_hat = lambda1_h_hat,
      lambda1_a_hat = lambda1_a_hat,
      lambda2_h_hat = lambda2_h_hat,
      lambda2_a_hat = lambda2_a_hat,
      Z1_grid  = Z1_grid, Z2_grid  = Z2_grid,
      grid_list = grid_list,
      home_i = home_i, away_i = away_i,
      teams  = teams
    )
    sim_time = Sys.time() - t0
    message("Sim time:", sim_time)

    summary_tab <- summarize_simulations(sim_array)
    idx <- order(summary_tab[, 2], decreasing = TRUE)

    rank_tab <- ranking_probability_table(sim_array)

    # Create folder safely
    folder_path <- file.path(
      cop_skellam2_path, paste0(leagues[j], seasons[i])
    )
    if (!dir.exists(folder_path)) dir.create(folder_path, recursive = TRUE)

    # Heatmap
    heatmap <- plot_heatmap(rank_tab[idx,])
    ggsave(
      file.path(folder_path, "ranking.pdf"),
      plot = heatmap,
      width = 10, height = 5.15, units = "in"
    )

    # Overall table
    overall_df <- read.csv(
      file.path(
        overview_path, paste0(leagues[j], seasons[i], ".csv")
      )
    )

    # Points plot
    points_with_errobars <- plot_points_with_errobars(
      summary_tab, overall_df
    )
    ggsave(
      file.path(folder_path, "points.pdf"),
      plot = points_with_errobars,
      width = 10, height = 5.15, units = "in"
    )

    # GF / GA plot
    goals_plot <- plot_gf_ga(sim_array, overall_df)
    ggsave(
      file.path(folder_path, "goals.pdf"),
      plot = goals_plot,
      width = 10, height = 5.15, units = "in"
    )

    ## ======================
    ##  Error metrics + save
    ## ======================

    # Expected vs actual points
    pred_df <- as.data.frame(summary_tab[, c("Team", "Exp-P")])
    colnames(pred_df) <- c("Squad", "Exp_P")
    pred_df$Exp_P <- as.numeric(pred_df$Exp_P)

    mse_df <- merge(overall_df, pred_df, by = "Squad")

    points_mse <- mean((mse_df$Exp_P - mse_df$Pts)^2)
    points_ae  <- mean(abs(mse_df$Exp_P - mse_df$Pts))

    # GF MSE / AE
    gf_tab <- gf_df(sim_array, overall_df)
    gf_mse <- mean((gf_tab[, "Actual_GF"] - gf_tab[, "50%"])^2)
    gf_ae  <- mean(abs(gf_tab[, "Actual_GF"] - gf_tab[, "50%"]))

    # GA MSE / AE
    ga_tab <- ga_df(sim_array, overall_df)
    ga_mse <- mean((ga_tab[, "Actual_GA"] - ga_tab[, "50%"])^2)
    ga_ae  <- mean(abs(ga_tab[, "Actual_GA"] - ga_tab[, "50%"]))

    # Final results (BUG FIXED here: correct variables used)
    final_results <- list(
      Value         = cop_skellam2_model$model$value,
      special_param = tail(cop_skellam2_model$model$par, 3),
      Est_probs     = colSums(fitted_league_cop_skellam2[, 10:12]),
      MSE = cbind(Point = points_mse,
                  GF    = gf_mse,
                  GA    = ga_mse),
      AE  = cbind(Point = points_ae,
                  GF    = gf_ae,
                  GA    = ga_ae)
    )

    # Save & read back
    rds_name <- paste0(leagues[j], seasons[i], ".rds")
    saveRDS(
      final_results,
      file = file.path(folder_path, rds_name)
    )

    obj <- readRDS(file.path(folder_path, rds_name))
  }
}

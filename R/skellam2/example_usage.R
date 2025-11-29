devtools::load_all()

data_path <- "C:/Users/lntre/Desktop/thesis_data"
modelling_df_path <- file.path(data_path, "modelling")
overview_path <- file.path(data_path, "overview")

skellam2_path <- "C:/Users/lntre/Desktop/results/skellam2"

raw_leagues <- c("Premier League", "Serie A")
leagues <- c("pl_", "sa_")
seasons <- c("2020-2021","2021-2022","2022-2023","2023-2024","2024-2025")

std_prfx <- "modelling_df_"
j = 2 ; i = 1

for (j in seq_along(leagues)){

  for (i in seq_along(seasons)){

    message("Currently in league : ", raw_leagues[j], " and Season: ", seasons[i])

    modelling_file_path <- file.path(
      modelling_df_path, paste0(std_prfx, leagues[j], seasons[i], ".csv")
    )

    z_df <- read.csv(modelling_file_path) %>%
      mutate(Home = as.factor(Home), Away = as.factor(Away)) %>%
      select(Home, Away, Z)

    t0 = Sys.time()
    # fitted the non inflated model
    skellam2_model <- fit_skellam2(z_df = z_df)
    fitting_time = Sys.time() - t0

    message("Fitting time:", fitting_time)

    # Check the results:
    skellam2_results <- summarize_skellam2(model = skellam2_model, z_df)

    skellam2_fitted_league <- league_probabilities_skellam2(z_df, skellam2_model)

    B = 10
    teams <- levels(z_df$Home)
    nteams <- length(teams)

    message("Simulation:")
    t0 = Sys.time()
    sim_array <- simulate_season_skellam2(B, z_df = z_df, skellam2_model)
    sim_time = Sys.time() - t0
    message("Sim time:", sim_time)

    summary_tab <- summarize_simulations(sim_array)
    idx <- order(summary_tab[, 2], decreasing = TRUE)

    rank_tab <- ranking_probability_table(sim_array)

    # Create folder safely
    folder_path <- file.path(skellam2_path, paste0(leagues[j], seasons[i]))
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
      Value         = skellam2_model$model$value,
      special_param = tail(skellam2_model$model$par, 1),
      Est_probs     = colMeans(skellam2_fitted_league[, 4:6]),
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

    #obj <- readRDS(file.path(folder_path, rds_name))
  }
2}

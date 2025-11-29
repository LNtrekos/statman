# --- League function --- #
univariate_league <- function(Z, HG, AG, HT, AT) {

  hp <- 3*(Z > 0) + (Z == 0)
  ap <- 3*(Z < 0) + (Z == 0)

  points <- tapply(hp, HT, sum) + tapply(ap, AT, sum)
  GF <- tapply(HG, HT, sum) + tapply(AG, AT, sum)
  GA <- tapply(AG, HT, sum) + tapply(HG, AT, sum)

  pos <- nteams + 1 - rank(points, ties.method="max")

  list(points=points, pos=pos, GF=GF, GA=GA, GD=GF-GA)
}


# Summary table:
summarize_simulations <- function(sim_array){

  P <- sim_array[, , "P"]
  R <- sim_array[, , "R"]

  teams  <- dimnames(sim_array)[[2]]
  nteams <- length(teams)

  mean_P <- colMeans(P)
  sd_P   <- apply(P, 2, sd)

  q_P <- apply(P, 2, quantile, probs = c(0.025, 0.5, 0.975))
  q_R <- apply(R, 2, quantile, probs = c(0.025, 0.5, 0.975))

  summary.table <- matrix("", nrow = nteams, ncol = 9)

  colnames(summary.table) <- c(
    "Team", "Exp-P", "SD-P",
    "95%LB-P", "Med-P", "95%UB-P",
    "95%LB-R", "Med-R", "95%UB-R"
  )

  rownames(summary.table) <- 1:nteams

  summary.table[, "Team"]  <- teams
  summary.table[, "Exp-P"] <- round(mean_P, 1)
  summary.table[, "SD-P"]  <- round(sd_P, 1)

  summary.table[, c("95%LB-P","Med-P","95%UB-P")] <-
    t(round(q_P, 1))

  summary.table[, c("95%LB-R","Med-R","95%UB-R")] <-
    t(round(q_R, 1))

  return(summary.table)

}

ranking_probability_table <- function(sim_array){

  # Extract rankings from sim array
  R <- sim_array[, , "R"]
  teams <- dimnames(sim_array)[[2]]
  nteams <- length(teams)

  # Compute frequency tables for each team
  res2 <- apply(R, 2, table)

  # Initialize ranking probability table
  tabrankings <- matrix(0, nteams, nteams)
  rownames(tabrankings) <- teams
  colnames(tabrankings) <- 1:nteams

  # Fill probabilities
  for (i in 1:nteams) {
    idx <- as.numeric(names(res2[[i]]))
    tabrankings[i, idx] <- round(100 * res2[[i]] / B, 1)
  }

  tabrankings
}

plot_heatmap <- function(tabranking){

  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)

  tabranking_df <- tabranking %>%
    as.data.frame() %>%
    rownames_to_column("Team")


  df_long <- tabranking_df %>%
    pivot_longer(
      cols = -Team,
      names_to = "Rank",
      values_to = "Probability"
    ) %>%
    mutate(
      Rank = as.integer(Rank),
      Team = factor(Team, levels = rev(rownames(tabranking)))
    )


  ggplot(df_long, aes(x = Rank, y = Team, fill = Probability)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(Probability >= 1, sprintf("%.1f", Probability), "")),
              size = 2.5) +
    scale_x_continuous(breaks = 1:ncol(tabranking), position = "top") +
    scale_fill_gradient(low = "white", high = "darkred", name = "Probability (%)") +
    labs(
      title = "",
      x = "Final Rank",
      y = "Team"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0),
      panel.grid = element_blank(),
      legend.position = c(0.9, 0.8)
    )

}

plot_points_with_errobars <- function(summary_table, overall_df){
  summary_table = summary_tab

  plot_data <- as.data.frame(summary_table) %>%
    mutate(across(2:last_col(), as.numeric))

  plot_data$Pts <- overall_df$Pts[match(plot_data$Team, overall_df$Squad)]

  ggplot(plot_data, aes(y = reorder(Team, Pts))) +  # Reorder by points ascending
    geom_errorbar(aes(xmin = `95%LB-P`, xmax = `95%UB-P`),
                  width = 0.2, color = "gray60", linewidth = 0.7) +
    geom_point(aes(x = `Med-P`), size = 3, fill = "black", color = "black", shape = 21) +
    geom_point(aes(x = Pts), size = 3, fill = "white", color = "black", shape = 21) +
    geom_segment(aes(y = reorder(Team, Pts), yend = reorder(Team, Pts), x = `Med-P`, xend = Pts),
                 linetype = "dashed", alpha = 0.5, color = "gray40") +
    labs(
      title = "Observed vs Expected Points with 95% Confidence Intervals",
      y = "Team",
      x = "Points",
      caption = "White points: Observed points | Filled points: Expected points"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(hjust = 1),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    ) +
    scale_x_continuous(limits = c(0, NA))
}

gf_df <- function(sim_array, overall_df){

  # 1. Match model team order to league table order
  idx <- match(overall_df$Squad, dimnames(sim_array)[[2]])

  # 2. Extract GF simulations in correct team order
  sim_GF_ordered <- sim_array[, idx, "GF"]

  # 3. Compute quantiles
  GF_quants <- t(apply(sim_GF_ordered, 2, quantile,
                       probs = c(0.025, 0.5, 0.975)))

  # 4. Combine with actual GF
  result_GF <- cbind(
    GF_quants,
    Actual_GF = overall_df$GF
  )

  result_GF
}

ga_df <- function(sim_array, overall_df){

  # 1. Match model team order to league table order
  idx <- match(overall_df$Squad, dimnames(sim_array)[[2]])

  # 2. Extract GF simulations in correct team order
  sim_GA_ordered <- sim_array[, idx, "GA"]

  # 3. Compute quantiles
  GA_quants <- t(apply(sim_GA_ordered, 2, quantile,
                       probs = c(0.025, 0.5, 0.975)))

  # 4. Combine with actual GF
  result_GA <- cbind(
    GA_quants,
    Actual_GA = overall_df$GA
  )

  result_GA
}

plot_gf_ga <- function(sim_array, overall_df) {

  # --- Compute GF and GA tables ---
  gf_df <- as.data.frame(gf_df(sim_array, overall_df))
  ga_df <- as.data.frame(ga_df(sim_array, overall_df))

  # Add team names
  gf_df$Team <- rownames(gf_df)
  ga_df$Team <- rownames(ga_df)

  # --- Convert to long format ---
  gf_long <- gf_df %>%
    mutate(metric = "GF") %>%
    select(Team, lower = `2.5%`, median = `50%`, upper = `97.5%`,
           actual = Actual_GF, metric)

  ga_long <- ga_df %>%
    mutate(metric = "GA") %>%
    select(Team, lower = `2.5%`, median = `50%`, upper = `97.5%`,
           actual = Actual_GA, metric)

  # Combine into long dataset
  plot_df <- bind_rows(gf_long, ga_long)

  # Reverse order of teams (same order as GF table, reversed)
  plot_df$Team <- factor(plot_df$Team, levels = rev(gf_df$Team))

  # --- Plot ---
  p <- ggplot(plot_df, aes(y = Team)) +
    geom_errorbar(aes(xmin = lower, xmax = upper),
                  width = 0.25, color = "gray60") +
    geom_point(aes(x = median),
               size = 3, shape = 21, fill = "black", color = "black") +
    geom_point(aes(x = actual),
               size = 3, shape = 21, fill = "white", color = "black") +
    geom_segment(aes(x = median, xend = actual,
                     y = Team,  yend = Team),
                 linetype = "dashed", color = "gray40") +
    facet_wrap(~ metric, scales = "free_x") +
    labs(
      title = "Observed vs Expected Goals For and Against (95% CI)",
      x = "Goals",
      y = "Team",
      caption = "White: Actual | Black: Expected"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      strip.text = element_text(size = 16, face = "bold"),
      axis.text.y = element_text(size = 12),
      panel.grid.minor = element_blank()
    )

  return(p)
}

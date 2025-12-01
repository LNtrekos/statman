sample_Z_pairs_all_matches <- function(grid_list, Z1_grid, Z2_grid) {

  n <- length(grid_list)
  u <- runif(n)

  idx <- vapply(
    seq_len(n),
    function(j) {
      findInterval(u[j], grid_list[[j]]$cdf, rightmost.closed = TRUE)
    },
    integer(1)
  )

  cbind(Z1_grid[idx], Z2_grid[idx])
}

league_fast <- function(Z, HG1, AG1, HG2, AG2, home_i, away_i, nteams) {
  # Points
  hp <- 3 * (Z > 0) + (Z == 0)
  ap <- 3 * (Z < 0) + (Z == 0)

  points <- rowsum(hp, home_i)[, 1] + rowsum(ap, away_i)[, 1]

  # Goals for / against (combined from both Skellam components)
  GF <- rowsum(HG1 + HG2, home_i)[, 1] +
    rowsum(AG1 + AG2, away_i)[, 1]

  GA <- rowsum(AG1 + AG2, home_i)[, 1] +
    rowsum(HG1 + HG2, away_i)[, 1]

  pos <- nteams + 1 - rank(points, ties.method = "max")

  list(points = points, pos = pos, GF = GF, GA = GA, GD = GF - GA)
}

simulate_seasons_cop_skellam2 <- function(
    B, mus, lambda1_h_hat, lambda1_a_hat,
    lambda2_h_hat, lambda2_a_hat,
    Z1_grid, Z2_grid, grid_list,
    home_i, away_i, teams
) {
  nteams   <- length(teams)
  nmatches <- length(mus)
  metrics  <- c("P", "R", "GF", "GA", "GD")

  sim_array <- array(
    0, dim = c(B, nteams, length(metrics)),
    dimnames = list(NULL, teams, metrics)
  )

  for (b in seq_len(B)) {
    if (b %% 1000L == 0L) {
      cat("Simulation", b, "of", B, "\n")
    }

    # 1) Sample all (Z1,Z2) for all matches in this season
    Z_pairs <- sample_Z_pairs_all_matches(grid_list, Z1_grid, Z2_grid)
    Z_sim   <- Z_pairs[, 1] + Z_pairs[, 2]

    # 2) Sample Poisson goals
    HG1_sim <- rpois(nmatches, lambda1_h_hat)
    AG1_sim <- rpois(nmatches, lambda1_a_hat)
    HG2_sim <- rpois(nmatches, lambda2_h_hat)
    AG2_sim <- rpois(nmatches, lambda2_a_hat)

    # 3) Aggregate into table
    res <- league_fast(
      Z   = Z_sim, HG1 = HG1_sim, AG1 = AG1_sim, HG2 = HG2_sim,
      AG2 = AG2_sim, home_i = home_i, away_i = away_i,
      nteams = nteams
    )

    sim_array[b, , "P"]  <- res$points
    sim_array[b, , "R"]  <- res$pos
    sim_array[b, , "GF"] <- res$GF
    sim_array[b, , "GA"] <- res$GA
    sim_array[b, , "GD"] <- res$GD
  }

  sim_array
}

# Structural scope condition: the perceptual effect of class segregation
# should depend on how much of total wage variation lies between classes
# versus within them. Holding the between-class wage gap fixed at the
# main-text values (mu_L=100, mu_M=200, mu_H=300), this sweeps within-class
# dispersion (sigma) from small (classes are nearly homogeneous, so class
# membership carries almost all the wage information) to large (individual
# variation swamps the class gap), at both theta_H=0 (random mixing) and
# theta_H=5 (strong segregation, as in the main-text class-size sweep).
# The gap between the two curves is the segregation effect at each sigma;
# it should be largest when sigma is small and shrink as sigma grows.

suppressPackageStartupMessages({
  library(statnet); library(network); library(ergm); library(parallel)
  library(ggplot2); library(dplyr)
})

set.seed(20250906)
n_cores <- max(1, detectCores() - 1)
base_wage <- c(low = 100, mid = 200, high = 300)

ergm_control <- control.simulate.formula(MCMC.burnin = 200000, MCMC.interval = 1)

simulate_one_sigma <- function(N, sigma_log, homophily, density = 15 / (N - 1)) {
  low_n <- round(N / 3); mid_n <- round(N / 3); high_n <- N - low_n - mid_n
  class <- c(rep("low", low_n), rep("mid", mid_n), rep("high", high_n))

  logwage <- rnorm(N, mean = log(base_wage[class]) - sigma_log^2 / 2, sd = sigma_log)
  wage <- exp(logwage)

  e <- round((N * (N - 1) / 2) * density)
  net <- as.network(sna::rgnm(1, N, e, mode = "graph"), directed = FALSE)
  net %v% "class" <- class

  sim_net <- simulate(net ~ edges + nodematch("class", diff = FALSE),
                       constraints = ~degreedist, coef = c(0, homophily), nsim = 1,
                       control = ergm_control)
  adj <- as.matrix.network.adjacency(sim_net)

  perceived_ratio <- sapply(seq_len(N), function(i) {
    nbrs <- which(adj[i, ] == 1)
    k <- length(nbrs)
    if (k < 3) return(1)
    nbr_wage <- sort(wage[nbrs]); third <- floor(k / 3)
    mean(nbr_wage[(k - third + 1):k]) / mean(nbr_wage[1:third])
  })

  data.frame(sigma_log = sigma_log, homophily = homophily, perceived_ratio = perceived_ratio)
}

run_reps_parallel <- function(reps, ...) {
  out <- mclapply(seq_len(reps), function(r) simulate_one_sigma(...), mc.cores = n_cores)
  do.call(rbind, out)
}

se <- function(x) sd(x) / sqrt(length(x))

sigma_seq <- c(0.05, 0.10, 0.15, 0.25, 0.40, 0.60)
homophily_seq <- c(0, 5)
N1 <- 300; reps <- 50

cat("=== Structural scope condition: sigma sweep at theta_H in {0, 5} ===\n")
t0 <- Sys.time()
grid <- expand.grid(sigma_log = sigma_seq, homophily = homophily_seq)
sweep <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  run_reps_parallel(reps = reps, N = N1, sigma_log = grid$sigma_log[i], homophily = grid$homophily[i])
}))
cat("Runtime:", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "sec (",
    nrow(grid) * reps, "simulations )\n\n")

agg <- sweep %>%
  group_by(sigma_log, homophily) %>%
  summarise(ratio_mean = mean(perceived_ratio), ratio_se = se(perceived_ratio), .groups = "drop")

wide <- agg %>%
  tidyr::pivot_wider(id_cols = sigma_log, names_from = homophily, values_from = ratio_mean, names_prefix = "theta_") %>%
  mutate(segregation_gap = theta_5 - theta_0)

cat("Perceived ratio by sigma, at random mixing (theta_H=0) vs. strong segregation (theta_H=5):\n")
print(wide, n = 20)

write.csv(agg, "data/scope_condition_sigma_sweep.csv", row.names = FALSE)
cat("\nSaved data/scope_condition_sigma_sweep.csv\n")

p <- ggplot(agg, aes(x = sigma_log, y = ratio_mean, color = factor(homophily))) +
  geom_ribbon(aes(ymin = ratio_mean - 1.96 * ratio_se, ymax = ratio_mean + 1.96 * ratio_se, fill = factor(homophily)),
              alpha = 0.15, linewidth = 0) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
  scale_color_manual(values = c("0" = "grey50", "5" = "#d73027"),
                      labels = c(expression(paste("Random mixing (", theta[H], "=0)")),
                                 expression(paste("Strong segregation (", theta[H], "=5)")))) +
  scale_fill_manual(values = c("0" = "grey50", "5" = "#d73027"), guide = "none") +
  labs(title = "Structural scope condition: within- vs. between-class dispersion",
       subtitle = expression(paste("Between-class gap fixed (", mu[L], "=100, ", mu[M], "=200, ", mu[H],
                                    "=300); within-class dispersion (", sigma, ") varies")),
       x = expression(paste("Within-class dispersion (", sigma, ", log scale)")),
       y = "Perceived pay ratio", color = NULL) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10.5))

ggsave("results/scope_condition_sigma_sweep.png", p, width = 8.5, height = 5.5, dpi = 300, bg = "white")
cat("Saved results/scope_condition_sigma_sweep.png\n")

# Robustness check: does the compression pattern hold across a plausible
# range of personal-network sizes? Rather than an arbitrary bracket, this
# uses three of Dunbar's own named layers (Dunbar 1992, 1993): the "support
# clique" (~5, closest confidants), the "sympathy group" (~15, the main-text
# anchor), and the "band" (~50, a wider circle of regular contacts).
# N is held fixed at 300; density is varied so that expected degree hits
# each target (5, 15, 50). Ties depend on class alone (no continuous
# income-distance term), matching the main-text specification.

suppressPackageStartupMessages({
  library(statnet); library(network); library(ergm); library(parallel)
  library(ggplot2); library(dplyr); library(igraph)
})
set.seed(20250906)
n_cores <- max(1, detectCores() - 1)
base_wage <- c(low = 100, mid = 200, high = 300)
sigma_log <- 0.15

# Verified to reach target mixing patterns reliably at N=300 across this
# degree range (up to ~7,500 edges at degree=50); do not reduce without
# re-checking convergence (see 01_simulation_main.R header for diagnostics).
ergm_control <- control.simulate.formula(MCMC.burnin = 400000, MCMC.interval = 1)

simulate_one <- function(N, homophily, density) {
  low_n <- round(N / 3); mid_n <- round(N / 3); high_n <- N - low_n - mid_n
  class <- c(rep("low", low_n), rep("mid", mid_n), rep("high", high_n))
  logwage <- rnorm(N, mean = log(base_wage[class]) - sigma_log^2 / 2, sd = sigma_log)
  wage <- exp(logwage)
  e <- round((N * (N - 1) / 2) * density)
  net <- as.network(sna::rgnm(1, N, e, mode = "graph"), directed = FALSE)
  net %v% "class" <- class; net %v% "logwage" <- logwage
  sim_net <- simulate(net ~ edges + nodematch("class", diff = FALSE),
                       constraints = ~degreedist, coef = c(0, homophily), nsim = 1,
                       control = ergm_control)
  adj <- as.matrix.network.adjacency(sim_net)
  g <- graph_from_adjacency_matrix(adj, mode = "undirected")
  realized_seg <- tryCatch(assortativity_nominal(g, as.integer(factor(class)), directed = FALSE), error = function(e) NA_real_)

  ratio <- sapply(seq_len(N), function(i) {
    nbrs <- which(adj[i, ] == 1)
    nbr_wage <- wage[nbrs]; k <- length(nbr_wage)
    if (k < 3) return(1)
    s <- sort(nbr_wage); third <- floor(k / 3)
    mean(s[(k - third + 1):k]) / mean(s[1:third])
  })
  data.frame(ratio = ratio, homophily = homophily, degree_target = round(density * (N - 1)), realized_seg = realized_seg)
}

run_reps <- function(reps, ...) {
  out <- mclapply(seq_len(reps), function(r) simulate_one(...), mc.cores = n_cores)
  bad <- sapply(out, function(x) !is.data.frame(x))
  if (any(bad)) out[bad] <- lapply(seq_len(sum(bad)), function(i) simulate_one(...))
  do.call(rbind, out)
}

N <- 300
degree_targets <- c(5, 15, 50)   # Dunbar's support clique, sympathy group (main text), band
homophily_grid <- c(0, 1, 2, 3, 4, 5)
reps <- 30

t0 <- Sys.time()
results <- do.call(rbind, lapply(degree_targets, function(dg) {
  dens <- dg / (N - 1)
  do.call(rbind, lapply(homophily_grid, function(h) run_reps(reps, N = N, homophily = h, density = dens)))
}))
cat("Runtime:", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "sec (",
    length(degree_targets) * length(homophily_grid) * reps, "simulations )\n\n")

agg <- results %>%
  group_by(degree_target, homophily) %>%
  summarise(ratio_mean = mean(ratio), ratio_se = sd(ratio) / sqrt(n()),
            realized_seg_mean = mean(realized_seg, na.rm = TRUE), .groups = "drop")

write.csv(agg, "data/degree_sensitivity.csv", row.names = FALSE)
print(agg, n = 30)

p <- ggplot(agg, aes(x = homophily, y = ratio_mean, color = factor(degree_target))) +
  geom_ribbon(aes(ymin = ratio_mean - 1.96*ratio_se, ymax = ratio_mean + 1.96*ratio_se, fill = factor(degree_target)), alpha = 0.12, linewidth = 0) +
  geom_hline(yintercept = base_wage["high"]/base_wage["low"], linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
  scale_color_manual(values = c("5" = "#66c2a5", "15" = "#3288bd", "50" = "#5e4fa2"),
                      labels = c("5" = "Degree ~5 (support clique)", "15" = "Degree ~15 (sympathy group, main text)", "50" = "Degree ~50 (band)"),
                      name = NULL) +
  scale_fill_manual(values = c("5" = "#66c2a5", "15" = "#3288bd", "50" = "#5e4fa2"), guide = "none") +
  labs(title = "Robustness to network size (expected degree)",
       subtitle = "Bracketing the main-text estimate with two of Dunbar's (1992, 1993) other named layers",
       x = "Categorical class homophily (theta[H])", y = "Perceived pay ratio (all classes pooled)") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 13))

ggsave("results/sim_degree_sensitivity.png",
       p, width = 8, height = 5.5, dpi = 300, bg = "white")
cat("\nSaved sim_degree_sensitivity.png\n")

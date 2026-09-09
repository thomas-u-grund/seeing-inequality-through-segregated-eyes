# Extension: what happens if continuous income-distance homophily is layered
# on top of the main-text class-only model? theta_D=0 reproduces the
# main-text specification (ties depend on class alone); theta_D=-1,-2,-4
# progressively add a continuous income-proximity tie-formation preference
# on top of it, to check whether this additional, more demanding mechanism
# changes the qualitative compression pattern or merely deepens it.

suppressPackageStartupMessages({
  library(statnet); library(network); library(ergm); library(parallel)
  library(ggplot2); library(dplyr); library(igraph)
})
set.seed(20250906)
n_cores <- max(1, detectCores() - 1)
base_wage <- c(low = 100, mid = 200, high = 300)
sigma_log <- 0.15

# Verified to reach target mixing patterns reliably at N=300/~2,250 edges;
# do not reduce without re-checking convergence (see 01_simulation_main.R).
ergm_control <- control.simulate.formula(MCMC.burnin = 200000, MCMC.interval = 1)

simulate_one <- function(N, homophily, lambda_dist, density = 15 / (N - 1)) {
  low_n <- round(N / 3); mid_n <- round(N / 3); high_n <- N - low_n - mid_n
  class <- c(rep("low", low_n), rep("mid", mid_n), rep("high", high_n))
  logwage <- rnorm(N, mean = log(base_wage[class]) - sigma_log^2 / 2, sd = sigma_log)
  wage <- exp(logwage)
  e <- round((N * (N - 1) / 2) * density)
  net <- as.network(sna::rgnm(1, N, e, mode = "graph"), directed = FALSE)
  net %v% "class" <- class; net %v% "logwage" <- logwage
  sim_net <- simulate(net ~ edges + nodematch("class", diff = FALSE) + absdiff("logwage"),
                       constraints = ~degreedist, coef = c(0, homophily, lambda_dist), nsim = 1,
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
  data.frame(ratio = ratio, homophily = homophily, lambda_dist = lambda_dist, realized_seg = realized_seg)
}

run_reps <- function(reps, ...) {
  out <- mclapply(seq_len(reps), function(r) simulate_one(...), mc.cores = n_cores)
  bad <- sapply(out, function(x) inherits(x, "try-error") || is.character(x) || !is.data.frame(x))
  if (any(bad)) {
    cat("WARNING:", sum(bad), "of", reps, "worker(s) failed; retrying those serially.\n")
    out[bad] <- lapply(seq_len(sum(bad)), function(i) simulate_one(...))
  }
  do.call(rbind, out)
}

homophily_grid <- c(0, 1, 2, 3, 4, 5)
lambda_grid <- c(0, -1, -2, -4)
reps <- 30; N <- 300  # matches main simulation's N/density (degree ~15, Dunbar's sympathy group)

t0 <- Sys.time()
results <- do.call(rbind, lapply(lambda_grid, function(ld) {
  do.call(rbind, lapply(homophily_grid, function(h) run_reps(reps, N = N, homophily = h, lambda_dist = ld)))
}))
cat("Runtime:", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "sec (",
    length(homophily_grid) * length(lambda_grid) * reps, "simulations )\n\n")

agg <- results %>%
  group_by(lambda_dist, homophily) %>%
  summarise(ratio_mean = mean(ratio), ratio_se = sd(ratio) / sqrt(n()),
            realized_seg_mean = mean(realized_seg, na.rm = TRUE), .groups = "drop") %>%
  mutate(lambda_label = factor(paste0("theta[D]==", lambda_dist),
                                levels = paste0("theta[D]==", lambda_grid)))

write.csv(agg, "data/theta_d_sensitivity.csv", row.names = FALSE)
print(agg, n = 30)

p <- ggplot(agg, aes(x = homophily, y = ratio_mean, color = factor(lambda_dist))) +
  geom_ribbon(aes(ymin = ratio_mean - 1.96*ratio_se, ymax = ratio_mean + 1.96*ratio_se, fill = factor(lambda_dist)), alpha = 0.12, linewidth = 0) +
  geom_hline(yintercept = base_wage["high"]/base_wage["low"], linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
  scale_color_manual(values = c("0" = "#999999", "-1" = "#66c2a5", "-2" = "#3288bd", "-4" = "#5e4fa2"),
                      labels = c("0" = "theta[D] = 0 (main text: class only)", "-1" = "theta[D] = -1 (extension)", "-2" = "theta[D] = -2 (extension)", "-4" = "theta[D] = -4 (extension)"),
                      name = NULL) +
  scale_fill_manual(values = c("0" = "#999999", "-1" = "#66c2a5", "-2" = "#3288bd", "-4" = "#5e4fa2"), guide = "none") +
  labs(title = "Extension: adding continuous income-distance homophily (theta[D])",
       subtitle = "theta[D]=0 is the main-text class-only model; theta[D]<0 layers income proximity on top",
       x = "Categorical class homophily (theta[H])", y = "Perceived pay ratio (all classes pooled)") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 13))

ggsave("results/sim_thetaD_sensitivity.png",
       p, width = 8, height = 5.5, dpi = 300, bg = "white")
cat("\nSaved sim_thetaD_sensitivity.png\n")

# Final simulation for "Seeing Inequality Through Segregated Eyes"
#
# v6: perception is no longer class-conditional at all, following the
# original (2014) model this project is built on. Earlier versions asked
# each agent to estimate a specific class's mean from same-class contacts,
# which required a fallback rule for agents with zero contacts in that
# class -- an awkward, behaviorally strong assumption (why would someone
# with no personal contact in a class have literally no outside
# information about it?). The original model never had this problem
# because it never asked agents to estimate individual classes: each
# agent's perception is simply a statistic of the wage composition of
# their OWN local network, with no class filtering at all. We recover
# this here, generalized to handle within-class wage heterogeneity
# (necessary so that a single contact does not reveal a class's wage with
# certainty, see Section 3.1):
#   - Perceived pay ratio (H1, H3): the ratio of the mean wage among ego's
#     highest-earning third of contacts to their lowest-earning third,
#     ranked purely by income -- no class labels involved in the split.
#   - Perceived population mean (H2): the raw, unfiltered mean wage among
#     ego's contacts, compared against the true population mean.
# Both quantities are always well-defined for any agent with at least a
# handful of contacts (never zero at the degree levels used below), so
# there is no missing-information edge case to assume away. Ties still
# depend on class homophily (edges + nodematch("class")) as in v5; an
# extension layering continuous income-distance homophily on top is
# reported separately (theta_d_sensitivity.R).
#
# MEAN-PRESERVING lognormal wages (retained from v4). log(y) ~ N(log(mu_c), sigma^2) makes
#    mu_c the *median*, not the mean (E[y] = mu_c * exp(sigma^2/2)). Fixed:
#    log(y) ~ N(log(mu_c) - sigma^2/2, sigma^2), so E[y_i] = mu_c exactly.
#
# 3. CHANCE-CORRECTED realized segregation. Raw same-class tie share has a
#    baseline that mechanically shifts with class-size composition (e.g.
#    .333 under equal thirds vs .375 under a 50/25/25 split at pure random
#    mixing) -- not appropriate for the class-size sweep, where the whole
#    point is to hold "segregation" fixed while composition varies. Realized
#    segregation is now Newman's nominal assortativity coefficient (via
#    igraph), which is 0 at random mixing regardless of class-size
#    composition and 1 at complete segregation, by construction. Raw
#    same-class share is still reported alongside it for intuition.
#
# 4. N=300 / degree=15 (not N=3000 / degree=150). ergm's degree-preserving
#    (CondDegreeDist) MCMC needs a burnin scaling with the number of tied
#    pairs to actually reach the target mixing pattern; at N=3000/deg=150
#    (~225,000 edges) the default burnin (and even a 5e6-step explicit one,
#    ~3 min/simulation) is impractical to run at the replicate counts used
#    below. N=300/degree=15 (~2,250 edges) converges reliably in ~0.1-0.2s
#    with an explicit, verified burnin, at the SAME density (~0.05) used
#    throughout. Degree=15 is anchored to Dunbar's (1993) "sympathy group"
#    layer (his own estimate of a person's closest ~12-15 relationships)
#    rather than the broader ~150-tie "stable relationships" layer used in
#    an earlier version of this script; for the purpose of this model --
#    inferring others' incomes with some accuracy -- the smaller, closer-tie
#    layer is arguably the more defensible anchor in any case. A bracketing
#    robustness check across three of Dunbar's named layers (support clique
#    ~5, sympathy group ~15, band ~50) is reported in the Appendix.

suppressPackageStartupMessages({
  library(statnet); library(network); library(ergm); library(ineq)
  library(parallel); library(ggplot2); library(dplyr); library(igraph)
})

set.seed(20250906)
n_cores <- max(1, detectCores() - 1)

base_wage <- c(low = 100, mid = 200, high = 300)
sigma_log <- 0.15

# Verified (see script header) to reach target mixing patterns reliably at
# N=300/~2,250 edges; do not reduce without re-checking convergence.
ergm_control <- control.simulate.formula(MCMC.burnin = 200000, MCMC.interval = 1)

simulate_one <- function(N, low_share, mid_share, high_share, homophily, density = 15 / (N - 1)) {
  low_n  <- round(N * low_share)
  mid_n  <- round(N * mid_share)
  high_n <- N - low_n - mid_n
  class  <- c(rep("low", low_n), rep("mid", mid_n), rep("high", high_n))

  # mean-preserving lognormal: E[y_i] = base_wage[class] exactly
  logwage <- rnorm(N, mean = log(base_wage[class]) - sigma_log^2 / 2, sd = sigma_log)
  wage <- exp(logwage)

  # population-level analogue of the local ratio below (same top/bottom-third
  # construction, applied to the whole simulated population rather than to
  # ego's contacts) -- the appropriate "accurate general knowledge" benchmark
  # for the omega extension, since it is the same estimand as R_i^L rather
  # than the class-mean ratio mu_H/mu_L
  s_pop <- sort(wage); third_pop <- floor(N / 3)
  true_ratio_pop <- mean(s_pop[(N - third_pop + 1):N]) / mean(s_pop[1:third_pop])

  e <- round((N * (N - 1) / 2) * density)
  net <- as.network(sna::rgnm(1, N, e, mode = "graph"), directed = FALSE)
  net %v% "class" <- class
  net %v% "logwage" <- logwage

  sim_net <- simulate(net ~ edges + nodematch("class", diff = FALSE),
                       constraints = ~degreedist, coef = c(0, homophily), nsim = 1,
                       control = ergm_control)

  adj <- as.matrix.network.adjacency(sim_net)
  total_edges <- sum(adj[upper.tri(adj)])
  same_class_mat <- outer(class, class, "==")
  same_class_edges <- sum(adj[upper.tri(adj)] & same_class_mat[upper.tri(same_class_mat)])
  raw_same_class_share <- if (total_edges > 0) same_class_edges / total_edges else NA

  g <- graph_from_adjacency_matrix(adj, mode = "undirected")
  realized_segregation <- tryCatch(
    assortativity_nominal(g, as.integer(factor(class)), directed = FALSE),
    error = function(e) NA_real_
  )

  perceived <- lapply(seq_len(N), function(i) {
    nbrs <- which(adj[i, ] == 1)
    k <- length(nbrs)
    if (k == 0) {
      list(mean_wage = wage[i], iov = 0, ratio = 1)
    } else {
      nbr_wage <- wage[nbrs]
      own_nbr_mean <- mean(nbr_wage)
      # rank-based split of ego's own contacts into top/bottom thirds by
      # income -- no class labels used; undefined (neutral, ratio = 1) for
      # the handful of agents with fewer than 3 contacts
      ratio <- if (k >= 3) {
        s <- sort(nbr_wage); third <- floor(k / 3)
        mean(s[(k - third + 1):k]) / mean(s[1:third])
      } else 1
      list(mean_wage = own_nbr_mean,
           iov = if (k > 1) (sd(nbr_wage) / mean(nbr_wage))^2 else 0,
           ratio = ratio)
    }
  })

  perceived_mean  <- sapply(perceived, `[[`, "mean_wage")
  perceived_iov   <- sapply(perceived, `[[`, "iov")
  perceived_ratio <- sapply(perceived, `[[`, "ratio")

  k_perception <- 1.0
  k_selfinterest <- 0.9
  own_wage_z <- (wage - mean(wage)) / sd(wage)
  attitude_latent <- k_perception * log(perceived_ratio) - k_selfinterest * own_wage_z

  data.frame(
    class = class, wage = wage,
    perceived_mean = perceived_mean, true_mean = mean(wage), bias_mean = perceived_mean - mean(wage),
    perceived_iov = perceived_iov, true_iov = (sd(wage) / mean(wage))^2,
    bias_iov = perceived_iov - (sd(wage) / mean(wage))^2,
    perceived_ratio = perceived_ratio, attitude = attitude_latent,
    homophily = homophily, realized_segregation = realized_segregation,
    raw_same_class_share = raw_same_class_share,
    low_share = low_share, high_share = high_share,
    true_ratio_pop = true_ratio_pop
  )
}

run_reps_parallel <- function(reps, ...) {
  out <- mclapply(seq_len(reps), function(r) simulate_one(...), mc.cores = n_cores)
  do.call(rbind, out)
}

se <- function(x) sd(x) / sqrt(length(x))

cat("Using", n_cores, "cores.\n\n")

# ============================================================
# Sweep 1 (MAIN): segregation, random mixing to strong segregation
# ============================================================
cat("=== Sweep 1 (main): beta from 0 (random) to 5 (strong segregation) ===\n")
homophily_seq_main <- seq(0, 5, by = 0.25)
# N and density chosen so that expected degree ~15, matching Dunbar's (1993)
# "sympathy group" estimate of a person's closest circle of relationships.
# At N=300, a degree of 15 is 5 percent of the simulated population -- the
# same "small window on a large society" property as a larger N with a
# proportionally larger degree target, but at a small enough scale (~2,250
# realized edges) that the degree-preserving MCMC below converges reliably
# within a verified, explicit burnin (see script header).
N1 <- 300; reps1 <- 50

t0 <- Sys.time()
sweep1 <- do.call(rbind, lapply(homophily_seq_main, function(h) {
  run_reps_parallel(reps = reps1, N = N1, low_share = 1/3, mid_share = 1/3, high_share = 1/3, homophily = h)
}))
cat("Sweep 1 runtime:", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "sec (",
    length(homophily_seq_main) * reps1, "simulations )\n\n")

agg1 <- sweep1 %>%
  group_by(class, homophily) %>%
  summarise(ratio_mean = mean(perceived_ratio), ratio_se = se(perceived_ratio),
            attitude_mean = mean(attitude), attitude_se = se(attitude),
            bias_mean_m = mean(bias_mean), bias_iov_m = mean(bias_iov),
            realized_seg_mean = mean(realized_segregation, na.rm = TRUE),
            raw_share_mean = mean(raw_same_class_share, na.rm = TRUE), .groups = "drop")

write.csv(agg1, "data/simulation_segregation_sweep.csv", row.names = FALSE)

# ============================================================
# Sweep 1b (SUPPLEMENTARY): heterophily boundary check, beta -5 to 0
# ============================================================
cat("=== Sweep 1b (supplementary): beta from -5 (heterophily) to 0 ===\n")
homophily_seq_supp <- seq(-5, 0, by = 1)
reps1b <- 25
t0 <- Sys.time()
sweep1b <- do.call(rbind, lapply(homophily_seq_supp, function(h) {
  run_reps_parallel(reps = reps1b, N = N1, low_share = 1/3, mid_share = 1/3, high_share = 1/3, homophily = h)
}))
cat("Sweep 1b runtime:", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "sec\n\n")

agg1b <- sweep1b %>%
  group_by(class, homophily) %>%
  summarise(ratio_mean = mean(perceived_ratio), ratio_se = se(perceived_ratio),
            attitude_mean = mean(attitude), attitude_se = se(attitude),
            realized_seg_mean = mean(realized_segregation, na.rm = TRUE), .groups = "drop")
write.csv(agg1b, "data/simulation_heterophily_supplement.csv", row.names = FALSE)

# ============================================================
# Sweep 2: class size, high segregation fixed
# ============================================================
cat("=== Sweep 2: class size, high segregation fixed (beta=5) ===\n")
size_configs <- list(c(1/3,1/3,1/3), c(0.5,0.25,0.25), c(0.25,0.5,0.25), c(0.25,0.25,0.5))
reps2 <- 60

t0 <- Sys.time()
sweep2 <- do.call(rbind, lapply(size_configs, function(cfg) {
  run_reps_parallel(reps = reps2, N = N1, low_share = cfg[1], mid_share = cfg[2], high_share = cfg[3], homophily = 5)
}))
cat("Sweep 2 runtime:", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "sec (",
    length(size_configs) * reps2, "simulations )\n\n")

agg2 <- sweep2 %>%
  mutate(config = paste0("L", round(low_share*100), "/H", round(high_share*100))) %>%
  group_by(class, config, low_share, high_share) %>%
  summarise(bias_mean_m = mean(bias_mean), bias_mean_se = se(bias_mean),
            attitude_mean = mean(attitude), attitude_se = se(attitude),
            realized_seg_mean = mean(realized_segregation, na.rm = TRUE),
            raw_share_mean = mean(raw_same_class_share, na.rm = TRUE), .groups = "drop")

configs_order <- c("L25/H25", "L25/H50", "L33/H33", "L50/H25")
agg2$config <- factor(agg2$config, levels = configs_order)
write.csv(agg2, "data/simulation_classsize_sweep.csv", row.names = FALSE)

# ============================================================
# Figures (ggplot2, publication style) -- x-axis is chance-corrected
# realized segregation (assortativity)
# ============================================================
cols <- c(low = "#1b7837", mid = "#7570b3", high = "#d73027")
class_labels <- c(low = "Low class", mid = "Middle class", high = "High class")

# population-level benchmark: same top-third/bottom-third statistic as the
# individual-level perceived ratio, applied to the whole simulated
# population (matches R^G in the omega extension) -- not the class-mean
# ratio mu_H/mu_L, which is a different (if numerically close) quantity
true_ratio_pop_mean <- mean(sweep1$true_ratio_pop)

theme_pub <- theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, color = "grey88"),
        legend.position = "bottom", legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        strip.background = element_blank())

p1a <- ggplot(agg1, aes(x = realized_seg_mean, y = ratio_mean, color = class, fill = class)) +
  geom_ribbon(aes(ymin = ratio_mean - 1.96*ratio_se, ymax = ratio_mean + 1.96*ratio_se), alpha = 0.15, linewidth = 0) +
  geom_hline(yintercept = true_ratio_pop_mean, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = cols, labels = class_labels) +
  scale_fill_manual(values = cols, labels = class_labels) +
  labs(title = "Perceived inequality vs. realized segregation", x = "Realized segregation (nominal assortativity)", y = "Perceived pay ratio") +
  annotate("text", x = 0.02, y = true_ratio_pop_mean + 0.08, label = "population benchmark", color = "grey40", size = 3.2, hjust = 0) +
  theme_pub

p1b <- ggplot(agg1, aes(x = realized_seg_mean, y = attitude_mean, color = class, fill = class)) +
  geom_ribbon(aes(ymin = attitude_mean - 1.96*attitude_se, ymax = attitude_mean + 1.96*attitude_se), alpha = 0.15, linewidth = 0) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = cols, labels = class_labels) +
  scale_fill_manual(values = cols, labels = class_labels) +
  labs(title = "Attitude vs. realized segregation", x = "Realized segregation (nominal assortativity)", y = "Attitude (support for redistribution)") +
  theme_pub

combined1 <- gridExtra::grid.arrange(p1a, p1b, ncol = 2)
ggsave("results/sim_segregation_ratio_attitude.png",
       combined1, width = 11, height = 5, dpi = 300, bg = "white")

p1c <- ggplot(agg1b, aes(x = homophily, y = ratio_mean, color = class)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
  geom_hline(yintercept = true_ratio_pop_mean, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = cols, labels = class_labels) +
  labs(title = "Supplementary: heterophily boundary (beta < 0)", x = "Homophily coefficient (beta)", y = "Perceived pay ratio") +
  theme_pub
ggsave("results/sim_heterophily_supplement.png",
       p1c, width = 6, height = 4.5, dpi = 300, bg = "white")

p2a <- ggplot(agg2, aes(x = config, y = bias_mean_m, color = class, group = class)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
  geom_errorbar(aes(ymin = bias_mean_m - 1.96*bias_mean_se, ymax = bias_mean_m + 1.96*bias_mean_se), width = 0.08, linewidth = 0.6) +
  geom_line(aes(group = class), linewidth = 0.9) + geom_point(size = 2.2) +
  scale_color_manual(values = cols, labels = class_labels) +
  labs(title = "Minority amplification: perceptual bias", x = "Class-size configuration", y = "Perceived - true wage") +
  theme_pub

p2b <- ggplot(agg2, aes(x = config, y = attitude_mean, color = class, group = class)) +
  geom_errorbar(aes(ymin = attitude_mean - 1.96*attitude_se, ymax = attitude_mean + 1.96*attitude_se), width = 0.08, linewidth = 0.6) +
  geom_line(aes(group = class), linewidth = 0.9) + geom_point(size = 2.2) +
  scale_color_manual(values = cols, labels = class_labels) +
  labs(title = "Minority amplification: attitude", x = "Class-size configuration", y = "Attitude (support for redistribution)") +
  theme_pub

combined2 <- gridExtra::grid.arrange(p2a, p2b, ncol = 2)
ggsave("results/sim_classsize_bias_attitude.png",
       combined2, width = 11, height = 5, dpi = 300, bg = "white")

cat("\n=== Sweep 1 (main) summary, selected points ===\n")
print(agg1 %>% filter(homophily %in% c(0, 1, 2, 3, 4, 5)) %>% arrange(class, homophily) %>%
        select(class, homophily, realized_seg_mean, raw_share_mean, ratio_mean, attitude_mean), n = 30)
cat("\n=== Sweep 2 summary ===\n")
print(agg2 %>% arrange(class, config) %>% select(class, config, realized_seg_mean, raw_share_mean, bias_mean_m, attitude_mean), n = 30)

saveRDS(list(sweep1 = sweep1, agg1 = agg1, sweep1b = sweep1b, agg1b = agg1b, sweep2 = sweep2, agg2 = agg2),
        "data/simulation_results.rds")
cat("\nDone.\n")

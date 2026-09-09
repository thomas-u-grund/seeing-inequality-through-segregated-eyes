# Scope-condition extension (Section 4.3): what happens to the compression
# pattern when perceived inequality blends ego's local network experience
# with accurate, network-independent general knowledge of the true
# population ratio?
#
#   log(R_hat_i) = omega * log(R_i^L) + (1 - omega) * log(R^G)
#
# R_i^L (ego's local, network-based ratio) and R^G (the population-level
# top-third/bottom-third ratio) are both already computed and returned per
# agent by simulate_one() in 01_simulation_main.R, as "perceived_ratio" and
# "true_ratio_pop" respectively -- so this extension is a closed-form
# post-processing step on the raw Sweep 1 output saved in
# data/simulation_results.rds, not a new simulation. This reuses the exact
# same simulated societies (seed 20250906) as the main-text Figure 2, so
# the omega=1 column reproduces that figure's numbers exactly.

suppressPackageStartupMessages({ library(dplyr) })

res <- readRDS("data/simulation_results.rds")
sweep1 <- res$sweep1

se <- function(x) sd(x) / sqrt(length(x))

omega_seq <- c(0, 0.25, 0.5, 0.75, 1)

out <- do.call(rbind, lapply(omega_seq, function(w) {
  sweep1 %>%
    mutate(r_hat = exp(w * log(perceived_ratio) + (1 - w) * log(true_ratio_pop))) %>%
    group_by(homophily) %>%
    summarise(ratio_mean = mean(r_hat), ratio_se = se(r_hat),
              realized_seg_mean = mean(realized_segregation, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(alpha = w)
}))

out <- out %>% select(homophily, ratio_mean, ratio_se, realized_seg_mean, alpha)

cat("=== Omega extension: perceived ratio at homophily = 5, by omega ===\n")
print(out %>% filter(homophily == 5))

write.csv(out, "data/alpha_extension.csv", row.names = FALSE)
cat("\nSaved data/alpha_extension.csv (column 'alpha' = paper's omega).\n")

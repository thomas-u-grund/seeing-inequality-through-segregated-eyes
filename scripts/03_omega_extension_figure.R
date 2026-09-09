# Figure: perceived pay ratio against realized segregation at five weights
# (omega = 0, 0.25, 0.5, 0.75, 1) on local network experience relative to
# accurate general knowledge of the true population ratio (main text,
# Figure 2 / fig:omegaext). Reads data/alpha_extension.csv produced by
# 02_omega_extension.R (column "alpha" = paper's omega).

suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

d <- read.csv("data/alpha_extension.csv")
d$omega <- factor(d$alpha, levels = c(0, 0.25, 0.5, 0.75, 1))

p <- ggplot(d, aes(x = realized_seg_mean, y = ratio_mean, color = omega)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.4) +
  scale_color_brewer(palette = "Blues", direction = 1,
                      labels = c("0", "0.25", "0.5", "0.75", "1 (baseline)")) +
  labs(title = "Local vs. general-knowledge weighting",
       subtitle = expression(paste("Perceived ratio as ", omega, " (weight on local network) varies")),
       x = "Realized segregation (nominal assortativity)",
       y = "Perceived pay ratio", color = expression(omega)) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13))

ggsave("results/sim_alpha_extension.png", p, width = 7, height = 5, dpi = 300, bg = "white")
cat("Saved results/sim_alpha_extension.png\n")

suppressPackageStartupMessages({library(ggplot2); library(dplyr)})

d <- read.csv("data/country_slopes.csv")
d <- d %>% mutate(lo = b - 1.96*se, hi = b + 1.96*se, c_alphan = reorder(c_alphan, b))

p <- ggplot(d, aes(x = b, y = c_alphan)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, color = "#3288bd", linewidth = 0.6) +
  geom_point(color = "#3288bd", size = 2) +
  labs(x = "Country-specific exposure coefficient (H1)", y = NULL,
       title = "Country-specific exposure-perception slopes",
       subtitle = "Within-country OLS, no country fixed effects, own income and demographic controls included") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9.5, color = "grey30"))

ggsave("results/country_slopes.png",
       p, width = 7, height = 7, dpi = 300, bg = "white")
cat("Saved country_slopes.png\n")

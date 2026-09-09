# Two additional figures for the ESR manuscript:
# 1. Example networks at low/medium/high realized segregation, illustrating
#    what the homophily parameter actually does structurally.
# 2. A coefficient plot summarizing the empirical results (Tables 1-2)
#    visually, across all three occupational anchors.

suppressPackageStartupMessages({
  library(statnet); library(network); library(ergm); library(igraph)
  library(ggplot2); library(dplyr); library(tibble)
})
set.seed(42)

# ---------------------------------------------------------------
# Figure: example networks at three homophily levels
# ---------------------------------------------------------------
base_wage <- c(low = 100, mid = 200, high = 300)
sigma_log <- 0.15

draw_example_net <- function(homophily, N = 60, density = 0.06) {
  class <- c(rep("low", N/3), rep("mid", N/3), rep("high", N/3))
  logwage <- rnorm(N, mean = log(base_wage[class]) - sigma_log^2 / 2, sd = sigma_log)
  e <- round((N * (N - 1) / 2) * density)
  net <- as.network(sna::rgnm(1, N, e, mode = "graph"), directed = FALSE)
  net %v% "class" <- class
  net %v% "logwage" <- logwage
  sim_net <- simulate(net ~ edges + nodematch("class", diff = FALSE),
                       constraints = ~degreedist, coef = c(0, homophily), nsim = 1,
                       control = control.simulate.formula(MCMC.burnin = 200000, MCMC.interval = 1))
  adj <- as.matrix.network.adjacency(sim_net)
  g <- graph_from_adjacency_matrix(adj, mode = "undirected")
  V(g)$class <- class
  realized_seg <- tryCatch(assortativity_nominal(g, as.integer(factor(class)), directed = FALSE),
                            error = function(e) NA_real_)
  list(g = g, realized_seg = realized_seg)
}

png("results/sim_example_networks.png",
    width = 3300, height = 1150, res = 300)
par(mfrow = c(1, 3), mar = c(3, 1, 4, 1))
cols <- c(low = "#1b7837", mid = "#7570b3", high = "#d73027")

set.seed(1)
for (h in c(0, 2, 5)) {
  ex <- draw_example_net(h)
  set.seed(11)  # same layout seed across panels for visual comparability
  lay <- layout_with_fr(ex$g)
  plot(ex$g, layout = lay, vertex.color = cols[V(ex$g)$class], vertex.size = 6,
       vertex.label = NA, edge.color = "grey70", edge.width = 0.6,
       main = parse(text = sprintf("theta[H] == %d", h)))
  title(sub = sprintf("realized segregation = %.2f", ex$realized_seg), line = 1, cex.sub = 1.1)
}
dev.off()
cat("Saved sim_example_networks.png\n")

# ---------------------------------------------------------------
# Figure: coefficient plot of the final empirical results
# (composite-index H1 primary test + H3 path decomposition)
# ---------------------------------------------------------------
# All intervals use CR2 Satterthwaite corrections, applied consistently
# across every headline coefficient given the modest number of clusters
# (29): H1's exposure coefficient (df=18.0, t-crit=2.101; this interval
# includes zero) and the path-decomposition lambda/tau rows (df~=17.4-17.9,
# t-crit~=2.11-2.12), computed from the CR2 SEs reported in Table 2.
plot_df <- tribble(
  ~panel,                          ~term,                      ~b,      ~lo,        ~hi,
  "H1: Perceived\ninequality index", "Exposure",   0.00995,  -0.00082,  0.02072,
  "H2: \"Too large\"", "Perceived inequality (lambda)", 0.0747, 0.05269, 0.09671,
  "H2: \"Too large\"", "Exposure, net of perception (tau)", 0.0313, 0.01648, 0.04612,
  "H2: \"Gov't reduce\"", "Perceived inequality (lambda)", 0.0615, 0.03353, 0.08947,
  "H2: \"Gov't reduce\"", "Exposure, net of perception (tau)", 0.0214, 0.00889, 0.03391
) %>%
  mutate(panel = factor(panel, levels = c("H1: Perceived\ninequality index",
                                           "H2: \"Too large\"",
                                           "H2: \"Gov't reduce\"")),
         term = factor(term, levels = c("Exposure", "Perceived inequality (lambda)", "Exposure, net of perception (tau)")))

p <- ggplot(plot_df, aes(x = b, y = term, color = panel)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15, linewidth = 0.9) +
  geom_point(size = 3) +
  facet_wrap(~panel, scales = "free", ncol = 3) +
  scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3"), guide = "none") +
  labs(x = "Coefficient (95% CI)", y = NULL,
       title = "Primary test and path decomposition",
       subtitle = "All intervals CR2-corrected (29 clusters). H2: lambda (perception pathway) and tau (residual direct) from the path decomposition") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9.5, color = "grey30"),
        axis.text.y = element_text(size = 9.5))

ggsave("results/coef_plot_main_results.png",
       p, width = 11, height = 3.6, dpi = 300, bg = "white")
cat("Saved coef_plot_main_results.png\n")

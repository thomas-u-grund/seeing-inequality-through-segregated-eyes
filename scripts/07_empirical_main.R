# Final empirical analysis for the ESR manuscript -- third revision round.
#
# Changes in this round:
# 1. Own household income (within-country percentile, continuous) added to
#    EVERY main model (H1 and both H3 equations), matching the theoretical
#    attitude equation's self-interest term z_i. Tertile/quintile versions
#    reported as robustness.
# 2. CR2 small-sample corrected cluster-robust SEs (clubSandwich) for the
#    headline H1 and H3 results, given only 29 country clusters.
# 3. Cluster-bootstrap CI for the indirect effect (beta*lambda) in the path
#    decomposition, resampling countries.
# 4. Composite-index validation: inter-item correlations, Cronbach's alpha,
#    and a complete-cases-only (all 3 anchors) robustness check.
# 5. Country-specific exposure slopes (H1), for the appendix figure.

suppressPackageStartupMessages({
  library(haven); library(dplyr); library(sandwich); library(lmtest)
  library(clubSandwich); library(boot)
})

data_path <- "raw/ZA7600_v3-0-0.dta"

d <- read_dta(data_path,
              col_select = c(country, c_alphan, v11, v12, v14, v15, v21, v22,
                             v51, v52, WEIGHT, AGE, SEX, EDUCYRS, matches("_INC$")))

clean_earn  <- function(x) { x <- as.numeric(x); ifelse(x <= 0, NA, x) }
clean_scale <- function(x, max_valid) { x <- as.numeric(x); ifelse(x < 1 | x > max_valid, NA, x) }

inc_cols <- grep("_INC$", names(d), value = TRUE)
d <- d %>% mutate(across(all_of(inc_cols), ~ ifelse(as.numeric(.x) < 0, NA, as.numeric(.x))))
d$hh_income <- do.call(dplyr::coalesce, as.list(d[inc_cols]))

d <- d %>%
  mutate(
    doctor_pay = clean_earn(v11), minister_pay = clean_earn(v15), chairman_pay = clean_earn(v12),
    bottom_pay = clean_earn(v14),
    contact_poorer = clean_scale(v51, 7), contact_richer = clean_scale(v52, 7),
    exposure_avg = (contact_poorer + contact_richer) / 2,
    income_too_large  = 6 - clean_scale(v21, 5),
    gov_should_reduce = 6 - clean_scale(v22, 5),
    age = ifelse(as.numeric(AGE) > 0, as.numeric(AGE), NA),
    female = ifelse(as.numeric(SEX) == 2, 1, ifelse(as.numeric(SEX) == 1, 0, NA)),
    educyrs = ifelse(as.numeric(EDUCYRS) >= 0 & as.numeric(EDUCYRS) < 90, as.numeric(EDUCYRS), NA),
    wt = as.numeric(WEIGHT)
  ) %>%
  group_by(c_alphan) %>%
  mutate(income_pct = percent_rank(hh_income)) %>%
  ungroup() %>%
  mutate(income_tertile = case_when(income_pct <= 1/3 ~ 1, income_pct <= 2/3 ~ 2, TRUE ~ 3),
         income_quintile = pmin(5, floor(income_pct * 5) + 1))

trim_by_country <- function(df, var) {
  df %>% group_by(c_alphan) %>%
    mutate(lo = quantile({{ var }}, 0.02, na.rm = TRUE), hi = quantile({{ var }}, 0.98, na.rm = TRUE)) %>%
    ungroup() %>% mutate("{{var}}" := ifelse({{ var }} < lo | {{ var }} > hi, NA, {{ var }})) %>% select(-lo, -hi)
}
d <- trim_by_country(d, bottom_pay); d <- trim_by_country(d, doctor_pay)
d <- trim_by_country(d, minister_pay); d <- trim_by_country(d, chairman_pay)

d <- d %>%
  mutate(
    ratio_doctor   = ifelse(!is.na(doctor_pay)   & !is.na(bottom_pay) & doctor_pay/bottom_pay >= 1   & doctor_pay/bottom_pay <= 300,   log(doctor_pay/bottom_pay), NA),
    ratio_minister = ifelse(!is.na(minister_pay) & !is.na(bottom_pay) & minister_pay/bottom_pay >= 1 & minister_pay/bottom_pay <= 300, log(minister_pay/bottom_pay), NA),
    ratio_chairman = ifelse(!is.na(chairman_pay) & !is.na(bottom_pay) & chairman_pay/bottom_pay >= 1 & chairman_pay/bottom_pay <= 300, log(chairman_pay/bottom_pay), NA)
  ) %>%
  group_by(c_alphan) %>%
  mutate(z_doctor = as.numeric(scale(ratio_doctor)),
         z_minister = as.numeric(scale(ratio_minister)),
         z_chairman = as.numeric(scale(ratio_chairman))) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(composite_index = mean(c(z_doctor, z_minister, z_chairman), na.rm = TRUE),
         n_anchors = sum(!is.na(c(z_doctor, z_minister, z_chairman)))) %>%
  ungroup() %>%
  mutate(composite_index = ifelse(is.nan(composite_index) | n_anchors == 0, NA, composite_index))

# =====================================================================
# Composite-index validation
# =====================================================================
cat("=== Composite index validation ===\n")
item_mat <- d %>% select(z_doctor, z_minister, z_chairman) %>% as.matrix()
cat("Inter-item correlations (pairwise complete):\n")
print(round(cor(item_mat, use = "pairwise.complete.obs"), 3))

cronbach_alpha <- function(mat) {
  mat <- mat[complete.cases(mat), ]
  k <- ncol(mat)
  item_var <- sum(apply(mat, 2, var))
  total_var <- var(rowSums(mat))
  (k / (k - 1)) * (1 - item_var / total_var)
}
cat("Cronbach's alpha (complete cases, n=", sum(complete.cases(item_mat)), "):",
    round(cronbach_alpha(item_mat), 3), "\n\n")

clust_se <- function(model, cluster) coeftest(model, vcov = vcovCL(model, cluster = cluster))
clust_se_cr2 <- function(model, cluster_var_name, data) {
  coef_test(model, vcov = "CR2", cluster = data[[cluster_var_name]], test = "Satterthwaite")
}

# =====================================================================
# H1 (primary): composite index ~ exposure + own income + controls
# =====================================================================
d1 <- d %>% filter(!is.na(composite_index), !is.na(exposure_avg), !is.na(income_pct),
                    !is.na(age), !is.na(female), !is.na(educyrs))
cat("=== H1 (PRIMARY, with own income): N =", nrow(d1), "===\n")
m_h1 <- lm(composite_index ~ exposure_avg + income_pct + age + female + educyrs + factor(c_alphan), data = d1, weights = d1$wt)
cat("Country-clustered (CR1):\n"); print(clust_se(m_h1, d1$c_alphan)[c("exposure_avg","income_pct"), ])
cat("CR2 small-sample corrected:\n")
ct1 <- clust_se_cr2(m_h1, "c_alphan", d1)
print(ct1[rownames(ct1) %in% c("exposure_avg","income_pct"), ])

cat("\n--- Robustness: income as tertile/quintile factor instead of continuous percentile ---\n")
m_h1_tert <- lm(composite_index ~ exposure_avg + factor(income_tertile) + age + female + educyrs + factor(c_alphan), data = d1, weights = d1$wt)
cat("Tertile spec, exposure coef:\n"); print(clust_se(m_h1_tert, d1$c_alphan)["exposure_avg", ])
m_h1_quint <- lm(composite_index ~ exposure_avg + factor(income_quintile) + age + female + educyrs + factor(c_alphan), data = d1, weights = d1$wt)
cat("Quintile spec, exposure coef:\n"); print(clust_se(m_h1_quint, d1$c_alphan)["exposure_avg", ])

cat("\n--- Robustness: complete cases only (all 3 anchors answered) ---\n")
d1_complete <- d1 %>% filter(n_anchors == 3)
m_h1_complete <- lm(composite_index ~ exposure_avg + income_pct + age + female + educyrs + factor(c_alphan), data = d1_complete, weights = d1_complete$wt)
cat("N =", nrow(d1_complete), "\n")
print(clust_se(m_h1_complete, d1_complete$c_alphan)["exposure_avg", ])

# =====================================================================
# H1 robustness: individual anchors, now with income control
# =====================================================================
cat("\n=== H1 robustness: individual anchors (with income control, clustered SEs) ===\n")
for (anchor_ratio in c("ratio_doctor", "ratio_minister", "ratio_chairman")) {
  dd <- d %>% filter(!is.na(.data[[anchor_ratio]]), !is.na(exposure_avg), !is.na(income_pct), !is.na(age), !is.na(female), !is.na(educyrs))
  m <- lm(as.formula(paste0(anchor_ratio, " ~ exposure_avg + income_pct + age + female + educyrs + factor(c_alphan)")), data = dd, weights = dd$wt)
  cat(anchor_ratio, "(N=", nrow(dd), "):\n")
  print(clust_se(m, dd$c_alphan)["exposure_avg", ])
}

# =====================================================================
# H3 path decomposition, with own income in both equations
# =====================================================================
cat("\n\n=== H3 path decomposition (with own income) ===\n")
beta <- coef(m_h1)["exposure_avg"]

boot_results <- list()
for (att_var in c("income_too_large", "gov_should_reduce")) {
  d2 <- d1 %>% filter(!is.na(.data[[att_var]]))
  m_att <- lm(as.formula(paste0(att_var, " ~ composite_index + exposure_avg + income_pct + age + female + educyrs + factor(c_alphan)")),
              data = d2, weights = d2$wt)
  cat("--- Outcome:", att_var, "(N=", nrow(d2), ") ---\n")
  print(clust_se(m_att, d2$c_alphan)[c("composite_index", "exposure_avg", "income_pct"), ])
  cat("CR2 small-sample corrected:\n")
  ct_att <- clust_se_cr2(m_att, "c_alphan", d2)
  print(ct_att[rownames(ct_att) %in% c("composite_index", "exposure_avg", "income_pct"), ])
  lambda <- coef(m_att)["composite_index"]
  tau <- coef(m_att)["exposure_avg"]
  indirect <- beta * lambda
  cat("Indirect (beta*lambda):", round(indirect, 5), " | tau:", round(tau, 5),
      " | indirect share of (tau+indirect):", round(indirect / (indirect + tau) * 100, 1), "%\n\n")

  # cluster bootstrap for CI on beta*lambda (resample countries with replacement)
  countries <- unique(d2$c_alphan)
  set.seed(123)
  nboot <- 5000
  boot_vals <- numeric(nboot)
  for (b in seq_len(nboot)) {
    samp_countries <- sample(countries, length(countries), replace = TRUE)
    boot_df <- bind_rows(lapply(samp_countries, function(cc) d2 %>% filter(c_alphan == cc)))
    m1b <- tryCatch(lm(composite_index ~ exposure_avg + income_pct + age + female + educyrs + factor(c_alphan), data = boot_df, weights = boot_df$wt), error = function(e) NULL)
    m2b <- tryCatch(lm(as.formula(paste0(att_var, " ~ composite_index + exposure_avg + income_pct + age + female + educyrs + factor(c_alphan)")), data = boot_df, weights = boot_df$wt), error = function(e) NULL)
    if (!is.null(m1b) && !is.null(m2b) && "exposure_avg" %in% names(coef(m1b)) && "composite_index" %in% names(coef(m2b))) {
      boot_vals[b] <- coef(m1b)["exposure_avg"] * coef(m2b)["composite_index"]
    } else boot_vals[b] <- NA
  }
  ci <- quantile(boot_vals, c(0.025, 0.975), na.rm = TRUE)
  cat("Cluster-bootstrap 95% CI for indirect effect (", nboot, "reps,", sum(!is.na(boot_vals)), "valid):",
      round(ci[1], 5), "to", round(ci[2], 5), "\n\n")
  boot_results[[att_var]] <- ci
}

# =====================================================================
# Country-specific exposure slopes (H1), for appendix figure
# =====================================================================
cat("\n=== Country-specific H1 slopes (no country FE, within-country OLS) ===\n")
country_slopes <- d1 %>%
  group_by(c_alphan) %>%
  filter(n() >= 200) %>%
  group_modify(~ {
    m <- tryCatch(lm(composite_index ~ exposure_avg + income_pct + age + female + educyrs, data = .x, weights = .x$wt), error = function(e) NULL)
    if (is.null(m) || !("exposure_avg" %in% names(coef(m)))) return(tibble(b = NA, se = NA, n = nrow(.x)))
    s <- summary(m)$coefficients
    tibble(b = s["exposure_avg", 1], se = s["exposure_avg", 2], n = nrow(.x))
  }) %>%
  ungroup() %>%
  filter(!is.na(b))
print(country_slopes, n = 40)
write.csv(country_slopes, "data/country_slopes.csv", row.names = FALSE)

# =====================================================================
# Auxiliary class-conditional channel tests (unchanged design, now also
# with income control for consistency; class_group already conditions on
# income tertile, so this is closely related but reported separately since
# it is testing a different question -- the channel-specific asymmetry)
# =====================================================================
d0 <- d %>%
  mutate(class_group = case_when(income_pct <= 1/3 ~ "low", income_pct >= 2/3 ~ "high", TRUE ~ "mid"),
         class_group = factor(class_group, levels = c("low", "mid", "high")))

cat("\n\n=== Auxiliary class-channel tests (country-clustered SEs, income tertile is the moderator) ===\n")
for (anchor in c("doctor_pay", "minister_pay", "chairman_pay")) {
  dd <- d0 %>% filter(!is.na(.data[[anchor]]), !is.na(bottom_pay), !is.na(class_group), !is.na(income_pct),
                       !is.na(contact_richer), !is.na(contact_poorer), !is.na(age), !is.na(female), !is.na(educyrs))
  dd$log_top <- log(dd[[anchor]]); dd$log_bottom <- log(dd$bottom_pay)
  # class_group is the interaction moderator; income_pct is added alongside it
  # as a separate continuous control, since income varies considerably within
  # tertiles and the moderator alone does not fully absorb it
  m_top <- lm(log_top ~ contact_richer * class_group + income_pct + age + female + educyrs + factor(c_alphan), data = dd, weights = dd$wt)
  m_bot <- lm(log_bottom ~ contact_poorer * class_group + income_pct + age + female + educyrs + factor(c_alphan), data = dd, weights = dd$wt)
  cat("--- Anchor:", anchor, "(N=", nrow(dd), ") ---\n")
  cat("Top channel:\n"); print(clust_se(m_top, dd$c_alphan)[c("contact_richer","contact_richer:class_grouphigh","income_pct"), ])
  cat("Bottom channel:\n"); print(clust_se(m_bot, dd$c_alphan)[c("contact_poorer","contact_poorer:class_grouphigh","income_pct"), ])
}

saveRDS(d, "data/empirical_final_data.rds")
cat("\nDone.\n")

# ============================================================
# HARMONi-2 style simulation: PFS / OS / censoring imbalance
# ============================================================
# Goal:
#   Simulate a two-arm trial resembling the public discussion:
#   - Arm A: ivonescimab-like arm
#   - Arm B: pembrolizumab-like arm
#   - Strong PFS signal (target HR ~0.51)
#   - OS trend only (target HR ~0.77)
#   - Imbalanced early treatment discontinuation / censoring
#
# This is NOT a reconstruction of patient-level truth.
# It is a teaching/demo script to explore how:
#   1) hazard ratio,
#   2) follow-up maturity,
#   3) informative censoring imbalance
# can change apparent trial results.
#
# ============================================================

# -----------------------------
# 0) Packages
# -----------------------------
pkgs <- c("survival", "survminer", "dplyr", "ggplot2", "purrr", "tibble")
new_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs)) install.packages(new_pkgs)

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(purrr)
library(tibble)

set.seed(112)

# -----------------------------
# 1) Helper functions
# -----------------------------

# Piecewise enrollment over 10 months: Nov 2022 -> Aug 2023 style
simulate_enrollment <- function(n, accrual_months = 10) {
  runif(n, min = 0, max = accrual_months)
}

# Generate exponential event times from median
rexp_from_median <- function(n, median_months) {
  rate <- log(2) / median_months
  rexp(n, rate = rate)
}

# Administrative cutoff from first patient in to data cut date
# Example:
#   first patient in = month 0
#   last patient in by month 10
#   data cut at month 15 from FPI
# Then each patient has follow-up = cut_date - entry_time
make_followup_limit <- function(entry_time, cut_date_month) {
  pmax(cut_date_month - entry_time, 0)
}

# Build observed data under event time + censor time + admin cut
build_time_to_event <- function(event_time, censor_time, admin_limit) {
  obs_time <- pmin(event_time, censor_time, admin_limit)
  status   <- as.integer(event_time <= censor_time & event_time <= admin_limit)
  tibble(time = obs_time, status = status)
}

# Kaplan-Meier median helper
safe_km_median <- function(time, status) {
  fit <- survfit(Surv(time, status) ~ 1)
  tb  <- summary(fit)$table
  if (is.null(tb) || is.na(tb["median"])) return(NA_real_)
  as.numeric(tb["median"])
}

# Cox HR helper
safe_cox_hr <- function(time, status, arm) {
  fit <- coxph(Surv(time, status) ~ arm)
  s   <- summary(fit)
  tibble(
    hr       = unname(s$coefficients[1, "exp(coef)"]),
    lower95  = unname(s$conf.int[1, "lower .95"]),
    upper95  = unname(s$conf.int[1, "upper .95"]),
    pvalue   = unname(s$coefficients[1, "Pr(>|z|)"])
  )
}

# Count maturity as event proportion among all randomized patients
calc_maturity <- function(status) {
  mean(status)
}

# Simple one-trial summary
summarize_trial <- function(dat, endpoint_label = "PFS") {
  hr_row <- safe_cox_hr(dat$time, dat$status, dat$arm)

  med <- dat %>%
    group_by(arm) %>%
    summarise(
      n = n(),
      events = sum(status),
      censored = sum(status == 0),
      median = safe_km_median(time, status),
      .groups = "drop"
    )

  out <- med %>%
    mutate(endpoint = endpoint_label) %>%
    bind_cols(hr_row[rep(1, nrow(med)), ])

  out
}

# -----------------------------
# 2) Core simulation engine
# -----------------------------
simulate_trial <- function(
  n_armA = 198,
  n_armB = 200,

  # enrollment / cutoff
  accrual_months = 10,
  pfs_cut_month  = 15,    # roughly short follow-up / WCLC-like
  os_cut_month   = 30,    # longer follow-up / later interim-like

  # target underlying medians
  pfs_median_B = 5.8,
  pfs_hr_A_vs_B = 0.51,
  os_median_B  = 24,
  os_hr_A_vs_B = 0.77,

  # discontinuation/censoring imbalance
  # "early discontinuation" is modeled as a censoring time
  early_disc_A_n = 19,
  early_disc_B_n = 4,
  early_disc_window = c(1, 8),  # months after randomization

  # background non-informative dropout
  dropout_rate = 0.015
) {
  # Derived medians from target HR under exponential survival
  # HR = lambda_A / lambda_B = median_B / median_A
  pfs_median_A <- pfs_median_B / pfs_hr_A_vs_B
  os_median_A  <- os_median_B  / os_hr_A_vs_B

  # Enrollment
  entry_A <- simulate_enrollment(n_armA, accrual_months)
  entry_B <- simulate_enrollment(n_armB, accrual_months)

  # True event times
  true_pfs_A <- rexp_from_median(n_armA, pfs_median_A)
  true_pfs_B <- rexp_from_median(n_armB, pfs_median_B)

  true_os_A  <- rexp_from_median(n_armA, os_median_A)
  true_os_B  <- rexp_from_median(n_armB, os_median_B)

  # Background dropout times
  bg_censor_A <- rexp(n_armA, rate = dropout_rate)
  bg_censor_B <- rexp(n_armB, rate = dropout_rate)

  # Add arm-specific early discontinuation / censoring
  disc_censor_A <- rep(Inf, n_armA)
  disc_censor_B <- rep(Inf, n_armB)

  idxA <- sample(seq_len(n_armA), early_disc_A_n)
  idxB <- sample(seq_len(n_armB), early_disc_B_n)

  disc_censor_A[idxA] <- runif(early_disc_A_n, early_disc_window[1], early_disc_window[2])
  disc_censor_B[idxB] <- runif(early_disc_B_n, early_disc_window[1], early_disc_window[2])

  censor_A <- pmin(bg_censor_A, disc_censor_A)
  censor_B <- pmin(bg_censor_B, disc_censor_B)

  # Administrative follow-up
  admin_pfs_A <- make_followup_limit(entry_A, pfs_cut_month)
  admin_pfs_B <- make_followup_limit(entry_B, pfs_cut_month)
  admin_os_A  <- make_followup_limit(entry_A, os_cut_month)
  admin_os_B  <- make_followup_limit(entry_B, os_cut_month)

  # Observed PFS
  obs_pfs_A <- build_time_to_event(true_pfs_A, censor_A, admin_pfs_A)
  obs_pfs_B <- build_time_to_event(true_pfs_B, censor_B, admin_pfs_B)

  # Observed OS
  # In a standard OS analysis, discontinuing study treatment should not
  # automatically censor OS if death follow-up continues.
  # To illustrate both possibilities, we create:
  #   1) proper OS follow-up (background dropout only)
  #   2) "bad OS handling" where early discontinuation also censors OS
  obs_os_A_proper <- build_time_to_event(true_os_A, bg_censor_A, admin_os_A)
  obs_os_B_proper <- build_time_to_event(true_os_B, bg_censor_B, admin_os_B)

  obs_os_A_bad <- build_time_to_event(true_os_A, censor_A, admin_os_A)
  obs_os_B_bad <- build_time_to_event(true_os_B, censor_B, admin_os_B)

  # Assemble data
  pfs_dat <- bind_rows(
    tibble(
      id = paste0("A_", seq_len(n_armA)),
      arm = "ArmA_ivonescimab_like",
      entry = entry_A,
      true_event = true_pfs_A,
      censor = censor_A,
      time = obs_pfs_A$time,
      status = obs_pfs_A$status,
      endpoint = "PFS",
      early_disc = seq_len(n_armA) %in% idxA
    ),
    tibble(
      id = paste0("B_", seq_len(n_armB)),
      arm = "ArmB_pembro_like",
      entry = entry_B,
      true_event = true_pfs_B,
      censor = censor_B,
      time = obs_pfs_B$time,
      status = obs_pfs_B$status,
      endpoint = "PFS",
      early_disc = seq_len(n_armB) %in% idxB
    )
  )

  os_dat_proper <- bind_rows(
    tibble(
      id = paste0("A_", seq_len(n_armA)),
      arm = "ArmA_ivonescimab_like",
      entry = entry_A,
      true_event = true_os_A,
      censor = bg_censor_A,
      time = obs_os_A_proper$time,
      status = obs_os_A_proper$status,
      endpoint = "OS_proper_followup",
      early_disc = seq_len(n_armA) %in% idxA
    ),
    tibble(
      id = paste0("B_", seq_len(n_armB)),
      arm = "ArmB_pembro_like",
      entry = entry_B,
      true_event = true_os_B,
      censor = bg_censor_B,
      time = obs_os_B_proper$time,
      status = obs_os_B_proper$status,
      endpoint = "OS_proper_followup",
      early_disc = seq_len(n_armB) %in% idxB
    )
  )

  os_dat_bad <- bind_rows(
    tibble(
      id = paste0("A_", seq_len(n_armA)),
      arm = "ArmA_ivonescimab_like",
      entry = entry_A,
      true_event = true_os_A,
      censor = censor_A,
      time = obs_os_A_bad$time,
      status = obs_os_A_bad$status,
      endpoint = "OS_if_early_disc_also_censors",
      early_disc = seq_len(n_armA) %in% idxA
    ),
    tibble(
      id = paste0("B_", seq_len(n_armB)),
      arm = "ArmB_pembro_like",
      entry = entry_B,
      true_event = true_os_B,
      censor = censor_B,
      time = obs_os_B_bad$time,
      status = obs_os_B_bad$status,
      endpoint = "OS_if_early_disc_also_censors",
      early_disc = seq_len(n_armB) %in% idxB
    )
  )

  list(
    pfs = pfs_dat,
    os_proper = os_dat_proper,
    os_bad = os_dat_bad,
    assumptions = list(
      pfs_median_A = pfs_median_A,
      pfs_median_B = pfs_median_B,
      os_median_A  = os_median_A,
      os_median_B  = os_median_B,
      pfs_hr_target = pfs_hr_A_vs_B,
      os_hr_target  = os_hr_A_vs_B
    )
  )
}

# -----------------------------
# 3) Run one simulated trial
# -----------------------------
sim1 <- simulate_trial(
  n_armA = 198,
  n_armB = 200,
  accrual_months = 10,
  pfs_cut_month = 15,    # short follow-up
  os_cut_month  = 30,    # later interim-like
  pfs_median_B = 5.8,
  pfs_hr_A_vs_B = 0.51,
  os_median_B = 24,
  os_hr_A_vs_B = 0.77,
  early_disc_A_n = 19,
  early_disc_B_n = 4,
  early_disc_window = c(1, 8),
  dropout_rate = 0.015
)

# -----------------------------
# 4) Summaries
# -----------------------------
cat("\nAssumptions used:\n")
print(sim1$assumptions)

cat("\nPFS summary:\n")
print(summarize_trial(sim1$pfs, "PFS"))

cat("\nOS summary (proper OS follow-up):\n")
print(summarize_trial(sim1$os_proper, "OS_proper_followup"))

cat("\nOS summary (if early discontinuation also censors OS):\n")
print(summarize_trial(sim1$os_bad, "OS_if_early_disc_also_censors"))

cat("\nMaturity estimates:\n")
cat("PFS maturity =", round(100 * calc_maturity(sim1$pfs$status), 1), "%\n")
cat("OS proper maturity =", round(100 * calc_maturity(sim1$os_proper$status), 1), "%\n")
cat("OS bad maturity =", round(100 * calc_maturity(sim1$os_bad$status), 1), "%\n")

# -----------------------------
# 5) Kaplan-Meier plots
# -----------------------------
fit_pfs <- survfit(Surv(time, status) ~ arm, data = sim1$pfs)
fit_os1 <- survfit(Surv(time, status) ~ arm, data = sim1$os_proper)
fit_os2 <- survfit(Surv(time, status) ~ arm, data = sim1$os_bad)

p1 <- ggsurvplot(
  fit_pfs, data = sim1$pfs,
  risk.table = TRUE, conf.int = FALSE, pval = TRUE,
  title = "Simulated PFS",
  xlab = "Months", ylab = "PFS probability"
)

p2 <- ggsurvplot(
  fit_os1, data = sim1$os_proper,
  risk.table = TRUE, conf.int = FALSE, pval = TRUE,
  title = "Simulated OS (proper follow-up)",
  xlab = "Months", ylab = "OS probability"
)

p3 <- ggsurvplot(
  fit_os2, data = sim1$os_bad,
  risk.table = TRUE, conf.int = FALSE, pval = TRUE,
  title = "Simulated OS (if early discontinuation also censors OS)",
  xlab = "Months", ylab = "OS probability"
)

print(p1)
print(p2)
print(p3)

# -----------------------------
# 6) Show censoring imbalance directly
# -----------------------------
censor_summary <- bind_rows(sim1$pfs, sim1$os_proper, sim1$os_bad) %>%
  group_by(endpoint, arm) %>%
  summarise(
    n = n(),
    events = sum(status),
    censored = sum(status == 0),
    early_disc_n = sum(early_disc),
    early_disc_pct = round(100 * mean(early_disc), 1),
    censor_pct = round(100 * mean(status == 0), 1),
    .groups = "drop"
  )

cat("\nCensoring / early discontinuation summary:\n")
print(censor_summary)

# -----------------------------
# 7) Repeat many trials to study average behavior
# -----------------------------
run_many_trials <- function(
  n_sim = 500,
  n_armA = 198,
  n_armB = 200,
  pfs_cut_month = 15,
  os_cut_month = 30,
  pfs_median_B = 5.8,
  pfs_hr_A_vs_B = 0.51,
  os_median_B = 24,
  os_hr_A_vs_B = 0.77,
  early_disc_A_n = 19,
  early_disc_B_n = 4
) {
  map_dfr(seq_len(n_sim), function(i) {
    s <- simulate_trial(
      n_armA = n_armA,
      n_armB = n_armB,
      pfs_cut_month = pfs_cut_month,
      os_cut_month = os_cut_month,
      pfs_median_B = pfs_median_B,
      pfs_hr_A_vs_B = pfs_hr_A_vs_B,
      os_median_B = os_median_B,
      os_hr_A_vs_B = os_hr_A_vs_B,
      early_disc_A_n = early_disc_A_n,
      early_disc_B_n = early_disc_B_n
    )

    pfs_hr <- safe_cox_hr(s$pfs$time, s$pfs$status, s$pfs$arm)
    os1_hr <- safe_cox_hr(s$os_proper$time, s$os_proper$status, s$os_proper$arm)
    os2_hr <- safe_cox_hr(s$os_bad$time, s$os_bad$status, s$os_bad$arm)

    tibble(
      sim = i,
      pfs_hr = pfs_hr$hr,
      pfs_p  = pfs_hr$pvalue,
      pfs_maturity = calc_maturity(s$pfs$status),

      os_hr_proper = os1_hr$hr,
      os_p_proper  = os1_hr$pvalue,
      os_maturity_proper = calc_maturity(s$os_proper$status),

      os_hr_bad = os2_hr$hr,
      os_p_bad  = os2_hr$pvalue,
      os_maturity_bad = calc_maturity(s$os_bad$status)
    )
  })
}

res_many <- run_many_trials(
  n_sim = 500,
  pfs_cut_month = 15,
  os_cut_month = 30,
  pfs_median_B = 5.8,
  pfs_hr_A_vs_B = 0.51,
  os_median_B = 24,
  os_hr_A_vs_B = 0.77,
  early_disc_A_n = 19,
  early_disc_B_n = 4
)

cat("\nAverage simulation results across 500 runs:\n")
res_many %>%
  summarise(
    mean_pfs_hr = mean(pfs_hr),
    sig_pfs_0.05 = mean(pfs_p < 0.05),
    mean_pfs_maturity = mean(pfs_maturity),

    mean_os_hr_proper = mean(os_hr_proper),
    sig_os_0.05_proper = mean(os_p_proper < 0.05),
    mean_os_maturity_proper = mean(os_maturity_proper),

    mean_os_hr_bad = mean(os_hr_bad),
    sig_os_0.05_bad = mean(os_p_bad < 0.05),
    mean_os_maturity_bad = mean(os_maturity_bad)
  ) %>%
  print()

# -----------------------------
# 8) Visualize distribution of HRs
# -----------------------------
hr_long <- res_many %>%
  select(sim, pfs_hr, os_hr_proper, os_hr_bad) %>%
  tidyr::pivot_longer(-sim, names_to = "endpoint", values_to = "hr")

ggplot(hr_long, aes(x = hr)) +
  geom_histogram(bins = 35) +
  facet_wrap(~ endpoint, scales = "free_y") +
  geom_vline(xintercept = 1, linetype = 2) +
  labs(
    title = "Distribution of estimated hazard ratios across simulations",
    x = "Estimated HR",
    y = "Count"
  )

# -----------------------------
# 9) Sensitivity analysis:
#    What if censoring imbalance disappears?
# -----------------------------
compare_censoring_scenarios <- function(n_sim = 300) {
  scen1 <- run_many_trials(
    n_sim = n_sim,
    early_disc_A_n = 19,
    early_disc_B_n = 4
  ) %>%
    mutate(scenario = "Imbalanced censoring (19 vs 4)")

  scen2 <- run_many_trials(
    n_sim = n_sim,
    early_disc_A_n = 4,
    early_disc_B_n = 4
  ) %>%
    mutate(scenario = "Balanced censoring (4 vs 4)")

  bind_rows(scen1, scen2)
}

sens <- compare_censoring_scenarios(300)

sens_plot_dat <- sens %>%
  select(scenario, pfs_hr, os_hr_proper, os_hr_bad) %>%
  tidyr::pivot_longer(cols = c(pfs_hr, os_hr_proper, os_hr_bad),
                      names_to = "endpoint", values_to = "hr")

ggplot(sens_plot_dat, aes(x = scenario, y = hr)) +
  geom_boxplot() +
  facet_wrap(~ endpoint, scales = "free_y") +
  geom_hline(yintercept = 1, linetype = 2) +
  labs(
    title = "Sensitivity to censoring imbalance",
    x = "",
    y = "Estimated HR"
  ) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

# -----------------------------
# 10) Sensitivity analysis:
#     What if control OS is better than expected?
# -----------------------------
# This reflects the discussion point:
# if the control arm lives longer than assumed, the study may have less power.

compare_control_os <- function(control_os_values = c(18, 21, 24, 27), n_sim = 200) {
  map_dfr(control_os_values, function(medB) {
    run_many_trials(
      n_sim = n_sim,
      os_median_B = medB,
      os_hr_A_vs_B = 0.77,
      early_disc_A_n = 19,
      early_disc_B_n = 4
    ) %>%
      mutate(control_os_median = medB)
  })
}

os_sens <- compare_control_os(c(18, 21, 24, 27), 200)

os_sens %>%
  group_by(control_os_median) %>%
  summarise(
    mean_os_hr = mean(os_hr_proper),
    mean_os_maturity = mean(os_maturity_proper),
    power_like_p_lt_0.05 = mean(os_p_proper < 0.05),
    .groups = "drop"
  ) %>%
  print()

ggplot(os_sens, aes(x = factor(control_os_median), y = os_p_proper < 0.05)) +
  stat_summary(fun = mean, geom = "bar") +
  labs(
    title = "Apparent OS success rate vs control median OS",
    x = "Assumed control median OS (months)",
    y = "Fraction with p < 0.05"
  )

# -----------------------------
# 11) Approximate "39% maturity" calibration
# -----------------------------
# We can search for an OS data cut that yields around 39% OS maturity.

find_cut_for_target_maturity <- function(
  target = 0.39,
  grid = seq(22, 40, by = 0.5),
  n_try = 80,
  os_median_B = 24,
  os_hr_A_vs_B = 0.77
) {
  res <- map_dfr(grid, function(cutm) {
    vals <- replicate(n_try, {
      s <- simulate_trial(
        pfs_cut_month = 15,
        os_cut_month = cutm,
        os_median_B = os_median_B,
        os_hr_A_vs_B = os_hr_A_vs_B,
        early_disc_A_n = 19,
        early_disc_B_n = 4
      )
      calc_maturity(s$os_proper$status)
    })
    tibble(cut_month = cutm, mean_maturity = mean(vals))
  })

  res %>%
    mutate(diff = abs(mean_maturity - target)) %>%
    arrange(diff)
}

maturity_search <- find_cut_for_target_maturity(target = 0.39)
cat("\nClosest cut months to 39% OS maturity:\n")
print(head(maturity_search, 10))

# -----------------------------
# 12) Export one simulated dataset
# -----------------------------
write.csv(sim1$pfs, "simulated_pfs_dataset.csv", row.names = FALSE)
write.csv(sim1$os_proper, "simulated_os_proper_dataset.csv", row.names = FALSE)
write.csv(sim1$os_bad, "simulated_os_bad_dataset.csv", row.names = FALSE)

cat("\nFiles written:\n")
cat("  simulated_pfs_dataset.csv\n")
cat("  simulated_os_proper_dataset.csv\n")
cat("  simulated_os_bad_dataset.csv\n")

# ============================================================
# How to interpret this script
# ============================================================
# 1) PFS:
#    With median_B ~ 5.8 and target HR ~ 0.51, Arm A should look strongly better.
#
# 2) OS:
#    With target HR ~ 0.77 and immature follow-up, you often see a favorable trend
#    without a decisive p-value.
#
# 3) Censoring imbalance:
#    The "19 vs 4" early discontinuation mechanism lets you explore whether
#    unequal censoring can distort observed PFS / OS estimates.
#
# 4) Proper vs improper OS handling:
#    The OS_proper dataset assumes study-treatment discontinuation does NOT
#    automatically stop death follow-up.
#    The OS_bad dataset shows what happens if it does censor OS.
#
# 5) Sensitivity analysis:
#    You can vary:
#      - pfs_hr_A_vs_B
#      - os_hr_A_vs_B
#      - early_disc_A_n / early_disc_B_n
#      - os_median_B
#      - pfs_cut_month / os_cut_month
#
# This is the safest way to discuss the "story" quantitatively without pretending
# to know the real patient-level data.
# ============================================================

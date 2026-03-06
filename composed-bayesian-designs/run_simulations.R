#!/usr/bin/env Rscript
# ============================================================================
# Computational Supplement: Interactions in Composed Bayesian Trial Designs
# Supporting FDA Guidance Comments (Docket No. FDA-2025-D-3217)
# ============================================================================

set.seed(20260307)
options(scipen = 999)

fig_dir <- "figures"
dir.create(fig_dir, showWarnings = FALSE)

# ============================================================================
# HELPERS: Mixture prior posterior computations (Beta-Binomial conjugate)
# ============================================================================

#' Compute posterior mixture weight, posterior probability, and mixture-weighted
#' prior ESS for a two-component Beta mixture prior after observing y in n.
#' ESS = w_post*(a1+b1) + (1-w_post)*(a0+b0): the posterior-weighted average
#' of component prior sample sizes. This is one reasonable proxy; other
#' ESS definitions exist for dynamic borrowing schemes.
#'
#' Prior: w0 * Beta(a1, b1) + (1-w0) * Beta(a0, b0)
#' Posterior: w_post * Beta(a1+y, b1+n-y) + (1-w_post) * Beta(a0+y, b0+n-y)
mixture_posterior <- function(y, n, a1, b1, a0, b0, w0) {
  # Log marginal likelihoods for each component
  log_m1 <- lbeta(a1 + y, b1 + n - y) - lbeta(a1, b1)
  log_m0 <- lbeta(a0 + y, b0 + n - y) - lbeta(a0, b0)

  # Posterior weight on informative component (log-sum-exp for stability)
  log_w1 <- log(w0) + log_m1
  log_w0 <- log(1 - w0) + log_m0
  log_denom <- pmax(log_w1, log_w0) +
    log(exp(log_w1 - pmax(log_w1, log_w0)) +
        exp(log_w0 - pmax(log_w1, log_w0)))
  w_post <- exp(log_w1 - log_denom)

  # ESS from prior components (weighted)
  ess <- w_post * (a1 + b1) + (1 - w_post) * (a0 + b0)

  list(w_post = w_post, ess = ess,
       # Posterior parameters for each component
       a1_post = a1 + y, b1_post = b1 + n - y,
       a0_post = a0 + y, b0_post = b0 + n - y)
}

#' Posterior probability that p > threshold under mixture posterior
post_prob_gt <- function(threshold, post) {
  post$w_post * pbeta(threshold, post$a1_post, post$b1_post, lower.tail = FALSE) +
    (1 - post$w_post) * pbeta(threshold, post$a0_post, post$b0_post, lower.tail = FALSE)
}

# ============================================================================
# DESIGN PARAMETERS
# ============================================================================

# Historical data: 100 patients, 30 responders (rate = 0.30)
n_hist <- 100
y_hist <- 30
p_hist <- y_hist / n_hist  # 0.30

# Informative prior component: Beta(y_hist + 1, n_hist - y_hist + 1) = Beta(31, 71)
a_info <- y_hist + 1  # 31
b_info <- n_hist - y_hist + 1  # 71

# Weakly informative component: Beta(1, 1)
a_weak <- 1
b_weak <- 1

# Initial mixture weight on informative component
w_init <- 0.70

# Null hypothesis: p <= 0.30
p_null <- 0.30

# Sequential design: 3 looks
n_looks <- c(30, 60, 90)
n_max <- max(n_looks)

# Number of simulations
nsim <- 100000

cat("Design ESS from prior:", w_init * (a_info + b_info) + (1 - w_init) * (a_weak + b_weak), "\n")
cat("  Informative component ESS:", a_info + b_info, "\n")
cat("  Weakly informative ESS:", a_weak + b_weak, "\n\n")

# ============================================================================
# SIMULATION 1: Borrowing Weight Trajectories Under Conflict
# ============================================================================
cat("=== Simulation 1: Borrowing Weight Trajectories ===\n")

# Scenarios: true response rates that create varying conflict with historical
p_scenarios <- c(0.30, 0.25, 0.20, 0.15, 0.40, 0.50)
scenario_labels <- c("p=0.30 (no conflict)", "p=0.25 (mild)",
                      "p=0.20 (moderate)", "p=0.15 (severe)",
                      "p=0.40 (mild, opposite)", "p=0.50 (severe, opposite)")

n_scen <- length(p_scenarios)
nsim_traj <- 5000  # fewer sims needed for trajectory visualization

# Store median and IQR of weight at each look for each scenario
weight_summary <- array(NA, dim = c(n_scen, length(n_looks), 3),
                         dimnames = list(NULL, paste0("Look", 1:3),
                                         c("median", "q25", "q75")))
ess_summary <- array(NA, dim = c(n_scen, length(n_looks), 3),
                      dimnames = list(NULL, paste0("Look", 1:3),
                                      c("median", "q25", "q75")))

for (s in seq_along(p_scenarios)) {
  p_true <- p_scenarios[s]
  weight_mat <- matrix(NA, nsim_traj, length(n_looks))
  ess_mat <- matrix(NA, nsim_traj, length(n_looks))

  # Pre-generate patient-level data for coherent paths
  patient_outcomes <- matrix(rbinom(nsim_traj * n_max, 1, p_true),
                             nrow = nsim_traj, ncol = n_max)

  for (k in seq_along(n_looks)) {
    y_cum <- rowSums(patient_outcomes[, 1:n_looks[k], drop = FALSE])
    post <- mixture_posterior(y_cum, n_looks[k], a_info, b_info, a_weak, b_weak, w_init)
    weight_mat[, k] <- post$w_post
    ess_mat[, k] <- post$ess
  }

  for (k in seq_along(n_looks)) {
    weight_summary[s, k, ] <- quantile(weight_mat[, k], c(0.5, 0.25, 0.75))
    ess_summary[s, k, ] <- quantile(ess_mat[, k], c(0.5, 0.25, 0.75))
  }
}

# --- Figure 1: Borrowing Weight Trajectories ---
png(file.path(fig_dir, "fig1_weight_trajectories.png"), width = 900, height = 550, res = 120)
par(mar = c(5.5, 4.5, 2.5, 1))

cols <- c("#2166AC", "#4393C3", "#D6604D", "#B2182B", "#92C5DE", "#67001F")
plot(NULL, xlim = c(0.8, 3.2), ylim = c(0, 1),
     xlab = "", ylab = "Posterior Weight on Informative Component",
     main = "Borrowing Weight Across Sequential Looks",
     xaxt = "n", las = 1, cex.lab = 1.1)
axis(1, at = 1:3, labels = paste0("Look ", 1:3, "\n(n=", n_looks, ")"), padj = 0.5)
mtext("Interim Look", side = 1, line = 4.2, cex = 1.1)
abline(h = w_init, lty = 2, col = "gray50")
text(0.8, w_init + 0.03, expression(w[0] == 0.70), col = "gray50", cex = 0.8, adj = 0)

for (s in seq_along(p_scenarios)) {
  lines(1:3, weight_summary[s, , "median"], col = cols[s], lwd = 2.5, type = "b", pch = 16)
  for (k in 1:3) {
    segments(k, weight_summary[s, k, "q25"], k, weight_summary[s, k, "q75"],
             col = adjustcolor(cols[s], 0.4), lwd = 4)
  }
}

legend(1.2, 0.35, legend = scenario_labels[c(4, 3, 2, 1, 5, 6)],
       col = cols[c(4, 3, 2, 1, 5, 6)], lwd = 2, pch = 16, cex = 0.7, bty = "n")
dev.off()
cat("  Figure 1 saved.\n")

# --- Figure 2: Effective ESS Trajectories ---
png(file.path(fig_dir, "fig2_ess_trajectories.png"), width = 900, height = 550, res = 120)
par(mar = c(5.5, 4.5, 2.5, 1))

plot(NULL, xlim = c(0.8, 3.2), ylim = c(0, 100),
     xlab = "", ylab = "Mixture-Weighted Prior ESS",
     main = "Prior ESS Across Sequential Looks",
     xaxt = "n", las = 1, cex.lab = 1.1)
axis(1, at = 1:3, labels = paste0("Look ", 1:3, "\n(n=", n_looks, ")"), padj = 0.5)
mtext("Interim Look", side = 1, line = 4.2, cex = 1.1)
design_ess <- w_init * (a_info + b_info) + (1 - w_init) * (a_weak + b_weak)
abline(h = design_ess, lty = 2, col = "gray50")
text(1.1, design_ess + 2.5, paste0("Design ESS = ", round(design_ess, 1)),
     col = "gray50", cex = 0.8, adj = 0)

for (s in seq_along(p_scenarios)) {
  lines(1:3, ess_summary[s, , "median"], col = cols[s], lwd = 2.5, type = "b", pch = 16)
  for (k in 1:3) {
    segments(k, ess_summary[s, k, "q25"], k, ess_summary[s, k, "q75"],
             col = adjustcolor(cols[s], 0.4), lwd = 4)
  }
}

legend(1.2, 35, legend = scenario_labels[c(4, 3, 2, 1, 5, 6)],
       col = cols[c(4, 3, 2, 1, 5, 6)], lwd = 2, pch = 16, cex = 0.7, bty = "n")
dev.off()
cat("  Figure 2 saved.\n")


# ============================================================================
# SIMULATION 2: Calibration Gap — Component vs. Composed
# ============================================================================
cat("\n=== Simulation 2: Calibration Gap ===\n")

#' Simulate sequential trial with dynamic mixture prior
#' Returns: vector of reject (TRUE/FALSE) for each simulation
simulate_sequential <- function(nsim, p_true, n_looks, threshold,
                                 a1, b1, a0, b0, w0, p_null,
                                 dynamic = TRUE) {
  rejected <- logical(nsim)
  stopped <- logical(nsim)

  # Pre-generate patient-level outcomes for coherent path accumulation
  n_max <- max(n_looks)
  patient_outcomes <- matrix(rbinom(nsim * n_max, 1, p_true), nrow = nsim, ncol = n_max)

  y_cum_all <- matrix(NA_real_, nrow = nsim, ncol = length(n_looks))
  for (k in seq_along(n_looks)) {
    y_cum_all[, k] <- rowSums(patient_outcomes[, 1:n_looks[k], drop = FALSE])
  }

  for (k in seq_along(n_looks)) {
    active <- !stopped
    n_active <- sum(active)
    if (n_active == 0) break

    y_cum <- y_cum_all[active, k]

    if (dynamic) {
      post <- mixture_posterior(y_cum, n_looks[k], a1, b1, a0, b0, w0)
    } else {
      post <- list(
        w_post = rep(w0, n_active),
        a1_post = a1 + y_cum, b1_post = b1 + n_looks[k] - y_cum,
        a0_post = a0 + y_cum, b0_post = b0 + n_looks[k] - y_cum
      )
    }

    pp <- post_prob_gt(p_null, post)
    declare <- pp > threshold

    rejected[which(active)[declare]] <- TRUE
    stopped[which(active)[declare]] <- TRUE
  }
  rejected
}

# Two-stage calibration: coarse grid, then fine grid around best candidate
calibrate_threshold <- function(nsim, p_null, n_looks, a1, b1, a0, b0, w0,
                                 dynamic, alpha_target = 0.025) {
  # Stage 1: coarse grid
  grid1 <- seq(0.955, 0.999, by = 0.001)
  alpha1 <- sapply(grid1, function(c) {
    mean(simulate_sequential(nsim, p_null, n_looks, c,
      a1, b1, a0, b0, w0, p_null, dynamic = dynamic))
  })
  best1 <- grid1[which.min(abs(alpha1 - alpha_target))]

  # Stage 2: fine grid around coarse best
  grid2 <- seq(max(best1 - 0.002, 0.950), min(best1 + 0.002, 0.9999), by = 0.0002)
  alpha2 <- sapply(grid2, function(c) {
    mean(simulate_sequential(nsim, p_null, n_looks, c,
      a1, b1, a0, b0, w0, p_null, dynamic = dynamic))
  })
  best2 <- grid2[which.min(abs(alpha2 - alpha_target))]

  # Stage 3: independent verification on fresh samples to avoid selection bias
  verified <- mean(simulate_sequential(nsim, p_null, n_looks, best2,
    a1, b1, a0, b0, w0, p_null, dynamic = dynamic))
  mc_se <- sqrt(verified * (1 - verified) / nsim)

  list(threshold = best2, alpha = verified, mc_se = mc_se)
}

alpha_target <- 0.025

cat("  Calibrating threshold (dynamic, two-stage)...\n")
calib_dynamic <- calibrate_threshold(nsim, p_null, n_looks,
  a_info, b_info, a_weak, b_weak, w_init, dynamic = TRUE)
best_dynamic <- calib_dynamic$threshold
cat(sprintf("  Dynamic-calibrated threshold: %.4f\n", best_dynamic))
cat(sprintf("  Verified Type I error: %.4f ± %.4f\n", calib_dynamic$alpha, calib_dynamic$mc_se))

cat("  Calibrating threshold (static, two-stage)...\n")
calib_static <- calibrate_threshold(nsim, p_null, n_looks,
  a_info, b_info, a_weak, b_weak, w_init, dynamic = FALSE)
best_static <- calib_static$threshold
cat(sprintf("  Static-calibrated threshold: %.4f\n", best_static))
cat(sprintf("  Verified Type I error: %.4f ± %.4f\n", calib_static$alpha, calib_static$mc_se))

# Step 3: Cross-evaluate — use static-calibrated threshold with dynamic prior
# This is what happens when a sponsor calibrates assuming fixed borrowing
# but the actual trial uses dynamic borrowing
cat("\n  Cross-evaluation: static threshold applied to dynamic trial\n")

p_grid <- seq(0.10, 0.55, by = 0.025)
reject_dynamic_dynamic <- numeric(length(p_grid))  # correct: dynamic threshold, dynamic trial
reject_static_dynamic <- numeric(length(p_grid))   # mismatch: static threshold, dynamic trial
reject_static_static <- numeric(length(p_grid))     # component-level: static throughout

for (i in seq_along(p_grid)) {
  rej_dd <- simulate_sequential(nsim, p_grid[i], n_looks, best_dynamic,
                                 a_info, b_info, a_weak, b_weak, w_init, p_null,
                                 dynamic = TRUE)
  reject_dynamic_dynamic[i] <- mean(rej_dd)

  rej_sd <- simulate_sequential(nsim, p_grid[i], n_looks, best_static,
                                 a_info, b_info, a_weak, b_weak, w_init, p_null,
                                 dynamic = TRUE)
  reject_static_dynamic[i] <- mean(rej_sd)

  rej_ss <- simulate_sequential(nsim, p_grid[i], n_looks, best_static,
                                 a_info, b_info, a_weak, b_weak, w_init, p_null,
                                 dynamic = FALSE)
  reject_static_static[i] <- mean(rej_ss)

  if (i %% 5 == 0) cat("    p =", p_grid[i], "done\n")
}

# --- Figure 3: Operating Characteristic Curves ---
png(file.path(fig_dir, "fig3_oc_curves.png"), width = 900, height = 550, res = 120)
par(mar = c(4.5, 4.5, 3, 1))

plot(p_grid, reject_static_static, type = "l", lwd = 2.5, col = "#4393C3",
     ylim = c(0, 1), xlab = "True Response Rate",
     ylab = "Probability of Declaring Efficacy",
     main = "Operating Characteristics: Component-Level vs. Composed Calibration",
     las = 1, cex.lab = 1.05)
lines(p_grid, reject_static_dynamic, lwd = 2.5, col = "#D6604D", lty = 2)
lines(p_grid, reject_dynamic_dynamic, lwd = 2.5, col = "#2166AC", lty = 1)

abline(v = p_null, lty = 3, col = "gray50")
abline(h = alpha_target, lty = 3, col = "gray50")
text(p_null + 0.01, 0.08, expression(H[0]: ~ p == 0.30), cex = 0.8, adj = 0)
text(0.12, alpha_target + 0.03, expression(alpha == 0.025), cex = 0.8, col = "gray50")

legend("topleft", bty = "n", cex = 0.8,
       legend = c(
         paste0("Composed: dynamic prior + dynamic threshold (c=", formatC(best_dynamic, format = "f", digits = 4), ")"),
         paste0("Mismatch: dynamic prior + static threshold (c=", formatC(best_static, format = "f", digits = 4), ")"),
         paste0("Component: static prior + static threshold (c=", formatC(best_static, format = "f", digits = 4), ")")
       ),
       col = c("#2166AC", "#D6604D", "#4393C3"),
       lwd = 2.5, lty = c(1, 2, 1))
dev.off()
cat("  Figure 3 saved.\n")

# Report key Type I error values
cat("\n  Type I error at p =", p_null, ":\n")
idx_null <- which.min(abs(p_grid - p_null))
cat("    Composed (dynamic/dynamic):", reject_dynamic_dynamic[idx_null], "\n")
cat("    Mismatch (static thresh/dynamic prior):", reject_static_dynamic[idx_null], "\n")
cat("    Component (static/static):", reject_static_static[idx_null], "\n")


# ============================================================================
# SIMULATION 3: ESS Distribution and Power Impact Under Drift
# ============================================================================
cat("\n=== Simulation 3: ESS Distribution Under Drift ===\n")

# For different true rates, show the distribution of realized ESS at final analysis
drift_scenarios <- c(0.30, 0.25, 0.20, 0.15)
drift_labels <- c("No drift (p=0.30)", "Mild drift (p=0.25)",
                   "Moderate drift (p=0.20)", "Severe drift (p=0.15)")

nsim_ess <- 50000
ess_at_final <- matrix(NA, nsim_ess, length(drift_scenarios))

for (s in seq_along(drift_scenarios)) {
  y_final <- rbinom(nsim_ess, n_max, drift_scenarios[s])
  post <- mixture_posterior(y_final, n_max, a_info, b_info, a_weak, b_weak, w_init)
  ess_at_final[, s] <- post$ess
}

# --- Figure 4: ESS Distributions ---
png(file.path(fig_dir, "fig4_ess_distributions.png"), width = 900, height = 500, res = 120)
par(mar = c(4.5, 4.5, 3, 1))

ess_cols <- c("#2166AC", "#4393C3", "#D6604D", "#B2182B")
plot(NULL, xlim = c(0, 110), ylim = c(0, 0.20),
     xlab = "Effective Sample Size (Prior Contribution)",
     ylab = "Density",
     main = "Distribution of Realized ESS at Final Analysis (n=90)",
     las = 1, cex.lab = 1.05)

for (s in seq_along(drift_scenarios)) {
  d <- density(ess_at_final[, s], from = 0, to = 110, bw = 2)
  polygon(d$x, d$y, col = adjustcolor(ess_cols[s], 0.2), border = ess_cols[s], lwd = 2)
}

abline(v = design_ess, lty = 2, col = "gray30", lwd = 1.5)
text(design_ess + 1.5, 0.19, paste0("Design ESS = ", round(design_ess, 0)),
     cex = 0.8, adj = 0, col = "gray30")

legend(1, 0.2, bty = "n", cex = 0.8,
       legend = drift_labels, fill = adjustcolor(ess_cols, 0.3),
       border = ess_cols)

# Add summary stats
for (s in seq_along(drift_scenarios)) {
  med_ess <- median(ess_at_final[, s])
  cat("  ", drift_labels[s], ": median ESS =", round(med_ess, 1),
      ", IQR = [", round(quantile(ess_at_final[, s], 0.25), 1), ",",
      round(quantile(ess_at_final[, s], 0.75), 1), "]\n")
}
dev.off()
cat("  Figure 4 saved.\n")


# ============================================================================
# SIMULATION 4: Power Impact — Designed vs. Actual
# ============================================================================
cat("\n=== Simulation 4: Power Under Drift Scenarios ===\n")

# The sponsor designs the trial assuming no drift (historical control rate = 0.30)
# and calculates power for detecting a treatment effect (alternative p = 0.45)
# But what if the null landscape has drifted?

# We evaluate power at p_alternative across different "true null" baselines
# (i.e., the population's actual base rate has drifted from the historical 0.30)

p_alt_grid <- seq(0.30, 0.60, by = 0.02)
drift_nulls <- c(0.30, 0.25, 0.20)
drift_null_labels <- c("Null=0.30 (designed)", "Null=0.25 (mild drift)", "Null=0.20 (moderate drift)")

power_dynamic <- matrix(NA, length(p_alt_grid), length(drift_nulls))
power_static <- matrix(NA, length(p_alt_grid), length(drift_nulls))

for (d in seq_along(drift_nulls)) {
  for (i in seq_along(p_alt_grid)) {
    # Dynamic borrowing with composed threshold
    rej_dyn <- simulate_sequential(nsim, p_alt_grid[i], n_looks, best_dynamic,
                                    a_info, b_info, a_weak, b_weak, w_init, drift_nulls[d],
                                    dynamic = TRUE)
    power_dynamic[i, d] <- mean(rej_dyn)

    # Static borrowing (component-level calibration)
    rej_stat <- simulate_sequential(nsim, p_alt_grid[i], n_looks, best_static,
                                     a_info, b_info, a_weak, b_weak, w_init, drift_nulls[d],
                                     dynamic = FALSE)
    power_static[i, d] <- mean(rej_stat)
  }
  cat("  Drift null =", drift_nulls[d], "done\n")
}

# --- Figure 5: Power Curves Under Drift ---
png(file.path(fig_dir, "fig5_power_curves.png"), width = 900, height = 550, res = 120)
par(mar = c(4.5, 4.5, 3, 1))

drift_cols <- c("#2166AC", "#D6604D", "#B2182B")
plot(NULL, xlim = range(p_alt_grid), ylim = c(0, 1),
     xlab = "True Response Rate", ylab = "Power (Probability of Declaring Efficacy)",
     main = "Power: Composed (Dynamic) vs. Component (Static) Under Null Drift",
     las = 1, cex.lab = 1.05)

# Draw solid lines (dynamic) first as background
for (d in seq_along(drift_nulls)) {
  lines(p_alt_grid, power_dynamic[, d], col = drift_cols[d], lwd = 2.5)
}

# Draw dashed lines (static) on top so they are visible
for (d in seq_along(drift_nulls)) {
  lines(p_alt_grid, power_static[, d], col = drift_cols[d], lwd = 2.5, lty = 3)
}

abline(h = 0.80, lty = 3, col = "gray60")
text(0.31, 0.83, "80% power", cex = 0.75, col = "gray50")

legend(0.5, 0.35, bty = "n", cex = 0.75,
       legend = c(drift_null_labels, "",
                  paste0(drift_null_labels, " (static)")),
       col = c(drift_cols, NA, drift_cols),
       lwd = c(rep(2.5, 3), NA, rep(2.5, 3)),
       lty = c(rep(1, 3), NA, rep(3, 3)),
       title = "Solid = dynamic, dotted = static")
dev.off()
cat("  Figure 5 saved.\n")


# ============================================================================
# SIMULATION 5: Alpha Spent at Each Look — The Core Interaction
# ============================================================================
cat("\n=== Simulation 5: Alpha Spending by Look ===\n")

# Show how much alpha is actually spent at each interim under dynamic vs static
# This is the sharpest demonstration of the interaction

compute_alpha_by_look <- function(nsim, p_true, n_looks, threshold,
                                   a1, b1, a0, b0, w0, p_null, dynamic) {
  alpha_per_look <- numeric(length(n_looks))
  stopped <- logical(nsim)

  # Pre-generate coherent patient-level data
  n_max <- max(n_looks)
  patient_outcomes <- matrix(rbinom(nsim * n_max, 1, p_true), nrow = nsim, ncol = n_max)
  y_cum_all <- matrix(NA_real_, nrow = nsim, ncol = length(n_looks))
  for (k in seq_along(n_looks)) {
    y_cum_all[, k] <- rowSums(patient_outcomes[, 1:n_looks[k], drop = FALSE])
  }

  for (k in seq_along(n_looks)) {
    active <- !stopped
    n_active <- sum(active)
    if (n_active == 0) break

    y_cum <- y_cum_all[active, k]

    if (dynamic) {
      post <- mixture_posterior(y_cum, n_looks[k], a1, b1, a0, b0, w0)
    } else {
      post <- list(
        w_post = rep(w0, n_active),
        a1_post = a1 + y_cum, b1_post = b1 + n_looks[k] - y_cum,
        a0_post = a0 + y_cum, b0_post = b0 + n_looks[k] - y_cum
      )
    }

    pp <- post_prob_gt(p_null, post)
    declare <- pp > threshold

    alpha_per_look[k] <- sum(declare) / nsim
    stopped[which(active)[declare]] <- TRUE
  }

  alpha_per_look
}

alpha_dynamic <- compute_alpha_by_look(nsim, p_null, n_looks, best_dynamic,
                                        a_info, b_info, a_weak, b_weak, w_init, p_null,
                                        dynamic = TRUE)

alpha_static <- compute_alpha_by_look(nsim, p_null, n_looks, best_static,
                                       a_info, b_info, a_weak, b_weak, w_init, p_null,
                                       dynamic = FALSE)

# What happens when you apply the static threshold to the dynamic trial
alpha_mismatch <- compute_alpha_by_look(nsim, p_null, n_looks, best_static,
                                         a_info, b_info, a_weak, b_weak, w_init, p_null,
                                         dynamic = TRUE)

cat("  Alpha by look (dynamic threshold + dynamic prior):", round(alpha_dynamic, 5), "| total:", round(sum(alpha_dynamic), 5), "\n")
cat("  Alpha by look (static threshold + static prior):", round(alpha_static, 5), "| total:", round(sum(alpha_static), 5), "\n")
cat("  Alpha by look (static threshold + dynamic prior):", round(alpha_mismatch, 5), "| total:", round(sum(alpha_mismatch), 5), "\n")

# --- Figure 6: Alpha Spending by Look ---
png(file.path(fig_dir, "fig6_alpha_spending.png"), width = 900, height = 550, res = 120)
par(mar = c(6.5, 4.5, 3, 1))

bar_data <- rbind(alpha_dynamic, alpha_mismatch, alpha_static)
bar_cols <- c("#2166AC", "#D6604D", "#4393C3")

bp <- barplot(bar_data, beside = TRUE, col = bar_cols,
        names.arg = paste0("Look ", 1:3, " (n=", n_looks, ")"),
        ylab = expression("Type I Error (" * alpha * " spent)"),
        main = expression("Alpha Spending by Interim Look Under " * H[0] * ": p = 0.30"),
        las = 1, ylim = c(0, max(bar_data) * 1.4), cex.lab = 1.05,
        border = NA)

# Add value labels
for (i in 1:nrow(bar_data)) {
  text(bp[i, ], bar_data[i, ] + max(bar_data) * 0.04,
       formatC(bar_data[i, ], format = "f", digits = 4), cex = 0.7, col = bar_cols[i])
}

# Legend near Look 1 (topleft)
legend("topleft", bty = "n", cex = 0.7,
       legend = c(
         paste0("Composed: dynamic prior + dynamic threshold (c=", formatC(best_dynamic, format = "f", digits = 4), ")"),
         paste0("Mismatch: dynamic prior + static threshold (c=", formatC(best_static, format = "f", digits = 4), ")"),
         paste0("Component: static prior + static threshold (c=", formatC(best_static, format = "f", digits = 4), ")")
       ),
       fill = bar_cols, border = NA)

# Footer totals on separate lines
mtext(paste0("Totals \u2014 Composed: ", formatC(sum(alpha_dynamic), format = "f", digits = 4),
             "    Mismatch: ", formatC(sum(alpha_mismatch), format = "f", digits = 4),
             "    Component: ", formatC(sum(alpha_static), format = "f", digits = 4)),
      side = 1, line = 3, cex = 0.75)
dev.off()
cat("  Figure 6 saved.\n")

# ============================================================================
# SUMMARY TABLE
# ============================================================================
cat("\n=== Summary Table ===\n")
cat(sprintf("%-45s  %10s  %10s  %10s\n", "Metric", "Composed", "Mismatch", "Component"))
cat(paste(rep("-", 80), collapse = ""), "\n")
cat(sprintf("%-45s  %10.4f  %10.4f  %10.4f\n",
            "Threshold (c)",
            best_dynamic, best_static, best_static))
cat(sprintf("%-45s  %10.4f  %10.4f  %10.4f\n",
            "Type I error (overall)",
            sum(alpha_dynamic), sum(alpha_mismatch), sum(alpha_static)))
cat(sprintf("%-45s  %10.4f  %10.4f  %10.4f\n",
            "  at Look 1 (n=30)",
            alpha_dynamic[1], alpha_mismatch[1], alpha_static[1]))
cat(sprintf("%-45s  %10.4f  %10.4f  %10.4f\n",
            "  at Look 2 (n=60)",
            alpha_dynamic[2], alpha_mismatch[2], alpha_static[2]))
cat(sprintf("%-45s  %10.4f  %10.4f  %10.4f\n",
            "  at Look 3 (n=90)",
            alpha_dynamic[3], alpha_mismatch[3], alpha_static[3]))

# Power at p=0.45 (designed alternative)
idx_45 <- which.min(abs(p_alt_grid - 0.45))
cat(sprintf("%-45s  %10.4f  %10s  %10.4f\n",
            "Power at p=0.45, null=0.30",
            power_dynamic[idx_45, 1], "  -  ", power_static[idx_45, 1]))
idx_45_drift <- which.min(abs(p_alt_grid - 0.45))
cat(sprintf("%-45s  %10.4f  %10s  %10.4f\n",
            "Power at p=0.45, null=0.20 (drift)",
            power_dynamic[idx_45_drift, 3], "  -  ", power_static[idx_45_drift, 3]))

cat("\n=== All figures saved to", fig_dir, "===\n")

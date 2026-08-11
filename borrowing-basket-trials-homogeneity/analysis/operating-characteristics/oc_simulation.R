# =====================================================================
# oc_simulation.R
# Reimplements the operating-characteristics story of Zhou & Ji 2024
# (Tables 5 & 6): how borrowing trades power against type I error across
# homogeneous and mixed basket scenarios.
#
# Thresholds are estimated from 2,500 calibration trials per scenario and
# evaluated on 5,000 separate trials per scenario. Wilson 95% Monte Carlo
# intervals accompany every reported probability. The shared simulated
# counts are saved so the Python validation uses identical inputs.
#
# Base R only. Run: Rscript oc_simulation.R
# =====================================================================

expit <- plogis
logit <- qlogis
l0 <- logit(0.20)
N <- 20L
NB <- 4L
R_CAL <- 2500L
R_EVAL <- 5000L
N_ITER <- 8000L
BURN <- 3000L

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) {
  script_path <- gsub("~\\+~", " ", sub("^--file=", "", file_arg[1]))
  dirname(normalizePath(script_path))
} else {
  getwd()
}

## True response rates per scenario (reference 20% everywhere)
SCEN <- list(
  "1" = c(.20, .20, .20, .20),
  "2" = c(.35, .35, .35, .35),
  "3" = c(.20, .35, .35, .35),
  "4" = c(.20, .20, .35, .35),
  "5" = c(.10, .20, .30, .40),
  "6" = c(.20, .20, .20, .35)
)
promising <- lapply(SCEN, function(p) which(p > 0.20))
nullset <- lapply(SCEN, function(p) which(p <= 0.20))
A_SIGMA <- c(Moderate = 3.0, Strong = 0.3)

binom_ll <- function(y, gamma) {
  p <- pmin(pmax(expit(l0 + gamma), 1e-12), 1 - 1e-12)
  y * log(p) + (N - y) * log(1 - p)
}

## ---- Prior I: adaptive integration over the full real line --------
pp_noborrow <- function() {
  vapply(0:N, function(yy) {
    kernel <- function(gamma) {
      p <- expit(l0 + gamma)
      dnorm(gamma, 0, 100) * p^yy * (1 - p)^(N - yy)
    }
    numerator <- integrate(kernel, 0, Inf, rel.tol = 1e-10,
                           subdivisions = 2000L)$value
    denominator <- integrate(kernel, -Inf, Inf, rel.tol = 1e-10,
                             subdivisions = 2000L)$value
    numerator / denominator
  }, numeric(1))
}
PP_NB <- pp_noborrow()

## ---- Priors II/III: vectorized Metropolis-within-Gibbs ------------
run_hier <- function(y, A_sigma, n_iter = N_ITER, burn = BURN) {
  nr <- nrow(y)
  gamma <- matrix(0, nr, NB)
  mu <- rep(0, nr)
  sigma <- rep(0.5, nr)
  ll <- binom_ll(y, gamma)
  step_gamma <- 0.5
  step_log_sigma <- 0.35
  posterior_positive <- matrix(0, nr, NB)
  gamma_accept <- 0
  sigma_accept <- 0
  kept <- 0L

  for (it in seq_len(n_iter)) {
    proposal <- gamma + matrix(rnorm(nr * NB, 0, step_gamma), nr, NB)
    ll_proposal <- binom_ll(y, proposal)
    lp_old <- -0.5 * ((gamma - mu) / sigma)^2
    lp_new <- -0.5 * ((proposal - mu) / sigma)^2
    accept <- matrix(log(runif(nr * NB)), nr, NB) <
      (ll_proposal + lp_new) - (ll + lp_old)
    gamma_accept <- gamma_accept + sum(accept)
    gamma[accept] <- proposal[accept]
    ll[accept] <- ll_proposal[accept]

    precision <- 1 / 1e4 + NB / sigma^2
    mu <- rnorm(nr, (rowSums(gamma) / sigma^2) / precision,
                sqrt(1 / precision))

    log_sigma <- log(sigma)
    proposal_log_sigma <- log_sigma + rnorm(nr, 0, step_log_sigma)
    proposal_sigma <- exp(proposal_log_sigma)
    log_target <- function(s, log_s) {
      -NB * log(s) - 0.5 * rowSums((gamma - mu)^2) / s^2 -
        s^2 / (2 * A_sigma^2) + log_s
    }
    accept_sigma <- log(runif(nr)) <
      log_target(proposal_sigma, proposal_log_sigma) -
      log_target(sigma, log_sigma)
    sigma_accept <- sigma_accept + sum(accept_sigma)
    sigma[accept_sigma] <- proposal_sigma[accept_sigma]

    if (it > burn) {
      kept <- kept + 1L
      posterior_positive <- posterior_positive + (gamma > 0)
    }
  }

  list(
    probability = posterior_positive / kept,
    gamma_acceptance = gamma_accept / (n_iter * nr * NB),
    sigma_acceptance = sigma_accept / (n_iter * nr)
  )
}

## ---- Independent calibration and evaluation data ------------------
simulate_counts <- function(n_trials, probabilities, seed) {
  set.seed(seed)
  matrix(
    rbinom(n_trials * NB, N, rep(probabilities, each = n_trials)),
    n_trials,
    NB
  )
}

cal_y <- list()
eval_y <- list()
input_rows <- list()
for (s in names(SCEN)) {
  cal_y[[s]] <- simulate_counts(R_CAL, SCEN[[s]], 11000L + as.integer(s))
  eval_y[[s]] <- simulate_counts(R_EVAL, SCEN[[s]], 21000L + as.integer(s))
  input_rows[[paste0("cal.", s)]] <- data.frame(
    phase = "calibration", scenario = as.integer(s), replicate = seq_len(R_CAL),
    b1 = cal_y[[s]][, 1], b2 = cal_y[[s]][, 2],
    b3 = cal_y[[s]][, 3], b4 = cal_y[[s]][, 4]
  )
  input_rows[[paste0("eval.", s)]] <- data.frame(
    phase = "evaluation", scenario = as.integer(s), replicate = seq_len(R_EVAL),
    b1 = eval_y[[s]][, 1], b2 = eval_y[[s]][, 2],
    b3 = eval_y[[s]][, 3], b4 = eval_y[[s]][, 4]
  )
}
write.csv(do.call(rbind, input_rows),
          file.path(script_dir, "oc_simulation_inputs.csv"), row.names = FALSE)

## Posterior Pr(gamma_j > 0) is computed once, then split by phase.
prob_cal <- list()
prob_eval <- list()
diagnostic_rows <- list()
for (s in names(SCEN)) {
  y_all <- rbind(cal_y[[s]], eval_y[[s]])
  nr_all <- nrow(y_all)
  no_borrow <- matrix(PP_NB[as.integer(y_all) + 1L], nr_all, NB)
  prob_cal[[paste0(s, ".No")]] <- no_borrow[seq_len(R_CAL), , drop = FALSE]
  prob_eval[[paste0(s, ".No")]] <- no_borrow[R_CAL + seq_len(R_EVAL), , drop = FALSE]

  for (lab in c("Moderate", "Strong")) {
    set.seed(if (lab == "Moderate") 31000L + as.integer(s) else 41000L + as.integer(s))
    fit <- run_hier(y_all, A_SIGMA[lab])
    prob_cal[[paste0(s, ".", lab)]] <-
      fit$probability[seq_len(R_CAL), , drop = FALSE]
    prob_eval[[paste0(s, ".", lab)]] <-
      fit$probability[R_CAL + seq_len(R_EVAL), , drop = FALSE]
    diagnostic_rows[[paste(s, lab, sep = ".")]] <- data.frame(
      scenario = as.integer(s), borrow = lab,
      trials = nr_all, iterations = N_ITER, burn_in = BURN,
      retained_draws_per_trial = N_ITER - BURN,
      gamma_acceptance = fit$gamma_acceptance,
      sigma_acceptance = fit$sigma_acceptance
    )
  }
  cat("posterior complete for scenario ", s, "\n", sep = "")
}
write.csv(do.call(rbind, diagnostic_rows),
          file.path(script_dir, "oc_diagnostics.csv"), row.names = FALSE)

## ---- Threshold calibration on calibration trials only -------------
empirical_cutoff <- function(values, target = 0.05) {
  as.numeric(quantile(values, 1 - target, type = 1, names = FALSE))
}
Q <- list(weak = c(), strong = c())
for (lab in c("No", "Moderate", "Strong")) {
  Q$weak[lab] <- empirical_cutoff(apply(prob_cal[[paste0("1.", lab)]], 1, max))
  candidates <- c()
  for (s in names(SCEN)) {
    nulls <- nullset[[s]]
    if (length(nulls)) {
      candidates <- c(
        candidates,
        empirical_cutoff(apply(prob_cal[[paste0(s, ".", lab)]][, nulls, drop = FALSE], 1, max))
      )
    }
  }
  Q$strong[lab] <- max(candidates)
}

threshold_rows <- do.call(rbind, lapply(c("weak", "strong"), function(control) {
  do.call(rbind, lapply(c("No", "Moderate", "Strong"), function(lab) {
    data.frame(control = control, borrow = lab, threshold = Q[[control]][lab],
               calibration_trials_per_scenario = R_CAL,
               target_fwer = 0.05)
  }))
}))
write.csv(threshold_rows, file.path(script_dir, "oc_thresholds.csv"), row.names = FALSE)
cat("\nCalibrated thresholds q (independent calibration set):\n")
print(round(rbind(weak = Q$weak, strong = Q$strong), 4))

## ---- Evaluation metrics with Wilson 95% Monte Carlo intervals -----
wilson <- function(indicator, z = qnorm(0.975)) {
  n_trials <- length(indicator)
  estimate <- mean(indicator)
  denominator <- 1 + z^2 / n_trials
  center <- (estimate + z^2 / (2 * n_trials)) / denominator
  half_width <- z * sqrt(estimate * (1 - estimate) / n_trials +
                         z^2 / (4 * n_trials^2)) / denominator
  c(estimate = estimate, low = max(0, center - half_width),
    high = min(1, center + half_width)) * 100
}

metrics <- function(reject, scenario) {
  active <- promising[[scenario]]
  nulls <- nullset[[scenario]]
  basket <- lapply(seq_len(NB), function(j) wilson(reject[, j]))
  fwer <- if (length(nulls)) wilson(apply(reject[, nulls, drop = FALSE], 1, any)) else rep(NA_real_, 3)
  fwpd <- if (length(active)) wilson(apply(reject[, active, drop = FALSE], 1, any)) else rep(NA_real_, 3)
  fwpc <- if (length(active)) wilson(apply(reject[, active, drop = FALSE], 1, all)) else rep(NA_real_, 3)
  list(basket = basket, FWER = fwer, FWPD = fwpd, FWPC = fwpc)
}

RES <- list()
for (control in c("weak", "strong")) {
  cat(sprintf("\n===== %s FWER TARGET: INDEPENDENT EVALUATION =====\n", toupper(control)))
  cat(sprintf("%-5s %-10s %-28s %8s %8s %8s\n",
              "Scen", "Borrow", "%reject (b1,b2,b3,b4)", "FWER", "FWP-D", "FWP-C"))
  for (s in names(SCEN)) for (lab in c("No", "Moderate", "Strong")) {
    reject <- prob_eval[[paste0(s, ".", lab)]] > Q[[control]][lab]
    result <- metrics(reject, s)
    RES[[paste(control, s, lab, sep = ".")]] <- result
    basket_est <- vapply(result$basket, `[[`, numeric(1), "estimate")
    fmt <- function(x) if (is.na(x)) "   -   " else sprintf("%6.1f", x)
    cat(sprintf("S%-4s %-10s %-28s %8s %8s %8s\n", s, lab,
                paste(sprintf("%4.1f", basket_est), collapse = ","),
                fmt(result$FWER["estimate"]), fmt(result$FWPD["estimate"]),
                fmt(result$FWPC["estimate"])))
  }
}

result_rows <- do.call(rbind, lapply(names(RES), function(key) {
  parts <- strsplit(key, ".", fixed = TRUE)[[1]]
  result <- RES[[key]]
  row <- data.frame(
    control = parts[1], scenario = as.integer(parts[2]), borrow = parts[3],
    threshold = Q[[parts[1]]][parts[3]], calibration_trials = R_CAL,
    evaluation_trials = R_EVAL
  )
  for (j in seq_len(NB)) {
    row[[paste0("b", j)]] <- result$basket[[j]]["estimate"]
    row[[paste0("b", j, "_low")]] <- result$basket[[j]]["low"]
    row[[paste0("b", j, "_high")]] <- result$basket[[j]]["high"]
  }
  for (metric in c("FWER", "FWPD", "FWPC")) {
    row[[metric]] <- result[[metric]]["estimate"]
    row[[paste0(metric, "_low")]] <- result[[metric]]["low"]
    row[[paste0(metric, "_high")]] <- result[[metric]]["high"]
  }
  row
}))
numeric_columns <- vapply(result_rows, is.numeric, logical(1))
result_rows[numeric_columns] <- lapply(result_rows[numeric_columns], round, 4)
write.csv(result_rows, file.path(script_dir, "oc_results.csv"), row.names = FALSE)

## ---- Figures: published values and independent reimplementation ---
get_metric <- function(control, scenario, metric) {
  vapply(c("No", "Moderate", "Strong"), function(lab) {
    RES[[paste(control, scenario, lab, sep = ".")]][[metric]]["estimate"]
  }, numeric(1))
}
get_interval <- function(control, scenario, metric) {
  rbind(
    low = vapply(c("No", "Moderate", "Strong"), function(lab) {
      RES[[paste(control, scenario, lab, sep = ".")]][[metric]]["low"]
    }, numeric(1)),
    high = vapply(c("No", "Moderate", "Strong"), function(lab) {
      RES[[paste(control, scenario, lab, sep = ".")]][[metric]]["high"]
    }, numeric(1))
  )
}

teal <- "#3E7C7B"; rust <- "#C8613A"; gold <- "#B08A3E"; grey <- "#9A948C"
plot_comparison <- function(path, homogeneous, mixed_fwer, one_active,
                            source_note, intervals = NULL) {
  png(path, width = 1950, height = 680, res = 150)
  op <- par(mfrow = c(1, 3), mar = c(5, 5, 4.5, 2), family = "sans")
  bar_panel <- function(values, colour, main, subtitle, ymax, ylab,
                        reference = NA, interval = NULL) {
    bars <- barplot(values, col = colour, border = NA, ylim = c(0, ymax),
                    names.arg = c("No", "Moderate", "Strong"),
                    ylab = ylab, xlab = "Degree of borrowing")
    if (!is.null(interval)) {
      arrows(bars, interval["low", ], bars, interval["high", ],
             angle = 90, code = 3, length = 0.04, col = "#20201E", lwd = 1.1)
    }
    labels <- ifelse(abs(values - round(values)) < 0.05,
                     sprintf("%.0f", values), sprintf("%.1f", values))
    label_height <- if (is.null(interval)) values else interval["high", ]
    text(bars, label_height, labels, pos = 3, xpd = NA)
    if (!is.na(reference)) {
      abline(h = reference, lty = 2, col = grey)
      text(max(bars), reference, "5% target", pos = 3, col = grey, cex = 0.8)
    }
    title(main = main, adj = 0, line = 2.4, cex.main = 1.05)
    mtext(subtitle, adj = 0, line = 0.7, col = grey, cex = 0.72)
  }
  bar_panel(homogeneous, teal, "Homogeneous: borrowing wins",
            paste("Scenario 2 - all baskets active -", source_note),
            100, "Disjunctive power (%)",
            interval = if (is.null(intervals)) NULL else intervals$homogeneous)
  bar_panel(mixed_fwer, rust, "Mixed: type I error inflates",
            paste("Scenario 3 - one inactive basket -", source_note),
            40, "Family-wise type I error (%)", reference = 5,
            interval = if (is.null(intervals)) NULL else intervals$mixed_fwer)
  bar_panel(one_active, gold, "Clamp the error, lose the power",
            paste("Scenario 6 - one active basket -", source_note),
            40, "Disjunctive power (%)",
            interval = if (is.null(intervals)) NULL else intervals$one_active)
  par(op)
  dev.off()
}

plot_comparison(
  file.path(script_dir, "figure_oc_borrowing_reimplementation.png"),
  get_metric("weak", "2", "FWPD"),
  get_metric("weak", "3", "FWER"),
  get_metric("strong", "6", "FWPD"),
  "independent evaluation (n=5,000)",
  intervals = list(
    homogeneous = get_interval("weak", "2", "FWPD"),
    mixed_fwer = get_interval("weak", "3", "FWER"),
    one_active = get_interval("strong", "6", "FWPD")
  )
)

plot_comparison(
  file.path(script_dir, "figure_oc_borrowing.png"),
  c(69.6, 87.0, 90.3),
  c(1.1, 8.7, 33.9),
  c(27.0, 16.9, 4.0),
  "published values"
)

cat("\nWrote figures, oc_results.csv, oc_thresholds.csv, oc_diagnostics.csv, and oc_simulation_inputs.csv\n")

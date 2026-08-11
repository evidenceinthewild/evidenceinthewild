# =====================================================================
# imatinib_shrinkage.R
# Reimplements the shrinkage illustration (Zhou & Ji 2024, Figure 2):
# stratified Clopper-Pearson intervals vs. a Bayesian hierarchical
# logit-increment model, on the 10-basket imatinib sarcoma trial.
#
# Four independent chains are used. The script reports acceptance rates and
# classical Gelman-Rubin R-hat values, and writes machine-readable results.
#
# Base R only. Run: Rscript imatinib_shrinkage.R
# =====================================================================

expit <- plogis
logit <- qlogis

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) {
  script_path <- gsub("~\\+~", " ", sub("^--file=", "", file_arg[1]))
  dirname(normalizePath(script_path))
} else {
  getwd()
}

## ---- Table 1 data (imatinib) --------------------------------------
subtype <- c("Angiosarcoma", "Ewing", "Fibrosarcoma", "Leiomyosarcoma",
             "Liposarcoma", "MFH", "Osteosarcoma", "MPNST",
             "Rhabdomyosarcoma", "Synovial")
y <- c(2, 0, 1, 6, 7, 3, 5, 1, 0, 3)
n <- c(15, 13, 12, 28, 29, 29, 26, 5, 2, 20)
J <- length(y)
pi0 <- 0.30
l0 <- logit(pi0)
overall <- sum(y) / sum(n)

## ---- Stratified analysis: y/n with 95% Clopper-Pearson exact CI ---
a <- 0.05
strat_pt <- y / n
cp_lo <- qbeta(a / 2, y, n - y + 1); cp_lo[y == 0] <- 0
cp_hi <- qbeta(1 - a / 2, y + 1, n - y); cp_hi[y == n] <- 1

## ---- Bayesian hierarchical model (Figure 2 specification) ---------
## gamma_j = logit(pi_j) - logit(pi0); gamma_j ~ N(mu, sigma^2)
## mu ~ N(0, 100^2); sigma ~ Half-N(3)
loglik <- function(g) {
  p <- pmin(pmax(expit(l0 + g), 1e-12), 1 - 1e-12)
  y * log(p) + (n - y) * log(1 - p)
}

n_chains <- 4L
n_iter <- 40000L
burn <- 10000L
A_sig <- 3.0
step_g <- 0.6
step_ls <- 0.4
chain_seeds <- 20260707L + seq_len(n_chains) - 1L

run_chain <- function(seed) {
  set.seed(seed)
  gamma <- rep(0, J); mu <- 0; sigma <- 1
  ll <- loglik(gamma)
  keep_pi <- matrix(NA_real_, n_iter - burn, J)
  keep_sigma <- numeric(n_iter - burn)
  gamma_accept <- 0
  sigma_accept <- 0
  k <- 0L

  for (it in seq_len(n_iter)) {
    prop <- gamma + rnorm(J, 0, step_g)
    llp <- loglik(prop)
    lp_old <- -0.5 * ((gamma - mu) / sigma)^2
    lp_new <- -0.5 * ((prop - mu) / sigma)^2
    accept <- log(runif(J)) < (llp + lp_new) - (ll + lp_old)
    gamma_accept <- gamma_accept + sum(accept)
    gamma[accept] <- prop[accept]
    ll[accept] <- llp[accept]

    prec <- 1 / 1e4 + J / sigma^2
    mu <- rnorm(1, (sum(gamma) / sigma^2) / prec, sqrt(1 / prec))

    ls <- log(sigma); lsp <- ls + rnorm(1, 0, step_ls); sp <- exp(lsp)
    log_target <- function(s, log_s) {
      -J * log(s) - 0.5 * sum((gamma - mu)^2) / s^2 -
        s^2 / (2 * A_sig^2) + log_s
    }
    accept_sigma <- log(runif(1)) <
      log_target(sp, lsp) - log_target(sigma, ls)
    if (accept_sigma) sigma <- sp
    sigma_accept <- sigma_accept + accept_sigma

    if (it > burn) {
      k <- k + 1L
      keep_pi[k, ] <- expit(l0 + gamma)
      keep_sigma[k] <- sigma
    }
  }

  list(
    pi = keep_pi,
    sigma = keep_sigma,
    gamma_accept = gamma_accept / (n_iter * J),
    sigma_accept = sigma_accept / n_iter
  )
}

chains <- lapply(chain_seeds, run_chain)
draws <- do.call(rbind, lapply(chains, `[[`, "pi"))

gelman_rhat <- function(chain_values) {
  chain_values <- lapply(chain_values, as.numeric)
  n_draws <- min(vapply(chain_values, length, integer(1)))
  chain_values <- lapply(chain_values, function(x) x[seq_len(n_draws)])
  within <- mean(vapply(chain_values, var, numeric(1)))
  between <- n_draws * var(vapply(chain_values, mean, numeric(1)))
  variance_plus <- ((n_draws - 1) / n_draws) * within + between / n_draws
  sqrt(variance_plus / within)
}

rhat_pi <- vapply(seq_len(J), function(j) {
  gelman_rhat(lapply(chains, function(x) x$pi[, j]))
}, numeric(1))
rhat_sigma <- gelman_rhat(lapply(chains, `[[`, "sigma"))

post_mean <- colMeans(draws)
post_lo <- apply(draws, 2, quantile, 0.025)
post_hi <- apply(draws, 2, quantile, 0.975)
w_strat <- mean(cp_hi - cp_lo)
w_borrow <- mean(post_hi - post_lo)
width_reduction <- 1 - w_borrow / w_strat

## ---- Reports -------------------------------------------------------
display_tab <- data.frame(
  subtype, y, n,
  strat = round(strat_pt * 100, 1),
  strat_CI = sprintf("[%.0f, %.0f]", cp_lo * 100, cp_hi * 100),
  borrow = round(post_mean * 100, 1),
  borrow_CI = sprintf("[%.0f, %.0f]", post_lo * 100, post_hi * 100)
)
cat(sprintf("Overall observed response rate: %.1f%%\n", overall * 100))
print(display_tab, row.names = FALSE)
cat(sprintf("\nMean 95%% interval width  stratified: %.1f pts   borrowing: %.1f pts\n",
            w_strat * 100, w_borrow * 100))
cat(sprintf("Intervals are %.0f%% narrower on average under borrowing.\n",
            width_reduction * 100))
cat(sprintf("Four-chain diagnostics: max R-hat %.4f; gamma acceptance %.1f%%; sigma acceptance %.1f%%\n",
            max(c(rhat_pi, rhat_sigma)),
            mean(vapply(chains, `[[`, numeric(1), "gamma_accept")) * 100,
            mean(vapply(chains, `[[`, numeric(1), "sigma_accept")) * 100))

results <- data.frame(
  subtype, y, n,
  stratified_rate = strat_pt,
  stratified_ci_low = cp_lo,
  stratified_ci_high = cp_hi,
  posterior_mean = post_mean,
  posterior_ci_low = post_lo,
  posterior_ci_high = post_hi,
  rhat = rhat_pi
)
write.csv(results, file.path(script_dir, "imatinib_results.csv"), row.names = FALSE)

diagnostics <- data.frame(
  implementation = "R",
  chains = n_chains,
  iterations_per_chain = n_iter,
  burn_in_per_chain = burn,
  retained_draws = nrow(draws),
  mean_gamma_acceptance = mean(vapply(chains, `[[`, numeric(1), "gamma_accept")),
  mean_sigma_acceptance = mean(vapply(chains, `[[`, numeric(1), "sigma_accept")),
  max_response_rate_rhat = max(rhat_pi),
  sigma_rhat = rhat_sigma,
  interval_width_reduction = width_reduction,
  r_version = R.version.string
)
write.csv(diagnostics, file.path(script_dir, "imatinib_diagnostics.csv"), row.names = FALSE)

## ---- Forest plot ---------------------------------------------------
figure_path <- file.path(script_dir, "figure_imatinib_shrinkage.png")
png(figure_path, width = 1300, height = 850, res = 150)
op <- par(mar = c(4.5, 10, 3, 2))
yy <- rev(seq_len(J)); off <- 0.16
plot(NA, xlim = c(0, 90), ylim = c(0.5, J + 0.5), yaxt = "n",
     xlab = "Response rate (%)", ylab = "",
     main = "Borrowing shrinks basket estimates toward 15.6%")
abline(v = overall * 100, lty = 2, col = "grey60")
rust <- "#C8613A"; teal <- "#3E7C7B"
for (i in seq_len(J)) {
  segments(cp_lo[i] * 100, yy[i] + off, cp_hi[i] * 100, yy[i] + off,
           col = rust, lwd = 2)
  points(strat_pt[i] * 100, yy[i] + off, pch = 19, col = rust, cex = 0.9)
  segments(post_lo[i] * 100, yy[i] - off, post_hi[i] * 100, yy[i] - off,
           col = teal, lwd = 2)
  points(post_mean[i] * 100, yy[i] - off, pch = 18, col = teal, cex = 1.1)
}
axis(2, at = yy, labels = sprintf("%s (%d)", subtype, n), las = 1, cex.axis = 0.8)
legend("bottomright", bty = "n", pch = c(19, 18), col = c(rust, teal),
       legend = c("Stratified (Clopper-Pearson 95% CI)",
                  "Bayesian hierarchical (95% CrI)"))
par(op); dev.off()
cat("Wrote figure_imatinib_shrinkage.png, imatinib_results.csv, and imatinib_diagnostics.csv\n")

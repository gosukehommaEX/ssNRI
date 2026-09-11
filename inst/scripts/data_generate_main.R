# =============================================================================
# data_generate_main.R
#
# Computation stage for the main-text figures and table of:
#   "A Cautionary Note on Sample Size Inflation for Trials with
#    Non-Responder Imputation"
#
# All numerical computation happens here and the results are written to
# "data/*.rds". The companion script table_and_figure_manuscript.R reads those
# files and renders the figures and the LaTeX table without recomputing
# anything, so labels, scales and captions can be adjusted freely.
#
# Everything is deterministic: exact power and exact type I error rate come
# from complete enumeration, not simulation, so no seed is involved and
# separating computation from rendering cannot change a number.
#
# Usage
#   Set the working directory to this script's location (inst/scripts), then
#   run the whole file. The package must be loaded first:
#     library(ssNRI)               # installed package, or
#     devtools::load_all()         # from the package root during development
#
# Outputs (relative to the working directory)
#   data/fig1_data.rds
#   data/fig2_data.rds
#   data/fig3_data.rds
#   data/table2_data.rds
#   data/application_data.rds
#   data/data_generate_main.log
#   data/sessionInfo_data_generate_main.txt
#
# Changes from the submitted version
#   Figure 1 now starts at omega = 0.05 rather than omega = 0. At omega = 0 the
#   completer count is degenerate and the exact power and type I error rate
#   behave discontinuously, which Reviewer 1 asked to remove (R1-m5).
#   Table 2 is reduced to the design assumptions of the original trial, with
#   the wider sensitivity grid moved to the supplementary material (AE-M4), and
#   gains the NRI-scale probabilities and effect (R2-M2a, R2-m4).
# =============================================================================

library(ssNRI)

# --------------------------------------------------------------------------- #
# Output location (relative paths only)
# --------------------------------------------------------------------------- #
data_dir <- "data"
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

log_file <- file.path(data_dir, "data_generate_main.log")
cat(sprintf("data_generate_main.R log -- %s\n", format(Sys.time())),
    file = log_file, append = FALSE)

log_msg <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n", sep = "", file = log_file, append = TRUE)
  message(msg)
}

timed <- function(label, expr) {
  t <- system.time(val <- force(expr))["elapsed"]
  log_msg("%-34s %8.1f sec", label, t)
  val
}

# --------------------------------------------------------------------------- #
# Fixed design parameters shared by every main-text display
# --------------------------------------------------------------------------- #
ALPHA <- 0.025
TARGET_POWER <- 0.80
R_ALLOC <- 1

# --------------------------------------------------------------------------- #
# Figure 1: validity of the simple inflation method under a CC analysis with
# within-group independence between response and dropout
# --------------------------------------------------------------------------- #
FIG1_DELTA <- c(0.15, 0.25)
FIG1_PI1 <- c(0.35, 0.45, 0.55)
FIG1_OMEGA <- seq(0.05, 0.30, by = 0.05)

fig1_grid <- expand.grid(
  delta_val = FIG1_DELTA,
  pi1_idx = seq_along(FIG1_PI1),
  omega = FIG1_OMEGA,
  stringsAsFactors = FALSE
)

fig1_df <- timed("Figure 1 (CC validity)", {
  rows <- lapply(seq_len(nrow(fig1_grid)), function(k) {
    pi1_val <- FIG1_PI1[fig1_grid$pi1_idx[k]]
    pi0_val <- pi1_val - fig1_grid$delta_val[k]
    omega <- fig1_grid$omega[k]

    ss <- sample_size_cc(
      pi1 = pi1_val, pi0 = pi0_val,
      omega1 = omega, omega0 = omega,
      r = R_ALLOC, alpha = ALPHA, target_power = TARGET_POWER
    )
    perf <- power_cc(
      n1 = ss$n1_adjust, n0 = ss$n0_adjust,
      pi1 = pi1_val, pi0 = pi0_val,
      omega1 = omega, omega0 = omega,
      rho1 = 0, rho0 = 0, alpha = ALPHA
    )
    data.frame(
      case = paste0("Case_", LETTERS[fig1_grid$pi1_idx[k]], "_p1_", pi1_val),
      case_label = LETTERS[fig1_grid$pi1_idx[k]],
      pi1_value = pi1_val,
      delta_val = fig1_grid$delta_val[k],
      omega = omega,
      n_total = ss$n1_adjust + ss$n0_adjust,
      power = perf$power,
      type1_error = perf$type1_error,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
})

saveRDS(
  list(
    df = fig1_df,
    delta = FIG1_DELTA,
    pi1 = FIG1_PI1,
    omega_values = FIG1_OMEGA,
    r = R_ALLOC,
    alpha = ALPHA,
    target_power = TARGET_POWER
  ),
  file.path(data_dir, "fig1_data.rds")
)

# --------------------------------------------------------------------------- #
# Figure 2: proposed NRI formula against the simple inflation method, as a
# function of the response-dropout correlation
# --------------------------------------------------------------------------- #
FIG2_DELTA <- 0.15
FIG2_PI1 <- c(0.35, 0.45, 0.55)
FIG2_OMEGA <- c(0.05, 0.10, 0.20)
FIG2_RHO_STEP <- 0.01

# Correlation grid: the range feasible in both groups, per (pi1, omega) cell
fig2_tasks <- list()
k <- 1L
for (i in seq_along(FIG2_PI1)) {
  pi1_val <- FIG2_PI1[i]
  pi0_val <- pi1_val - FIG2_DELTA
  for (omega in FIG2_OMEGA) {
    b1 <- rho_bounds(pi = pi1_val, omega = omega)
    b0 <- rho_bounds(pi = pi0_val, omega = omega)
    rho_min <- max(b1$rho_lower, b0$rho_lower)
    rho_max <- min(b1$rho_upper, b0$rho_upper)
    if (rho_min > rho_max) {
      log_msg("No feasible correlation range for Case %s, omega = %.2f",
              LETTERS[i], omega)
      next
    }
    for (rho in seq(rho_min, rho_max, by = FIG2_RHO_STEP)) {
      fig2_tasks[[k]] <- list(pi1_idx = i, pi1_val = pi1_val, pi0_val = pi0_val,
                              omega = omega, rho = rho)
      k <- k + 1L
    }
  }
}
log_msg("Figure 2 correlation grid: %d points", length(fig2_tasks))

fig2_df <- timed("Figure 2 (NRI vs simple inflation)", {
  rows <- lapply(fig2_tasks, function(task) {
    ss_nri <- sample_size_nri(
      pi1 = task$pi1_val, pi0 = task$pi0_val,
      omega1 = task$omega, omega0 = task$omega,
      rho1 = task$rho, rho0 = task$rho,
      r = R_ALLOC, alpha = ALPHA, target_power = TARGET_POWER
    )
    ss_cc <- sample_size_cc(
      pi1 = task$pi1_val, pi0 = task$pi0_val,
      omega1 = task$omega, omega0 = task$omega,
      r = R_ALLOC, alpha = ALPHA, target_power = TARGET_POWER
    )
    N_proposed <- ss_nri$n_total
    N_simple <- ss_cc$n_total_adjust

    N1_proposed <- ceiling(N_proposed * R_ALLOC / (1 + R_ALLOC))
    N0_proposed <- N_proposed - N1_proposed
    N1_simple <- ceiling(N_simple * R_ALLOC / (1 + R_ALLOC))
    N0_simple <- N_simple - N1_simple

    # Both sample sizes are evaluated under the NRI analysis, which is the
    # analysis the trial will actually run
    ep <- power_nri(
      n1 = N1_proposed, n0 = N0_proposed,
      pi1 = task$pi1_val, pi0 = task$pi0_val,
      omega1 = task$omega, omega0 = task$omega,
      rho1 = task$rho, rho0 = task$rho, alpha = ALPHA
    )
    es <- power_nri(
      n1 = N1_simple, n0 = N0_simple,
      pi1 = task$pi1_val, pi0 = task$pi0_val,
      omega1 = task$omega, omega0 = task$omega,
      rho1 = task$rho, rho0 = task$rho, alpha = ALPHA
    )

    data.frame(
      rho = task$rho,
      case = paste0("Case_", LETTERS[task$pi1_idx], "_p1_", task$pi1_val),
      pi1_value = task$pi1_val,
      omega = task$omega,
      N_proposed = N_proposed,
      N_Simple = N_simple,
      power_proposed = ep$power,
      power_Simple = es$power,
      type1_proposed = ep$type1_error,
      type1_Simple = es$type1_error,
      pi1_NRI = ep$pi1_NRI,
      pi0_NRI = ep$pi0_NRI,
      effect_NRI = ep$effect_NRI,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
})

saveRDS(
  list(
    df = fig2_df,
    delta = FIG2_DELTA,
    pi1 = FIG2_PI1,
    omega_values = FIG2_OMEGA,
    rho_step = FIG2_RHO_STEP,
    r = R_ALLOC,
    alpha = ALPHA,
    target_power = TARGET_POWER
  ),
  file.path(data_dir, "fig2_data.rds")
)

# --------------------------------------------------------------------------- #
# Figure 3: relative efficiency over the (pi1, omega) plane
# --------------------------------------------------------------------------- #
FIG3_DELTA <- c(0.05, 0.10, 0.15, 0.20)
FIG3_PI1 <- seq(0.30, 0.90, by = 0.01)
FIG3_OMEGA <- seq(0.05, 0.50, by = 0.01)
FIG3_CONTOUR_STEP <- 0.1

fig3_grid <- expand.grid(
  delta = FIG3_DELTA,
  pi1 = FIG3_PI1,
  omega = FIG3_OMEGA,
  stringsAsFactors = FALSE
)

fig3_df <- timed("Figure 3 (relative efficiency)", {
  n_cc <- numeric(nrow(fig3_grid))
  n_nri <- numeric(nrow(fig3_grid))
  for (k in seq_len(nrow(fig3_grid))) {
    pi1_val <- fig3_grid$pi1[k]
    pi0_val <- pi1_val - fig3_grid$delta[k]
    omega <- fig3_grid$omega[k]
    n_cc[k] <- sample_size_cc(
      pi1 = pi1_val, pi0 = pi0_val, omega1 = omega, omega0 = omega,
      r = R_ALLOC, alpha = ALPHA, target_power = TARGET_POWER
    )$n_total_adjust
    n_nri[k] <- sample_size_nri(
      pi1 = pi1_val, pi0 = pi0_val, omega1 = omega, omega0 = omega,
      rho1 = 0, rho0 = 0,
      r = R_ALLOC, alpha = ALPHA, target_power = TARGET_POWER
    )$n_total
  }
  data.frame(
    delta = fig3_grid$delta,
    pi1 = fig3_grid$pi1,
    omega = fig3_grid$omega,
    n_cc = n_cc,
    n_nri = n_nri,
    re = n_cc / n_nri,
    delta_label = paste0("delta[Full] == ", fig3_grid$delta),
    stringsAsFactors = FALSE
  )
})

saveRDS(
  list(
    df = fig3_df,
    delta = FIG3_DELTA,
    pi1_values = FIG3_PI1,
    omega_values = FIG3_OMEGA,
    contour_step = FIG3_CONTOUR_STEP,
    r = R_ALLOC,
    alpha = ALPHA,
    target_power = TARGET_POWER
  ),
  file.path(data_dir, "fig3_data.rds")
)

# --------------------------------------------------------------------------- #
# Table 2: application to the cytisine smoking cessation trial
#
# Reduced to the design assumptions the original trial actually used. The wider
# grid of pi1 and omega values is now supplementary material. The NRI-scale
# response probabilities and effect are reported alongside the sample sizes, so
# that the reader can see why the required sample size differs. Each correlation
# is also reported as the probability of responding among dropouts, gamma, which
# is what the design assumption amounts to in the two groups.
# --------------------------------------------------------------------------- #
APP_PI1 <- 0.47
APP_PI0 <- 0.41
APP_OMEGA <- 0.10
APP_N_REPORTED <- 2388   # total enrolled under the simple inflation method

table2_df <- timed("Table 2 (cytisine application)", {
  b1 <- rho_bounds(pi = APP_PI1, omega = APP_OMEGA)
  b0 <- rho_bounds(pi = APP_PI0, omega = APP_OMEGA)
  rho_l <- max(b1$rho_lower, b0$rho_lower)
  rho_u <- min(b1$rho_upper, b0$rho_upper)

  ss_cc <- sample_size_cc(
    pi1 = APP_PI1, pi0 = APP_PI0, omega1 = APP_OMEGA, omega0 = APP_OMEGA,
    r = R_ALLOC, alpha = ALPHA, target_power = TARGET_POWER
  )
  N_cc <- ss_cc$n_total_adjust

  rows <- lapply(seq_len(3), function(k) {
    rho_k <- c(rho_l, 0, rho_u)[k]
    ss_nri <- sample_size_nri(
      pi1 = APP_PI1, pi0 = APP_PI0, omega1 = APP_OMEGA, omega0 = APP_OMEGA,
      rho1 = rho_k, rho0 = rho_k,
      r = R_ALLOC, alpha = ALPHA, target_power = TARGET_POWER
    )
    data.frame(
      pi1 = APP_PI1,
      pi0 = APP_PI0,
      omega = APP_OMEGA,
      rho = rho_k,
      rho_label = c("L", "0", "U")[k],
      gamma1 = rho_to_gamma(pi = APP_PI1, omega = APP_OMEGA, rho = rho_k),
      gamma0 = rho_to_gamma(pi = APP_PI0, omega = APP_OMEGA, rho = rho_k),
      pi1_NRI = ss_nri$pi1_NRI,
      pi0_NRI = ss_nri$pi0_NRI,
      effect_full = ss_nri$effect_full,
      effect_NRI = ss_nri$effect_NRI,
      N_cc = N_cc,
      N_nri = ss_nri$n_total,
      RE = N_cc / ss_nri$n_total,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
})

saveRDS(
  list(
    df = table2_df,
    pi1 = APP_PI1,
    pi0 = APP_PI0,
    omega = APP_OMEGA,
    n_reported = APP_N_REPORTED,
    r = R_ALLOC,
    alpha = ALPHA,
    target_power = TARGET_POWER
  ),
  file.path(data_dir, "table2_data.rds")
)

# --------------------------------------------------------------------------- #
# Consistency check: the simple inflation method must reproduce the sample size
# the original trial reported, which is what identifies the formula it used
# --------------------------------------------------------------------------- #
if (table2_df$N_cc[1] != APP_N_REPORTED) {
  log_msg("WARNING: N_CC = %d does not match the reported %d",
          table2_df$N_cc[1], APP_N_REPORTED)
} else {
  log_msg("N_CC reproduces the sample size reported by Dogar et al. (%d)",
          APP_N_REPORTED)
}

# --------------------------------------------------------------------------- #
# What the trial actually enrolled and observed (AE-M5)
#
# The associate editor asked for the trial's results, not only its design. The
# quantities computed here are design-stage: the exact power the trial had
# against the effect it was designed to detect, under an NRI analysis, at the
# number of patients it planned and at the number it actually enrolled. They
# are not post hoc power, which would be a function of the observed effect and
# would carry no information beyond the p value.
#
# The observed quantities are recorded verbatim from the published report so
# that number_check.R can compare them with what the article states.
# --------------------------------------------------------------------------- #
APP_N1_ENROLLED <- 1239
APP_N0_ENROLLED <- 1233
APP_N1_PLANNED  <- 1194
APP_N0_PLANNED  <- 1194

# Reported results at six months, biochemically verified continuous abstinence
APP_OBS <- list(
  x1 = 401, x0 = 366,
  rate1 = 0.324, rate0 = 0.297,
  risk_difference = 0.0268,
  ci_lower = -0.0096, ci_upper = 0.0633,
  relative_risk = 1.09, rr_lower = 0.97, rr_upper = 1.23,
  p_value = 0.114,
  # Patients with follow-up data available; the rest were counted as failures
  n1_followed = 1142, n0_followed = 1130
)
APP_OBS$omega1_observed <- 1 - APP_OBS$n1_followed / APP_N1_ENROLLED
APP_OBS$omega0_observed <- 1 - APP_OBS$n0_followed / APP_N0_ENROLLED

application_power <- timed("Application (realized design power)", {
  sizes <- list(
    planned  = c(APP_N1_PLANNED,  APP_N0_PLANNED),
    enrolled = c(APP_N1_ENROLLED, APP_N0_ENROLLED),
    required = c(ceiling(table2_df$N_nri[table2_df$rho_label == "0"] / 2),
                 ceiling(table2_df$N_nri[table2_df$rho_label == "0"] / 2))
  )
  rows <- lapply(names(sizes), function(nm) {
    n <- sizes[[nm]]
    res <- power_nri(
      n1 = n[1], n0 = n[2], pi1 = APP_PI1, pi0 = APP_PI0,
      omega1 = APP_OMEGA, omega0 = APP_OMEGA,
      rho1 = 0, rho0 = 0, alpha = ALPHA
    )
    data.frame(
      scenario = nm, n1 = n[1], n0 = n[2], n_total = n[1] + n[2],
      pi1_NRI = res$pi1_NRI, pi0_NRI = res$pi0_NRI,
      effect_NRI = res$effect_NRI, power = res$power,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
})

saveRDS(
  list(power = application_power, observed = APP_OBS,
       pi1 = APP_PI1, pi0 = APP_PI0, omega = APP_OMEGA,
       n_planned = APP_N1_PLANNED + APP_N0_PLANNED,
       n_enrolled = APP_N1_ENROLLED + APP_N0_ENROLLED,
       alpha = ALPHA, target_power = TARGET_POWER),
  file.path(data_dir, "application_data.rds")
)

log_msg("Exact NRI power at the planned size  (%d) = %.4f",
        application_power$n_total[1], application_power$power[1])
log_msg("Exact NRI power at the enrolled size (%d) = %.4f",
        application_power$n_total[2], application_power$power[2])
log_msg("Exact NRI power at the required size (%d) = %.4f",
        application_power$n_total[3], application_power$power[3])

# --------------------------------------------------------------------------- #
# Record the execution environment for reproducibility
# --------------------------------------------------------------------------- #
writeLines(
  capture.output(sessionInfo()),
  file.path(data_dir, "sessionInfo_data_generate_main.txt")
)

log_msg("All intermediate data written to '%s/'.", data_dir)

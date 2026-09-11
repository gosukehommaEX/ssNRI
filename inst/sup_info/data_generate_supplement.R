# =============================================================================
# data_generate_supplement.R
#
# Computation stage for the supplementary material of:
#   "A Cautionary Note on Sample Size Inflation for Trials with
#    Non-Responder Imputation"
#
# Every display answers a specific referee comment:
#
#   Table S1   sensitivity grid over pi1 and omega, moved out of the main text
#              at the associate editor's request                      (AE-M4)
#   Figure S1  unequal dropout probabilities under the alternative  (R2-M2b)
#   Figure S2  unequal dropout probabilities under the full-data null,
#              where the NRI estimand is not null even though the latent
#              response probabilities are equal          (R2-M2b, AE on MAR)
#   Figure S3  unequal response-dropout correlations across groups  (R1-M3a)
#   Figure S4  the region where the NRI effect exceeds the full-data effect,
#              as the two dropout probabilities vary                 (R2-M1)
#   Figure S5  the same ordering question as the two correlations vary
#                                                          (R2-M1, R1-M3a)
#   Figure S6  conservative dropout assumptions as an alternative way to
#              protect against underpowering            (R1-M4, AE-m15)
#   Table S2   allocation ratios other than one                     (R1-M3b)
#   Table S3   type I error rate at several null nuisance parameters (R2-m1)
#   Table S4   unequal dropout probabilities with per-group inflation, a
#              counterexample to the necessity of the three conditions
#                                                            (R2-m3, AE-M2)
#   Figure S7  partial identification of the latent response probability from
#              a historical rate, carried through to the sample size
#                                                            (R1-M2, R2-M3)
#
# Usage
#   Set the working directory to this script's location (inst/sup_info), then
#   run the whole file, with the package loaded:
#     library(ssNRI)               # or devtools::load_all() from the root
#
# Outputs (relative to the working directory)
#   data/tableS1_data.rds ... data/tableS4_data.rds
#   data/figS1_data.rds   ... data/figS7_data.rds
#   data/data_generate_supplement.log
#   data/sessionInfo_data_generate_supplement.txt
# =============================================================================

library(ssNRI)

data_dir <- "data"
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

log_file <- file.path(data_dir, "data_generate_supplement.log")
cat(sprintf("data_generate_supplement.R log -- %s\n", format(Sys.time())),
    file = log_file, append = FALSE)

log_msg <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n", sep = "", file = log_file, append = TRUE)
  message(msg)
}

timed <- function(label, expr) {
  t <- system.time(val <- force(expr))["elapsed"]
  log_msg("%-40s %8.1f sec", label, t)
  val
}

# --------------------------------------------------------------------------- #
# Shared parameters. The reference design matches Case B of Figure 2, so that
# the supplementary displays connect directly to the main text.
# --------------------------------------------------------------------------- #
ALPHA <- 0.025
TARGET_POWER <- 0.80
REF_PI1 <- 0.45
REF_PI0 <- 0.30
REF_OMEGA <- 0.10

# Feasible correlation range in both groups at a given pair of probabilities
common_rho <- function(pi1, pi0, omega1, omega0 = omega1) {
  b1 <- rho_bounds(pi = pi1, omega = omega1)
  b0 <- rho_bounds(pi = pi0, omega = omega0)
  c(max(b1$rho_lower, b0$rho_lower), min(b1$rho_upper, b0$rho_upper))
}

# =============================================================================
# Table S1: sensitivity grid for the cytisine application (AE-M4)
#
# This is the grid that was Table 2 in the submitted version. The main text now
# reports only the design assumptions the original trial used; the wider grid
# lives here, with the NRI-scale quantities added (R2-M2a, R2-m4) and each
# correlation also reported as the dropout response probability gamma (AE-M6).
# =============================================================================
S1_PI1 <- rep(c(0.42, 0.47, 0.52), each = 3)
S1_PI0 <- rep(c(0.36, 0.41, 0.46), each = 3)
S1_OMEGA <- rep(c(0.05, 0.10, 0.15), times = 3)

tableS1_df <- timed("Table S1 (application sensitivity)", {
  rows <- lapply(seq_along(S1_PI1), function(i) {
    pi1_i <- S1_PI1[i]
    pi0_i <- S1_PI0[i]
    om_i <- S1_OMEGA[i]
    rng <- common_rho(pi1_i, pi0_i, om_i)

    N_cc <- sample_size_cc(
      pi1 = pi1_i, pi0 = pi0_i, omega1 = om_i, omega0 = om_i,
      r = 1, alpha = ALPHA, target_power = TARGET_POWER
    )$n_total_adjust

    sub <- lapply(seq_len(3), function(k) {
      rho_k <- c(rng[1], 0, rng[2])[k]
      ss <- sample_size_nri(
        pi1 = pi1_i, pi0 = pi0_i, omega1 = om_i, omega0 = om_i,
        rho1 = rho_k, rho0 = rho_k,
        r = 1, alpha = ALPHA, target_power = TARGET_POWER
      )
      data.frame(
        pi1 = pi1_i, pi0 = pi0_i, omega = om_i,
        rho = rho_k, rho_label = c("L", "0", "U")[k],
        gamma1 = rho_to_gamma(pi = pi1_i, omega = om_i, rho = rho_k),
        gamma0 = rho_to_gamma(pi = pi0_i, omega = om_i, rho = rho_k),
        pi1_NRI = ss$pi1_NRI, pi0_NRI = ss$pi0_NRI,
        effect_full = ss$effect_full, effect_NRI = ss$effect_NRI,
        N_cc = N_cc, N_nri = ss$n_total, RE = N_cc / ss$n_total,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, sub)
  })
  do.call(rbind, rows)
})
saveRDS(list(df = tableS1_df, alpha = ALPHA, target_power = TARGET_POWER, r = 1),
        file.path(data_dir, "tableS1_data.rds"))

# =============================================================================
# Figure S1: unequal dropout probabilities under the alternative (R2-M2b)
#
# The sample size formulas were derived allowing omega1 to differ from omega0,
# but the submitted numerical study always set them equal. Here omega1 is held
# at three values and omega0 is swept.
# =============================================================================
S2_OMEGA1 <- c(0.05, 0.10, 0.20)
S2_OMEGA0 <- seq(0.02, 0.35, by = 0.01)

figS1_df <- timed("Figure S1 (unequal dropout, alternative)", {
  grid <- expand.grid(omega1 = S2_OMEGA1, omega0 = S2_OMEGA0)
  rows <- lapply(seq_len(nrow(grid)), function(k) {
    om1 <- grid$omega1[k]
    om0 <- grid$omega0[k]

    ss_nri <- sample_size_nri(
      pi1 = REF_PI1, pi0 = REF_PI0, omega1 = om1, omega0 = om0,
      rho1 = 0, rho0 = 0, r = 1, alpha = ALPHA, target_power = TARGET_POWER
    )
    ss_cc <- sample_size_cc(
      pi1 = REF_PI1, pi0 = REF_PI0, omega1 = om1, omega0 = om0,
      r = 1, alpha = ALPHA, target_power = TARGET_POWER
    )

    n1_nri <- ss_nri$n1
    n0_nri <- ss_nri$n0
    n1_cc <- ss_cc$n1_adjust
    n0_cc <- ss_cc$n0_adjust

    ep <- power_nri(n1 = n1_nri, n0 = n0_nri, pi1 = REF_PI1, pi0 = REF_PI0,
                    omega1 = om1, omega0 = om0, alpha = ALPHA)
    es <- power_nri(n1 = n1_cc, n0 = n0_cc, pi1 = REF_PI1, pi0 = REF_PI0,
                    omega1 = om1, omega0 = om0, alpha = ALPHA)

    data.frame(
      omega1 = om1, omega0 = om0,
      N_proposed = ss_nri$n_total, N_Simple = ss_cc$n_total_adjust,
      power_proposed = ep$power, power_Simple = es$power,
      effect_NRI = ep$effect_NRI, effect_full = ep$effect_full,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
})
saveRDS(
  list(df = figS1_df, pi1 = REF_PI1, pi0 = REF_PI0, omega1_values = S2_OMEGA1,
       alpha = ALPHA, target_power = TARGET_POWER),
  file.path(data_dir, "figS1_data.rds")
)

# =============================================================================
# Figure S2: unequal dropout probabilities under the full-data null (R2-M2b)
#
# With pi1 = pi0 the full-data effect is zero, but the NRI effect is not unless
# the dropout probabilities are also equal. The rejection probability of the
# NRI test therefore departs from alpha as the dropout probabilities separate.
# This is not a failure of size control: it is a correctly sized test of a
# different null hypothesis, and the figure exists to make that explicit.
# =============================================================================
S3_PI <- c(0.30, 0.50, 0.70)
S3_OMEGA1 <- 0.10
S3_OMEGA0 <- seq(0.02, 0.35, by = 0.01)
S3_N <- 250   # per group, fixed, so that only the estimand changes

figS2_df <- timed("Figure S2 (unequal dropout, full-data null)", {
  grid <- expand.grid(pi = S3_PI, omega0 = S3_OMEGA0)
  rows <- lapply(seq_len(nrow(grid)), function(k) {
    pi_k <- grid$pi[k]
    om0 <- grid$omega0[k]
    res <- power_nri(
      n1 = S3_N, n0 = S3_N, pi1 = pi_k, pi0 = pi_k,
      omega1 = S3_OMEGA1, omega0 = om0, alpha = ALPHA
    )
    data.frame(
      pi = pi_k, omega1 = S3_OMEGA1, omega0 = om0, n_per_group = S3_N,
      effect_NRI = res$effect_NRI,
      # Rejection probability at the true parameter values, which under
      # pi1 = pi0 is the rejection probability under the full-data null
      reject_full_null = res$power,
      # Size of the same test against the NRI null, for contrast
      size_nri_null = res$type1_error,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
})
saveRDS(
  list(df = figS2_df, pi_values = S3_PI, omega1 = S3_OMEGA1,
       n_per_group = S3_N, alpha = ALPHA),
  file.path(data_dir, "figS2_data.rds")
)

# =============================================================================
# Figure S3: unequal correlations across groups (R1-M3a)
#
# Relative efficiency over the rectangle of correlations feasible in the two
# groups. The two correlations are separate parameters of separate groups, so
# the whole rectangle is attainable.
# =============================================================================
figS3_df <- timed("Figure S3 (unequal correlations)", {
  b1 <- rho_bounds(pi = REF_PI1, omega = REF_OMEGA)
  b0 <- rho_bounds(pi = REF_PI0, omega = REF_OMEGA)
  rho1_seq <- seq(b1$rho_lower, b1$rho_upper, length.out = 61)
  rho0_seq <- seq(b0$rho_lower, b0$rho_upper, length.out = 61)
  grid <- expand.grid(rho1 = rho1_seq, rho0 = rho0_seq)

  N_cc <- sample_size_cc(
    pi1 = REF_PI1, pi0 = REF_PI0, omega1 = REF_OMEGA, omega0 = REF_OMEGA,
    r = 1, alpha = ALPHA, target_power = TARGET_POWER
  )$n_total_adjust

  n_nri <- numeric(nrow(grid))
  eff <- numeric(nrow(grid))
  for (k in seq_len(nrow(grid))) {
    res <- tryCatch(
      sample_size_nri(
        pi1 = REF_PI1, pi0 = REF_PI0, omega1 = REF_OMEGA, omega0 = REF_OMEGA,
        rho1 = grid$rho1[k], rho0 = grid$rho0[k],
        r = 1, alpha = ALPHA, target_power = TARGET_POWER
      ),
      error = function(e) NULL
    )
    n_nri[k] <- if (is.null(res)) NA_real_ else res$n_total
    eff[k] <- if (is.null(res)) NA_real_ else res$effect_NRI
  }
  data.frame(
    rho1 = grid$rho1, rho0 = grid$rho0,
    N_cc = N_cc, N_nri = n_nri, re = N_cc / n_nri, effect_NRI = eff,
    stringsAsFactors = FALSE
  )
})
# The rectangle above contains points with rho1 different from rho0, and the
# region where the relative efficiency exceeds one lies entirely among them.
# The diagonal is computed separately so that the caption can state what
# happens when the two groups share a correlation, which is the assumption the
# main text works under.
figS3_diag <- timed("Figure S3 diagonal (rho1 = rho0)", {
  rng <- common_rho(REF_PI1, REF_PI0, REF_OMEGA)
  rho_seq <- seq(rng[1], rng[2], length.out = 61)
  N_cc <- sample_size_cc(
    pi1 = REF_PI1, pi0 = REF_PI0, omega1 = REF_OMEGA, omega0 = REF_OMEGA,
    r = 1, alpha = ALPHA, target_power = TARGET_POWER
  )$n_total_adjust
  rows <- lapply(rho_seq, function(rr) {
    ss <- sample_size_nri(
      pi1 = REF_PI1, pi0 = REF_PI0, omega1 = REF_OMEGA, omega0 = REF_OMEGA,
      rho1 = rr, rho0 = rr, r = 1, alpha = ALPHA, target_power = TARGET_POWER
    )
    data.frame(rho = rr, N_cc = N_cc, N_nri = ss$n_total,
               re = N_cc / ss$n_total, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
})

saveRDS(
  list(df = figS3_df, diagonal = figS3_diag,
       pi1 = REF_PI1, pi0 = REF_PI0, omega = REF_OMEGA,
       alpha = ALPHA, target_power = TARGET_POWER),
  file.path(data_dir, "figS3_data.rds")
)

# =============================================================================
# Figures S4 and S5: where the NRI effect exceeds the full-data effect (R2-M1)
#
# Reviewer 2 asked under what circumstances the NRI estimand is smaller than
# the full-data estimand, and whether it can ever be larger. Two maps answer
# this: one over the two dropout probabilities at independence, one over the
# two correlations at a common dropout probability. The zero contour of
# delta_NRI - delta_Full separates the two regimes. Both are stored in one
# file, distinguished by the panel column.
# =============================================================================
figS45_df <- timed("Figures S4 and S5 (effect ordering)", {
  delta_full <- REF_PI1 - REF_PI0

  om_seq <- seq(0.01, 0.40, by = 0.005)
  panel_a <- expand.grid(x = om_seq, y = om_seq)
  # At rho = 0 the correlation terms drop out of delta_NRI
  panel_a$effect_NRI <- REF_PI1 * (1 - panel_a$x) - REF_PI0 * (1 - panel_a$y)
  panel_a$diff <- panel_a$effect_NRI - delta_full
  panel_a$panel <- "dropout"

  b1 <- rho_bounds(pi = REF_PI1, omega = REF_OMEGA)
  b0 <- rho_bounds(pi = REF_PI0, omega = REF_OMEGA)
  panel_b <- expand.grid(
    x = seq(b1$rho_lower, b1$rho_upper, length.out = 81),
    y = seq(b0$rho_lower, b0$rho_upper, length.out = 81)
  )
  s1b <- sqrt(REF_PI1 * (1 - REF_PI1) * REF_OMEGA * (1 - REF_OMEGA))
  s0b <- sqrt(REF_PI0 * (1 - REF_PI0) * REF_OMEGA * (1 - REF_OMEGA))
  panel_b$effect_NRI <- REF_PI1 * (1 - REF_OMEGA) - REF_PI0 * (1 - REF_OMEGA) -
    panel_b$x * s1b + panel_b$y * s0b
  panel_b$diff <- panel_b$effect_NRI - delta_full
  panel_b$panel <- "correlation"

  out <- rbind(panel_a, panel_b)
  out$delta_full <- delta_full
  out
})
saveRDS(
  list(df = figS45_df, pi1 = REF_PI1, pi0 = REF_PI0, omega = REF_OMEGA),
  file.path(data_dir, "figS45_data.rds")
)

# =============================================================================
# Figure S6: conservative dropout assumptions (R1-M4, AE-m15)
#
# A practical alternative to a model-consistent sample size is to keep the
# simple inflation method and assume a dropout probability larger than the one
# expected. The figure reports the assumption needed to reach the NRI sample
# size, and the cost of that assumption when the true dropout probability turns
# out to be the expected one.
# =============================================================================
S5_OMEGA <- seq(0.05, 0.30, by = 0.01)

figS6_df <- timed("Figure S6 (conservative dropout)", {
  n_full <- sample_size_cc(
    pi1 = REF_PI1, pi0 = REF_PI0, r = 1, alpha = ALPHA,
    target_power = TARGET_POWER
  )$n_total

  rows <- lapply(S5_OMEGA, function(omega) {
    n_nri <- sample_size_nri(
      pi1 = REF_PI1, pi0 = REF_PI0, omega1 = omega, omega0 = omega,
      rho1 = 0, rho0 = 0, r = 1, alpha = ALPHA, target_power = TARGET_POWER
    )$n_total
    n_cc <- sample_size_cc(
      pi1 = REF_PI1, pi0 = REF_PI0, omega1 = omega, omega0 = omega,
      r = 1, alpha = ALPHA, target_power = TARGET_POWER
    )$n_total_adjust

    # Smallest assumed dropout probability whose inflated sample size reaches
    # the NRI requirement
    cand <- seq(omega, 0.95, by = 0.0005)
    inflated <- vapply(cand, function(w) {
      2 * ceiling((n_full / 2) / (1 - w))
    }, numeric(1))
    hit <- which(inflated >= n_nri)
    omega_star <- if (length(hit) == 0) NA_real_ else cand[hit[1]]

    data.frame(
      omega = omega,
      n_full = n_full,
      n_cc = n_cc,
      n_nri = n_nri,
      omega_star = omega_star,
      inflation_gap = omega_star - omega,
      excess_if_true_omega = n_nri - n_cc,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
})
saveRDS(
  list(df = figS6_df, pi1 = REF_PI1, pi0 = REF_PI0,
       alpha = ALPHA, target_power = TARGET_POWER),
  file.path(data_dir, "figS6_data.rds")
)

# =============================================================================
# Table S2: allocation ratios other than one (R1-M3b)
# =============================================================================
S6_R <- c(0.5, 1, 2, 3)
S6_SETTINGS <- list(
  list(pi1 = 0.45, pi0 = 0.30),
  list(pi1 = 0.55, pi0 = 0.40)
)
S6_OMEGA <- c(0.10, 0.20)

tableS2_df <- timed("Table S2 (allocation ratio)", {
  rows <- list()
  k <- 1L
  for (s in S6_SETTINGS) {
    for (omega in S6_OMEGA) {
      for (r in S6_R) {
        ss_cc <- sample_size_cc(
          pi1 = s$pi1, pi0 = s$pi0, omega1 = omega, omega0 = omega,
          r = r, alpha = ALPHA, target_power = TARGET_POWER
        )
        ss_nri <- sample_size_nri(
          pi1 = s$pi1, pi0 = s$pi0, omega1 = omega, omega0 = omega,
          rho1 = 0, rho0 = 0, r = r, alpha = ALPHA,
          target_power = TARGET_POWER
        )
        rows[[k]] <- data.frame(
          pi1 = s$pi1, pi0 = s$pi0, omega = omega, r = r,
          N_cc = ss_cc$n_total_adjust, N_nri = ss_nri$n_total,
          RE = ss_cc$n_total_adjust / ss_nri$n_total,
          stringsAsFactors = FALSE
        )
        k <- k + 1L
      }
    }
  }
  do.call(rbind, rows)
})
saveRDS(list(df = tableS2_df, alpha = ALPHA, target_power = TARGET_POWER),
        file.path(data_dir, "tableS2_data.rds"))

# =============================================================================
# Table S3: type I error rate at several null nuisance parameters (R2-m1)
#
# The size of a two-sample binomial test depends on the common response
# probability assumed under the null. The submitted version used the
# sample-size weighted average; this table reports the size across a range of
# null values, for both analyses, at the sample size the NRI formula gives.
# =============================================================================
tableS3_df <- timed("Table S3 (null nuisance parameter)", {
  ss <- sample_size_nri(
    pi1 = REF_PI1, pi0 = REF_PI0, omega1 = REF_OMEGA, omega0 = REF_OMEGA,
    rho1 = 0, rho0 = 0, r = 1, alpha = ALPHA, target_power = TARGET_POWER
  )
  n1 <- ss$n1
  n0 <- ss$n0
  pooled_nri <- (n1 * ss$pi1_NRI + n0 * ss$pi0_NRI) / (n1 + n0)
  pooled_cc <- (n1 * REF_PI1 + n0 * REF_PI0) / (n1 + n0)

  nri_nulls <- sort(unique(c(round(pooled_nri, 6), round(ss$pi0_NRI, 6),
                             round(ss$pi1_NRI, 6), 0.10, 0.20, 0.30, 0.50)))
  cc_nulls <- sort(unique(c(round(pooled_cc, 6), REF_PI0, REF_PI1,
                            0.10, 0.20, 0.50)))

  nri_rows <- lapply(nri_nulls, function(pi_bar) {
    res <- power_nri(n1 = n1, n0 = n0, pi1 = REF_PI1, pi0 = REF_PI0,
                     omega1 = REF_OMEGA, omega0 = REF_OMEGA,
                     alpha = ALPHA, pi_null = pi_bar)
    data.frame(analysis = "NRI", pi_null = pi_bar,
               is_default = isTRUE(all.equal(pi_bar, round(pooled_nri, 6))),
               type1_error = res$type1_error, stringsAsFactors = FALSE)
  })
  cc_rows <- lapply(cc_nulls, function(pi_bar) {
    res <- power_cc(n1 = n1, n0 = n0, pi1 = REF_PI1, pi0 = REF_PI0,
                    omega1 = REF_OMEGA, omega0 = REF_OMEGA,
                    alpha = ALPHA, pi_null = pi_bar)
    data.frame(analysis = "CC", pi_null = pi_bar,
               is_default = isTRUE(all.equal(pi_bar, round(pooled_cc, 6))),
               type1_error = res$type1_error, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, c(nri_rows, cc_rows))
  attr(out, "n1") <- n1
  attr(out, "n0") <- n0
  out
})
saveRDS(
  list(df = tableS3_df,
       n1 = attr(tableS3_df, "n1"), n0 = attr(tableS3_df, "n0"),
       pi1 = REF_PI1, pi0 = REF_PI0, omega = REF_OMEGA, alpha = ALPHA),
  file.path(data_dir, "tableS3_data.rds")
)

# =============================================================================
# Table S4: the three conditions are sufficient but not necessary (R2-m3, AE-M2)
#
# The submitted version moved between "if" and "only if" when describing when
# the simple inflation method is valid. A single counterexample settles the
# question. The dropout probabilities are made as unequal as is plausible, so
# condition (iii) fails, while the analysis is still CC and dropout is still
# independent of the latent response in each group. Inflating each group by its
# own dropout probability then delivers the target power exactly, which shows
# that equal dropout probabilities are not necessary.
#
# This also answers the associate editor's question about whether MCAR already
# implies equal dropout probabilities. It does not: dropout that depends on the
# randomized group but on nothing else is independent of the outcome within
# each group, and a CC analysis remains valid, yet omega1 differs from omega0.
# =============================================================================
S4_OMEGA1 <- c(0.10, 0.05, 0.20, 0.02, 0.35, 0.40)
S4_OMEGA0 <- c(0.10, 0.20, 0.05, 0.35, 0.02, 0.05)

tableS4_df <- timed("Table S4 (unequal dropout counterexample)", {
  rows <- lapply(seq_along(S4_OMEGA1), function(i) {
    w1 <- S4_OMEGA1[i]
    w0 <- S4_OMEGA0[i]
    ss <- sample_size_cc(
      pi1 = REF_PI1, pi0 = REF_PI0, omega1 = w1, omega0 = w0,
      r = 1, alpha = ALPHA, target_power = TARGET_POWER
    )
    res <- power_cc(
      n1 = ss$n1_adjust, n0 = ss$n0_adjust,
      pi1 = REF_PI1, pi0 = REF_PI0, omega1 = w1, omega0 = w0,
      rho1 = 0, rho0 = 0, alpha = ALPHA
    )
    data.frame(
      omega1 = w1, omega0 = w0,
      n_complete = ss$n1,
      n1 = ss$n1_adjust, n0 = ss$n0_adjust,
      n_total = ss$n_total_adjust,
      power = res$power,
      type1_error = res$type1_error,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
})
saveRDS(
  list(df = tableS4_df, pi1 = REF_PI1, pi0 = REF_PI0,
       alpha = ALPHA, target_power = TARGET_POWER),
  file.path(data_dir, "tableS4_data.rds")
)

# =============================================================================
# Figure S7: what a historical rate does and does not pin down (R1-M2, R2-M3)
#
# The cytisine trial took its control response probability from a previous
# trial in the same population, which reported that 254 of 620 randomized
# patients, 41 percent, were abstinent at six months, with abstinence verified
# biochemically and everything unverified counted as a failure. That reported
# number is a joint probability of responding and completing, not the latent
# response probability, so the latent value is identified only to an interval
# of width omega. The interval's position depends on the counting rule the
# historical trial used, which is the point R1-M2 raised about the best and
# worst case rules and the point R2-M3 raised about mixing an analysis model
# with a data-generating model.
#
# The figure carries the whole interval through to the planning question. Each
# curve spans its own rule's identified interval, so the horizontal extent of a
# curve is the interval itself. The two panels differ in how the treatment
# group is tied to the control group, which the historical data cannot settle
# either: a common correlation, as in Table 2, or a common probability of
# responding among dropouts.
# =============================================================================
S7_PI_OBS <- 0.41
S7_OMEGA <- 0.10
S7_DELTA <- 0.06
S7_N_HIST <- 620
S7_N_DROP_HIST <- round(S7_N_HIST * S7_OMEGA)
S7_N_POINTS <- 101

S7_RULES <- data.frame(
  method = c("nri", "cc", "wc"),
  group  = c("control", "control", "control"),
  label  = c("Dropouts counted as non-responders",
             "Complete case",
             "Dropouts counted as responders"),
  stringsAsFactors = FALSE
)

figS7_df <- timed("Figure S7 (partial identification)", {
  out <- lapply(seq_len(nrow(S7_RULES)), function(i) {
    meth <- S7_RULES$method[i]
    grp  <- S7_RULES$group[i]
    lb <- latent_bounds(
      pi_obs = S7_PI_OBS, N_randomized = S7_N_HIST,
      N_dropout = S7_N_DROP_HIST, method = meth, group = grp
    )
    grid <- seq(lb$pi_lower, lb$pi_upper, length.out = S7_N_POINTS)
    grid <- grid[grid > 0 & grid < 1]

    rows <- lapply(grid, function(pi0) {
      # Clamped only against floating point drift at the two endpoints,
      # where gamma is exactly 0 and exactly 1 by construction
      gamma0 <- min(max((pi0 - lb$pi_nri) / lb$omega, 0), 1)
      rho0 <- infer_rho(pi_obs = S7_PI_OBS, omega = lb$omega, pi = pi0,
                        method = meth, group = grp)$rho
      pi1 <- pi0 + S7_DELTA
      n_cc <- sample_size_cc(
        pi1 = pi1, pi0 = pi0, omega1 = lb$omega, omega0 = lb$omega,
        r = 1, alpha = ALPHA, target_power = TARGET_POWER
      )$n_total_adjust

      # Two ways of tying the treatment group to the control group
      conv <- list(
        "Common correlation" = rho0,
        "Common dropout response probability" =
          gamma_to_rho(pi = pi1, omega = lb$omega, gamma = gamma0)
      )
      do.call(rbind, lapply(names(conv), function(cn) {
        rho1 <- conv[[cn]]
        b1 <- rho_bounds(pi = pi1, omega = lb$omega)
        feasible <- rho1 >= b1$rho_lower - 1e-9 && rho1 <= b1$rho_upper + 1e-9
        n_nri <- if (feasible) {
          sample_size_nri(
            pi1 = pi1, pi0 = pi0, omega1 = lb$omega, omega0 = lb$omega,
            rho1 = rho1, rho0 = rho0,
            r = 1, alpha = ALPHA, target_power = TARGET_POWER
          )$n_total
        } else {
          NA_integer_
        }
        data.frame(
          rule = S7_RULES$label[i],
          convention = cn,
          pi0 = pi0, gamma0 = gamma0, rho1 = rho1, rho0 = rho0,
          gamma1 = rho_to_gamma(pi = pi1, omega = lb$omega, rho = rho1),
          feasible = feasible,
          n_cc = n_cc, n_nri = n_nri, RE = n_cc / n_nri,
          stringsAsFactors = FALSE
        )
      }))
    })
    dat <- do.call(rbind, rows)
    dat$pi_lower <- lb$pi_lower
    dat$pi_upper <- lb$pi_upper
    dat
  })
  do.call(rbind, out)
})
saveRDS(
  list(df = figS7_df, pi_obs = S7_PI_OBS, omega = S7_OMEGA,
       delta = S7_DELTA, n_hist = S7_N_HIST, alpha = ALPHA,
       target_power = TARGET_POWER, rules = S7_RULES),
  file.path(data_dir, "figS7_data.rds")
)

# --------------------------------------------------------------------------- #
# Record the execution environment for reproducibility
# --------------------------------------------------------------------------- #
writeLines(
  capture.output(sessionInfo()),
  file.path(data_dir, "sessionInfo_data_generate_supplement.txt")
)

log_msg("All supplementary data written to '%s/'.", data_dir)

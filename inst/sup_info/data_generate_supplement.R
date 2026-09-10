# =============================================================================
# data_generate_supplement.R
#
# Computation stage for the supplementary material of:
#   "A Cautionary Note on Sample Size Inflation for Trials with
#    Non-Responder Imputation"
#
# Every display answers a specific referee comment:
#
#   Table S1   sensitivity grid over p1 and omega, moved out of the main text
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
REF_P1 <- 0.45
REF_P2 <- 0.30
REF_OMEGA <- 0.10

# Feasible correlation range in both groups at a given pair of probabilities
common_rho <- function(p1, p2, omega1, omega2 = omega1) {
  b1 <- rho_bounds(p = p1, omega = omega1)
  b2 <- rho_bounds(p = p2, omega = omega2)
  c(max(b1$rho_lower, b2$rho_lower), min(b1$rho_upper, b2$rho_upper))
}

# =============================================================================
# Table S1: sensitivity grid for the cytisine application (AE-M4)
#
# This is the grid that was Table 2 in the submitted version. The main text now
# reports only the design assumptions the original trial used; the wider grid
# lives here, with the NRI-scale quantities added (R2-M2a, R2-m4) and each
# correlation also reported as the dropout response probability gamma (AE-M6).
# =============================================================================
S1_P1 <- rep(c(0.42, 0.47, 0.52), each = 3)
S1_P2 <- rep(c(0.36, 0.41, 0.46), each = 3)
S1_OMEGA <- rep(c(0.05, 0.10, 0.15), times = 3)

tableS1_df <- timed("Table S1 (application sensitivity)", {
  rows <- lapply(seq_along(S1_P1), function(i) {
    p1_i <- S1_P1[i]
    p2_i <- S1_P2[i]
    om_i <- S1_OMEGA[i]
    rng <- common_rho(p1_i, p2_i, om_i)

    N_cc <- sample_size_cc(
      p1 = p1_i, p2 = p2_i, omega1 = om_i, omega2 = om_i,
      r = 1, alpha = ALPHA, target_power = TARGET_POWER
    )$n_total_adjust

    sub <- lapply(seq_len(3), function(k) {
      rho_k <- c(rng[1], 0, rng[2])[k]
      ss <- sample_size_nri(
        p1 = p1_i, p2 = p2_i, omega1 = om_i, omega2 = om_i,
        rho1 = rho_k, rho2 = rho_k,
        r = 1, alpha = ALPHA, target_power = TARGET_POWER
      )
      data.frame(
        p1 = p1_i, p2 = p2_i, omega = om_i,
        rho = rho_k, rho_label = c("L", "0", "U")[k],
        gamma1 = rho_to_gamma(p = p1_i, omega = om_i, rho = rho_k),
        gamma2 = rho_to_gamma(p = p2_i, omega = om_i, rho = rho_k),
        p1_NRI = ss$p1_NRI, p2_NRI = ss$p2_NRI,
        effect_latent = ss$effect_latent, effect_NRI = ss$effect_NRI,
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
# The sample size formulas were derived allowing omega1 to differ from omega2,
# but the submitted numerical study always set them equal. Here omega1 is held
# at three values and omega2 is swept.
# =============================================================================
S2_OMEGA1 <- c(0.05, 0.10, 0.20)
S2_OMEGA2 <- seq(0.02, 0.35, by = 0.01)

figS1_df <- timed("Figure S1 (unequal dropout, alternative)", {
  grid <- expand.grid(omega1 = S2_OMEGA1, omega2 = S2_OMEGA2)
  rows <- lapply(seq_len(nrow(grid)), function(k) {
    om1 <- grid$omega1[k]
    om2 <- grid$omega2[k]

    ss_nri <- sample_size_nri(
      p1 = REF_P1, p2 = REF_P2, omega1 = om1, omega2 = om2,
      rho1 = 0, rho2 = 0, r = 1, alpha = ALPHA, target_power = TARGET_POWER
    )
    ss_cc <- sample_size_cc(
      p1 = REF_P1, p2 = REF_P2, omega1 = om1, omega2 = om2,
      r = 1, alpha = ALPHA, target_power = TARGET_POWER
    )

    n1_nri <- ss_nri$n1
    n2_nri <- ss_nri$n2
    n1_cc <- ss_cc$n1_adjust
    n2_cc <- ss_cc$n2_adjust

    ep <- power_nri(n1 = n1_nri, n2 = n2_nri, p1 = REF_P1, p2 = REF_P2,
                    omega1 = om1, omega2 = om2, alpha = ALPHA)
    es <- power_nri(n1 = n1_cc, n2 = n2_cc, p1 = REF_P1, p2 = REF_P2,
                    omega1 = om1, omega2 = om2, alpha = ALPHA)

    data.frame(
      omega1 = om1, omega2 = om2,
      N_proposed = ss_nri$n_total, N_Simple = ss_cc$n_total_adjust,
      power_proposed = ep$power, power_Simple = es$power,
      effect_NRI = ep$effect_NRI, effect_latent = ep$effect_latent,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
})
saveRDS(
  list(df = figS1_df, p1 = REF_P1, p2 = REF_P2, omega1_values = S2_OMEGA1,
       alpha = ALPHA, target_power = TARGET_POWER),
  file.path(data_dir, "figS1_data.rds")
)

# =============================================================================
# Figure S2: unequal dropout probabilities under the full-data null (R2-M2b)
#
# With p1 = p2 the full-data effect is zero, but the NRI effect is not unless
# the dropout probabilities are also equal. The rejection probability of the
# NRI test therefore departs from alpha as the dropout probabilities separate.
# This is not a failure of size control: it is a correctly sized test of a
# different null hypothesis, and the figure exists to make that explicit.
# =============================================================================
S3_P <- c(0.30, 0.50, 0.70)
S3_OMEGA1 <- 0.10
S3_OMEGA2 <- seq(0.02, 0.35, by = 0.01)
S3_N <- 250   # per group, fixed, so that only the estimand changes

figS2_df <- timed("Figure S2 (unequal dropout, full-data null)", {
  grid <- expand.grid(p = S3_P, omega2 = S3_OMEGA2)
  rows <- lapply(seq_len(nrow(grid)), function(k) {
    p_k <- grid$p[k]
    om2 <- grid$omega2[k]
    res <- power_nri(
      n1 = S3_N, n2 = S3_N, p1 = p_k, p2 = p_k,
      omega1 = S3_OMEGA1, omega2 = om2, alpha = ALPHA
    )
    data.frame(
      p = p_k, omega1 = S3_OMEGA1, omega2 = om2, n_per_group = S3_N,
      effect_NRI = res$effect_NRI,
      # Rejection probability at the true parameter values, which under
      # p1 = p2 is the rejection probability under the full-data null
      reject_full_null = res$power,
      # Size of the same test against the NRI null, for contrast
      size_nri_null = res$type1_error,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
})
saveRDS(
  list(df = figS2_df, p_values = S3_P, omega1 = S3_OMEGA1,
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
  b1 <- rho_bounds(p = REF_P1, omega = REF_OMEGA)
  b2 <- rho_bounds(p = REF_P2, omega = REF_OMEGA)
  rho1_seq <- seq(b1$rho_lower, b1$rho_upper, length.out = 61)
  rho2_seq <- seq(b2$rho_lower, b2$rho_upper, length.out = 61)
  grid <- expand.grid(rho1 = rho1_seq, rho2 = rho2_seq)

  N_cc <- sample_size_cc(
    p1 = REF_P1, p2 = REF_P2, omega1 = REF_OMEGA, omega2 = REF_OMEGA,
    r = 1, alpha = ALPHA, target_power = TARGET_POWER
  )$n_total_adjust

  n_nri <- numeric(nrow(grid))
  eff <- numeric(nrow(grid))
  for (k in seq_len(nrow(grid))) {
    res <- tryCatch(
      sample_size_nri(
        p1 = REF_P1, p2 = REF_P2, omega1 = REF_OMEGA, omega2 = REF_OMEGA,
        rho1 = grid$rho1[k], rho2 = grid$rho2[k],
        r = 1, alpha = ALPHA, target_power = TARGET_POWER
      ),
      error = function(e) NULL
    )
    n_nri[k] <- if (is.null(res)) NA_real_ else res$n_total
    eff[k] <- if (is.null(res)) NA_real_ else res$effect_NRI
  }
  data.frame(
    rho1 = grid$rho1, rho2 = grid$rho2,
    N_cc = N_cc, N_nri = n_nri, re = N_cc / n_nri, effect_NRI = eff,
    stringsAsFactors = FALSE
  )
})
# The rectangle above contains points with rho1 different from rho2, and the
# region where the relative efficiency exceeds one lies entirely among them.
# The diagonal is computed separately so that the caption can state what
# happens when the two groups share a correlation, which is the assumption the
# main text works under.
figS3_diag <- timed("Figure S3 diagonal (rho1 = rho2)", {
  rng <- common_rho(REF_P1, REF_P2, REF_OMEGA)
  rho_seq <- seq(rng[1], rng[2], length.out = 61)
  N_cc <- sample_size_cc(
    p1 = REF_P1, p2 = REF_P2, omega1 = REF_OMEGA, omega2 = REF_OMEGA,
    r = 1, alpha = ALPHA, target_power = TARGET_POWER
  )$n_total_adjust
  rows <- lapply(rho_seq, function(rr) {
    ss <- sample_size_nri(
      p1 = REF_P1, p2 = REF_P2, omega1 = REF_OMEGA, omega2 = REF_OMEGA,
      rho1 = rr, rho2 = rr, r = 1, alpha = ALPHA, target_power = TARGET_POWER
    )
    data.frame(rho = rr, N_cc = N_cc, N_nri = ss$n_total,
               re = N_cc / ss$n_total, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
})

saveRDS(
  list(df = figS3_df, diagonal = figS3_diag,
       p1 = REF_P1, p2 = REF_P2, omega = REF_OMEGA,
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
# delta_NRI - delta_Trad separates the two regimes. Both are stored in one
# file, distinguished by the panel column.
# =============================================================================
figS45_df <- timed("Figures S4 and S5 (effect ordering)", {
  delta_latent <- REF_P1 - REF_P2

  om_seq <- seq(0.01, 0.40, by = 0.005)
  panel_a <- expand.grid(x = om_seq, y = om_seq)
  # At rho = 0 the correlation terms drop out of delta_NRI
  panel_a$effect_NRI <- REF_P1 * (1 - panel_a$x) - REF_P2 * (1 - panel_a$y)
  panel_a$diff <- panel_a$effect_NRI - delta_latent
  panel_a$panel <- "dropout"

  b1 <- rho_bounds(p = REF_P1, omega = REF_OMEGA)
  b2 <- rho_bounds(p = REF_P2, omega = REF_OMEGA)
  panel_b <- expand.grid(
    x = seq(b1$rho_lower, b1$rho_upper, length.out = 81),
    y = seq(b2$rho_lower, b2$rho_upper, length.out = 81)
  )
  s1b <- sqrt(REF_P1 * (1 - REF_P1) * REF_OMEGA * (1 - REF_OMEGA))
  s2b <- sqrt(REF_P2 * (1 - REF_P2) * REF_OMEGA * (1 - REF_OMEGA))
  panel_b$effect_NRI <- REF_P1 * (1 - REF_OMEGA) - REF_P2 * (1 - REF_OMEGA) -
    panel_b$x * s1b + panel_b$y * s2b
  panel_b$diff <- panel_b$effect_NRI - delta_latent
  panel_b$panel <- "correlation"

  out <- rbind(panel_a, panel_b)
  out$delta_latent <- delta_latent
  out
})
saveRDS(
  list(df = figS45_df, p1 = REF_P1, p2 = REF_P2, omega = REF_OMEGA),
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
  n_trad <- sample_size_cc(
    p1 = REF_P1, p2 = REF_P2, r = 1, alpha = ALPHA,
    target_power = TARGET_POWER
  )$n_total

  rows <- lapply(S5_OMEGA, function(omega) {
    n_nri <- sample_size_nri(
      p1 = REF_P1, p2 = REF_P2, omega1 = omega, omega2 = omega,
      rho1 = 0, rho2 = 0, r = 1, alpha = ALPHA, target_power = TARGET_POWER
    )$n_total
    n_cc <- sample_size_cc(
      p1 = REF_P1, p2 = REF_P2, omega1 = omega, omega2 = omega,
      r = 1, alpha = ALPHA, target_power = TARGET_POWER
    )$n_total_adjust

    # Smallest assumed dropout probability whose inflated sample size reaches
    # the NRI requirement
    cand <- seq(omega, 0.95, by = 0.0005)
    inflated <- vapply(cand, function(w) {
      2 * ceiling((n_trad / 2) / (1 - w))
    }, numeric(1))
    hit <- which(inflated >= n_nri)
    omega_star <- if (length(hit) == 0) NA_real_ else cand[hit[1]]

    data.frame(
      omega = omega,
      n_trad = n_trad,
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
  list(df = figS6_df, p1 = REF_P1, p2 = REF_P2,
       alpha = ALPHA, target_power = TARGET_POWER),
  file.path(data_dir, "figS6_data.rds")
)

# =============================================================================
# Table S2: allocation ratios other than one (R1-M3b)
# =============================================================================
S6_R <- c(0.5, 1, 2, 3)
S6_SETTINGS <- list(
  list(p1 = 0.45, p2 = 0.30),
  list(p1 = 0.55, p2 = 0.40)
)
S6_OMEGA <- c(0.10, 0.20)

tableS2_df <- timed("Table S2 (allocation ratio)", {
  rows <- list()
  k <- 1L
  for (s in S6_SETTINGS) {
    for (omega in S6_OMEGA) {
      for (r in S6_R) {
        ss_cc <- sample_size_cc(
          p1 = s$p1, p2 = s$p2, omega1 = omega, omega2 = omega,
          r = r, alpha = ALPHA, target_power = TARGET_POWER
        )
        ss_nri <- sample_size_nri(
          p1 = s$p1, p2 = s$p2, omega1 = omega, omega2 = omega,
          rho1 = 0, rho2 = 0, r = r, alpha = ALPHA,
          target_power = TARGET_POWER
        )
        rows[[k]] <- data.frame(
          p1 = s$p1, p2 = s$p2, omega = omega, r = r,
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
    p1 = REF_P1, p2 = REF_P2, omega1 = REF_OMEGA, omega2 = REF_OMEGA,
    rho1 = 0, rho2 = 0, r = 1, alpha = ALPHA, target_power = TARGET_POWER
  )
  n1 <- ss$n1
  n2 <- ss$n2
  pooled_nri <- (n1 * ss$p1_NRI + n2 * ss$p2_NRI) / (n1 + n2)
  pooled_cc <- (n1 * REF_P1 + n2 * REF_P2) / (n1 + n2)

  nri_nulls <- sort(unique(c(round(pooled_nri, 6), round(ss$p2_NRI, 6),
                             round(ss$p1_NRI, 6), 0.10, 0.20, 0.30, 0.50)))
  cc_nulls <- sort(unique(c(round(pooled_cc, 6), REF_P2, REF_P1,
                            0.10, 0.20, 0.50)))

  nri_rows <- lapply(nri_nulls, function(p0) {
    res <- power_nri(n1 = n1, n2 = n2, p1 = REF_P1, p2 = REF_P2,
                     omega1 = REF_OMEGA, omega2 = REF_OMEGA,
                     alpha = ALPHA, p_null = p0)
    data.frame(analysis = "NRI", p_null = p0,
               is_default = isTRUE(all.equal(p0, round(pooled_nri, 6))),
               type1_error = res$type1_error, stringsAsFactors = FALSE)
  })
  cc_rows <- lapply(cc_nulls, function(p0) {
    res <- power_cc(n1 = n1, n2 = n2, p1 = REF_P1, p2 = REF_P2,
                    omega1 = REF_OMEGA, omega2 = REF_OMEGA,
                    alpha = ALPHA, p_null = p0)
    data.frame(analysis = "CC", p_null = p0,
               is_default = isTRUE(all.equal(p0, round(pooled_cc, 6))),
               type1_error = res$type1_error, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, c(nri_rows, cc_rows))
  attr(out, "n1") <- n1
  attr(out, "n2") <- n2
  out
})
saveRDS(
  list(df = tableS3_df,
       n1 = attr(tableS3_df, "n1"), n2 = attr(tableS3_df, "n2"),
       p1 = REF_P1, p2 = REF_P2, omega = REF_OMEGA, alpha = ALPHA),
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
# each group, and a CC analysis remains valid, yet omega1 differs from omega2.
# =============================================================================
S4_OMEGA1 <- c(0.10, 0.05, 0.20, 0.02, 0.35, 0.40)
S4_OMEGA2 <- c(0.10, 0.20, 0.05, 0.35, 0.02, 0.05)

tableS4_df <- timed("Table S4 (unequal dropout counterexample)", {
  rows <- lapply(seq_along(S4_OMEGA1), function(i) {
    w1 <- S4_OMEGA1[i]
    w2 <- S4_OMEGA2[i]
    ss <- sample_size_cc(
      p1 = REF_P1, p2 = REF_P2, omega1 = w1, omega2 = w2,
      r = 1, alpha = ALPHA, target_power = TARGET_POWER
    )
    res <- power_cc(
      n1 = ss$n1_adjust, n2 = ss$n2_adjust,
      p1 = REF_P1, p2 = REF_P2, omega1 = w1, omega2 = w2,
      rho1 = 0, rho2 = 0, alpha = ALPHA
    )
    data.frame(
      omega1 = w1, omega2 = w2,
      n_complete = ss$n1,
      n1 = ss$n1_adjust, n2 = ss$n2_adjust,
      n_total = ss$n_total_adjust,
      power = res$power,
      type1_error = res$type1_error,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
})
saveRDS(
  list(df = tableS4_df, p1 = REF_P1, p2 = REF_P2,
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
S7_P_OBS <- 0.41
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
      p_obs = S7_P_OBS, N_randomized = S7_N_HIST,
      N_dropout = S7_N_DROP_HIST, method = meth, group = grp
    )
    grid <- seq(lb$p_lower, lb$p_upper, length.out = S7_N_POINTS)
    grid <- grid[grid > 0 & grid < 1]

    rows <- lapply(grid, function(p2) {
      # Clamped only against floating point drift at the two endpoints,
      # where gamma is exactly 0 and exactly 1 by construction
      gamma2 <- min(max((p2 - lb$p_nri) / lb$omega, 0), 1)
      rho2 <- infer_rho(p_obs = S7_P_OBS, omega = lb$omega, p = p2,
                        method = meth, group = grp)$rho
      p1 <- p2 + S7_DELTA
      n_cc <- sample_size_cc(
        p1 = p1, p2 = p2, omega1 = lb$omega, omega2 = lb$omega,
        r = 1, alpha = ALPHA, target_power = TARGET_POWER
      )$n_total_adjust

      # Two ways of tying the treatment group to the control group
      conv <- list(
        "Common correlation" = rho2,
        "Common dropout response probability" =
          gamma_to_rho(p = p1, omega = lb$omega, gamma = gamma2)
      )
      do.call(rbind, lapply(names(conv), function(cn) {
        rho1 <- conv[[cn]]
        b1 <- rho_bounds(p = p1, omega = lb$omega)
        feasible <- rho1 >= b1$rho_lower - 1e-9 && rho1 <= b1$rho_upper + 1e-9
        n_nri <- if (feasible) {
          sample_size_nri(
            p1 = p1, p2 = p2, omega1 = lb$omega, omega2 = lb$omega,
            rho1 = rho1, rho2 = rho2,
            r = 1, alpha = ALPHA, target_power = TARGET_POWER
          )$n_total
        } else {
          NA_integer_
        }
        data.frame(
          rule = S7_RULES$label[i],
          convention = cn,
          p2 = p2, gamma2 = gamma2, rho1 = rho1, rho2 = rho2,
          gamma1 = rho_to_gamma(p = p1, omega = lb$omega, rho = rho1),
          feasible = feasible,
          n_cc = n_cc, n_nri = n_nri, RE = n_cc / n_nri,
          stringsAsFactors = FALSE
        )
      }))
    })
    dat <- do.call(rbind, rows)
    dat$p_lower <- lb$p_lower
    dat$p_upper <- lb$p_upper
    dat
  })
  do.call(rbind, out)
})
saveRDS(
  list(df = figS7_df, p_obs = S7_P_OBS, omega = S7_OMEGA,
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

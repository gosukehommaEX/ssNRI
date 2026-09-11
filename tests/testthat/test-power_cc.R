# Reference values in ref_power_cc.csv were produced by an independent
# implementation written in Python (SciPy), using the full double enumeration
# over completer counts and responder counts. They therefore check the R
# implementation against a different language and a different summation order,
# not merely against itself. The script that produces them is
# inst/reference/generate_reference_values.py; it is coded from the formulas in
# the article and does not load this package.

ref <- utils::read.csv(
  testthat::test_path("ref_power_cc.csv"),
  stringsAsFactors = FALSE
)

test_that("power_cc reproduces an independent Python enumeration", {
  for (i in seq_len(nrow(ref))) {
    r <- ref[i, ]
    out <- power_cc(
      n1 = r$n1, n0 = r$n0, pi1 = r$pi1, pi0 = r$pi0,
      omega1 = r$omega1, omega0 = r$omega0,
      rho1 = r$rho1, rho0 = r$rho0, alpha = r$alpha
    )
    info <- sprintf("row %d: n = (%d, %d), rho = (%g, %g), omega = (%g, %g)",
                    i, r$n1, r$n0, r$rho1, r$rho0, r$omega1, r$omega0)
    expect_equal(out$pi1_CC, r$pi1_CC, tolerance = 1e-12, info = info)
    expect_equal(out$pi0_CC, r$pi0_CC, tolerance = 1e-12, info = info)
    expect_equal(out$power, r$power, tolerance = 1e-10, info = info)
    expect_equal(out$type1_error, r$type1_error, tolerance = 1e-10, info = info)
  }
})

test_that("power_cc agrees with a direct double enumeration", {
  # Deliberately naive reference: enumerate every (m1, m0, x1, x0) cell
  naive_power_cc <- function(n1, n0, pi1, pi0, omega1, omega0,
                             rho1 = 0, rho0 = 0, alpha = 0.025) {
    s1 <- sqrt(pi1 * (1 - pi1) * omega1 * (1 - omega1))
    s0 <- sqrt(pi0 * (1 - pi0) * omega0 * (1 - omega0))
    pi1_00 <- (1 - pi1) * (1 - omega1) + rho1 * s1
    pi1_10 <- pi1 * (1 - omega1) - rho1 * s1
    pi0_00 <- (1 - pi0) * (1 - omega0) + rho0 * s0
    pi0_10 <- pi0 * (1 - omega0) - rho0 * s0
    pi1_CC <- pi1_10 / (pi1_00 + pi1_10)
    pi0_CC <- pi0_10 / (pi0_00 + pi0_10)
    pi_bar_CC <- (n1 * pi1_CC + n0 * pi0_CC) / (n1 + n0)
    z <- stats::qnorm(1 - alpha)
    # the null hypothesis concerns the response probabilities only, so the
    # completer counts keep the group-specific dropout probabilities
    pm1 <- stats::dbinom(0:n1, n1, 1 - omega1)
    pm0 <- stats::dbinom(0:n0, n0, 1 - omega0)
    pw <- 0
    t1 <- 0
    for (m1 in 1:n1) {
      for (m0 in 1:n0) {
        x1 <- 0:m1
        x0 <- 0:m0
        pp <- outer(x1, x0, function(a, b) (a + b) / (m1 + m0))
        se <- sqrt(pp * (1 - pp) * (1 / m1 + 1 / m0))
        se[se < 1e-10] <- 1e-10
        rej <- (outer(x1 / m1, rep(1, m0 + 1)) -
                  outer(rep(1, m1 + 1), x0 / m0)) / se > z
        pw <- pw + pm1[m1 + 1] * pm0[m0 + 1] *
          sum(outer(stats::dbinom(x1, m1, pi1_CC),
                    stats::dbinom(x0, m0, pi0_CC)) * rej)
        t1 <- t1 + pm1[m1 + 1] * pm0[m0 + 1] *
          sum(outer(stats::dbinom(x1, m1, pi_bar_CC),
                    stats::dbinom(x0, m0, pi_bar_CC)) * rej)
      }
    }
    c(power = pw, type1_error = t1)
  }

  cfg <- list(n1 = 22, n0 = 30, pi1 = 0.55, pi0 = 0.30,
              omega1 = 0.18, omega0 = 0.08, rho1 = 0.15, rho0 = -0.10,
              alpha = 0.025)
  fast <- do.call(power_cc, cfg)
  slow <- do.call(naive_power_cc, cfg)

  expect_equal(fast$power, unname(slow["power"]), tolerance = 1e-12)
  expect_equal(fast$type1_error, unname(slow["type1_error"]), tolerance = 1e-12)
})

test_that("tol = 0 is exact and a positive tol is a lower bound within eps", {
  cfg <- list(n1 = 80, n0 = 80, pi1 = 0.45, pi0 = 0.25,
              omega1 = 0.15, omega0 = 0.15, rho1 = 0, rho0 = 0)
  exact <- do.call(power_cc, cfg)
  eps <- 1e-10
  approx <- do.call(power_cc, c(cfg, list(tol = eps)))
  n_cells <- (cfg$n1 + 1) * (cfg$n0 + 1)

  expect_lte(approx$power, exact$power + 1e-15)
  expect_gte(approx$power, exact$power - n_cells * eps)
  expect_lte(approx$type1_error, exact$type1_error + 1e-15)
  expect_gte(approx$type1_error, exact$type1_error - n_cells * eps)
})

test_that("rho = 0 gives p_CC equal to the latent response probability", {
  out <- power_cc(n1 = 40, n0 = 40, pi1 = 0.6, pi0 = 0.4,
                  omega1 = 0.25, omega0 = 0.10)
  expect_equal(out$pi1_CC, 0.6, tolerance = 1e-12)
  expect_equal(out$pi0_CC, 0.4, tolerance = 1e-12)
  expect_equal(out$effect_CC, 0.2, tolerance = 1e-12)
})

test_that("power_cc rejects invalid input", {
  expect_error(power_cc(0, 40, 0.6, 0.4, 0.1, 0.1), "positive")
  expect_error(power_cc(40, 40, 1.2, 0.4, 0.1, 0.1), "pi1 must be")
  expect_error(power_cc(40, 40, 0.6, 0.4, 0.1, 0.1, tol = -1), "non-negative")
  expect_error(power_cc(40, 40, 0.6, 0.4, 0.1, 0.1, rho1 = 5), "rho1 must be")
})

test_that("pi_null defaults to the pooled value and can be set explicitly", {
  cfg <- list(n1 = 40, n0 = 40, pi1 = 0.6, pi0 = 0.4,
              omega1 = 0.15, omega0 = 0.15)
  default <- do.call(power_cc, cfg)
  expect_equal(default$pi_null, 0.5, tolerance = 1e-12)

  explicit <- do.call(power_cc, c(cfg, list(pi_null = 0.5)))
  expect_equal(explicit$type1_error, default$type1_error, tolerance = 1e-14)
  expect_equal(explicit$power, default$power, tolerance = 1e-14)
})

test_that("pi_null changes the size but never the power", {
  # Reference values from the independent Python enumeration, transcribed to
  # ten decimal places, so the comparison is absolute at that precision rather
  # than relative; see the same test in test-power_nri.R for the reasoning.
  cfg <- list(n1 = 40, n0 = 40, pi1 = 0.6, pi0 = 0.4,
              omega1 = 0.15, omega0 = 0.15)
  expected <- c("0.5" = 0.0256208522,
                "0.4" = 0.0257491919,
                "0.25" = 0.0258309009)
  base_power <- do.call(power_cc, cfg)$power
  for (nm in names(expected)) {
    out <- do.call(power_cc, c(cfg, list(pi_null = as.numeric(nm))))
    expect_lt(abs(out$type1_error - unname(expected[nm])), 1e-9,
              label = paste("|type1_error - reference| at pi_null =", nm))
    expect_equal(out$power, base_power, tolerance = 1e-14, info = nm)
  }
})

test_that("power_cc rejects an invalid pi_null", {
  cfg <- list(n1 = 30, n0 = 30, pi1 = 0.6, pi0 = 0.4,
              omega1 = 0.1, omega0 = 0.1)
  expect_error(do.call(power_cc, c(cfg, list(pi_null = 0))), "pi_null must be")
  expect_error(do.call(power_cc, c(cfg, list(pi_null = 1.5))), "pi_null must be")
})

# Reference values in ref_power_cc.csv were produced by an independent
# implementation written in Python (SciPy), using the full double enumeration
# over completer counts and responder counts. They therefore check the R
# implementation against a different language and a different summation order,
# not merely against itself.

ref <- utils::read.csv(
  testthat::test_path("ref_power_cc.csv"),
  stringsAsFactors = FALSE
)

test_that("power_cc reproduces an independent Python enumeration", {
  for (i in seq_len(nrow(ref))) {
    r <- ref[i, ]
    out <- power_cc(
      n1 = r$n1, n2 = r$n2, p1 = r$p1, p2 = r$p2,
      omega1 = r$omega1, omega2 = r$omega2,
      rho1 = r$rho1, rho2 = r$rho2, alpha = r$alpha
    )
    info <- sprintf("row %d: n = (%d, %d), rho = (%g, %g), omega = (%g, %g)",
                    i, r$n1, r$n2, r$rho1, r$rho2, r$omega1, r$omega2)
    expect_equal(out$p1_CC, r$p1_CC, tolerance = 1e-12, info = info)
    expect_equal(out$p2_CC, r$p2_CC, tolerance = 1e-12, info = info)
    expect_equal(out$power, r$power, tolerance = 1e-10, info = info)
    expect_equal(out$type1_error, r$type1_error, tolerance = 1e-10, info = info)
  }
})

test_that("power_cc agrees with a direct double enumeration", {
  # Deliberately naive reference: enumerate every (m1, m2, x1, x2) cell
  naive_power_cc <- function(n1, n2, p1, p2, omega1, omega2,
                             rho1 = 0, rho2 = 0, alpha = 0.025) {
    s1 <- sqrt(p1 * (1 - p1) * omega1 * (1 - omega1))
    s2 <- sqrt(p2 * (1 - p2) * omega2 * (1 - omega2))
    pi1_00 <- (1 - p1) * (1 - omega1) + rho1 * s1
    pi1_10 <- p1 * (1 - omega1) - rho1 * s1
    pi2_00 <- (1 - p2) * (1 - omega2) + rho2 * s2
    pi2_10 <- p2 * (1 - omega2) - rho2 * s2
    p1_CC <- pi1_10 / (pi1_00 + pi1_10)
    p2_CC <- pi2_10 / (pi2_00 + pi2_10)
    p0_CC <- (n1 * p1_CC + n2 * p2_CC) / (n1 + n2)
    d0 <- (n1 * omega1 + n2 * omega2) / (n1 + n2)
    z <- stats::qnorm(1 - alpha)
    pm1 <- stats::dbinom(0:n1, n1, 1 - omega1)
    pm2 <- stats::dbinom(0:n2, n2, 1 - omega2)
    qm1 <- stats::dbinom(0:n1, n1, 1 - d0)
    qm2 <- stats::dbinom(0:n2, n2, 1 - d0)
    pw <- 0
    t1 <- 0
    for (m1 in 1:n1) {
      for (m2 in 1:n2) {
        x1 <- 0:m1
        x2 <- 0:m2
        pp <- outer(x1, x2, function(a, b) (a + b) / (m1 + m2))
        se <- sqrt(pp * (1 - pp) * (1 / m1 + 1 / m2))
        se[se < 1e-10] <- 1e-10
        rej <- (outer(x1 / m1, rep(1, m2 + 1)) -
                  outer(rep(1, m1 + 1), x2 / m2)) / se > z
        pw <- pw + pm1[m1 + 1] * pm2[m2 + 1] *
          sum(outer(stats::dbinom(x1, m1, p1_CC),
                    stats::dbinom(x2, m2, p2_CC)) * rej)
        t1 <- t1 + qm1[m1 + 1] * qm2[m2 + 1] *
          sum(outer(stats::dbinom(x1, m1, p0_CC),
                    stats::dbinom(x2, m2, p0_CC)) * rej)
      }
    }
    c(power = pw, type1_error = t1)
  }

  cfg <- list(n1 = 22, n2 = 30, p1 = 0.55, p2 = 0.30,
              omega1 = 0.18, omega2 = 0.08, rho1 = 0.15, rho2 = -0.10,
              alpha = 0.025)
  fast <- do.call(power_cc, cfg)
  slow <- do.call(naive_power_cc, cfg)

  expect_equal(fast$power, unname(slow["power"]), tolerance = 1e-12)
  expect_equal(fast$type1_error, unname(slow["type1_error"]), tolerance = 1e-12)
})

test_that("tol = 0 is exact and a positive tol is a lower bound within eps", {
  cfg <- list(n1 = 80, n2 = 80, p1 = 0.45, p2 = 0.25,
              omega1 = 0.15, omega2 = 0.15, rho1 = 0, rho2 = 0)
  exact <- do.call(power_cc, cfg)
  eps <- 1e-10
  approx <- do.call(power_cc, c(cfg, list(tol = eps)))
  n_cells <- (cfg$n1 + 1) * (cfg$n2 + 1)

  expect_lte(approx$power, exact$power + 1e-15)
  expect_gte(approx$power, exact$power - n_cells * eps)
  expect_lte(approx$type1_error, exact$type1_error + 1e-15)
  expect_gte(approx$type1_error, exact$type1_error - n_cells * eps)
})

test_that("rho = 0 gives p_CC equal to the latent response probability", {
  out <- power_cc(n1 = 40, n2 = 40, p1 = 0.6, p2 = 0.4,
                  omega1 = 0.25, omega2 = 0.10)
  expect_equal(out$p1_CC, 0.6, tolerance = 1e-12)
  expect_equal(out$p2_CC, 0.4, tolerance = 1e-12)
  expect_equal(out$effect_CC, 0.2, tolerance = 1e-12)
})

test_that("power_cc rejects invalid input", {
  expect_error(power_cc(0, 40, 0.6, 0.4, 0.1, 0.1), "positive")
  expect_error(power_cc(40, 40, 1.2, 0.4, 0.1, 0.1), "p1 must be")
  expect_error(power_cc(40, 40, 0.6, 0.4, 0.1, 0.1, tol = -1), "non-negative")
  expect_error(power_cc(40, 40, 0.6, 0.4, 0.1, 0.1, rho1 = 5), "rho1 must be")
})

test_that("p_null defaults to the pooled value and can be set explicitly", {
  cfg <- list(n1 = 40, n2 = 40, p1 = 0.6, p2 = 0.4,
              omega1 = 0.15, omega2 = 0.15)
  default <- do.call(power_cc, cfg)
  expect_equal(default$p_null, 0.5, tolerance = 1e-12)

  explicit <- do.call(power_cc, c(cfg, list(p_null = 0.5)))
  expect_equal(explicit$type1_error, default$type1_error, tolerance = 1e-14)
  expect_equal(explicit$power, default$power, tolerance = 1e-14)
})

test_that("p_null changes the size but never the power", {
  # Reference values from the independent Python enumeration, transcribed to
  # ten decimal places, so the comparison is absolute at that precision rather
  # than relative; see the same test in test-power_nri.R for the reasoning.
  cfg <- list(n1 = 40, n2 = 40, p1 = 0.6, p2 = 0.4,
              omega1 = 0.15, omega2 = 0.15)
  expected <- c("0.5" = 0.0256208522,
                "0.4" = 0.0257491919,
                "0.25" = 0.0258309009)
  base_power <- do.call(power_cc, cfg)$power
  for (nm in names(expected)) {
    out <- do.call(power_cc, c(cfg, list(p_null = as.numeric(nm))))
    expect_lt(abs(out$type1_error - unname(expected[nm])), 1e-9,
              label = paste("|type1_error - reference| at p_null =", nm))
    expect_equal(out$power, base_power, tolerance = 1e-14, info = nm)
  }
})

test_that("power_cc rejects an invalid p_null", {
  cfg <- list(n1 = 30, n2 = 30, p1 = 0.6, p2 = 0.4,
              omega1 = 0.1, omega2 = 0.1)
  expect_error(do.call(power_cc, c(cfg, list(p_null = 0))), "p_null must be")
  expect_error(do.call(power_cc, c(cfg, list(p_null = 1.5))), "p_null must be")
})

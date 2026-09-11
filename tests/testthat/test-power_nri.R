# Reference values in ref_power_nri.csv were produced by an independent
# implementation written in Python (SciPy) that enumerates the full outcome
# matrix, so the check crosses both language and summation order.

ref_nri <- utils::read.csv(
  testthat::test_path("ref_power_nri.csv"),
  stringsAsFactors = FALSE
)

test_that("power_nri reproduces an independent Python enumeration", {
  for (i in seq_len(nrow(ref_nri))) {
    r <- ref_nri[i, ]
    out <- power_nri(
      n1 = r$n1, n0 = r$n0, pi1 = r$pi1, pi0 = r$pi0,
      omega1 = r$omega1, omega0 = r$omega0,
      rho1 = r$rho1, rho0 = r$rho0, alpha = r$alpha
    )
    info <- sprintf("row %d: n = (%d, %d), rho = (%g, %g), omega = (%g, %g)",
                    i, r$n1, r$n0, r$rho1, r$rho0, r$omega1, r$omega0)
    expect_equal(out$pi1_NRI, r$pi1_NRI, tolerance = 1e-12, info = info)
    expect_equal(out$pi0_NRI, r$pi0_NRI, tolerance = 1e-12, info = info)
    expect_equal(out$power, r$power, tolerance = 1e-10, info = info)
    expect_equal(out$type1_error, r$type1_error, tolerance = 1e-10, info = info)
  }
})

test_that("power_nri agrees with a direct matrix enumeration", {
  naive_power_nri <- function(n1, n0, pi1, pi0, omega1, omega0,
                              rho1 = 0, rho0 = 0, alpha = 0.025) {
    q1 <- pi1 * (1 - omega1) -
      rho1 * sqrt(pi1 * (1 - pi1) * omega1 * (1 - omega1))
    q0 <- pi0 * (1 - omega0) -
      rho0 * sqrt(pi0 * (1 - pi0) * omega0 * (1 - omega0))
    pi_bar <- (n1 * q1 + n0 * q0) / (n1 + n0)
    z <- stats::qnorm(1 - alpha)
    x1 <- 0:n1
    x0 <- 0:n0
    pp <- outer(x1, x0, function(a, b) (a + b) / (n1 + n0))
    se <- sqrt(pp * (1 - pp) * (1 / n1 + 1 / n0))
    se[se < 1e-10] <- 1e-10
    rej <- (outer(x1 / n1, rep(1, n0 + 1)) -
              outer(rep(1, n1 + 1), x0 / n0)) / se > z
    c(power = sum(outer(stats::dbinom(x1, n1, q1),
                        stats::dbinom(x0, n0, q0))[rej]),
      type1_error = sum(outer(stats::dbinom(x1, n1, pi_bar),
                              stats::dbinom(x0, n0, pi_bar))[rej]))
  }

  cfg <- list(n1 = 63, n0 = 47, pi1 = 0.52, pi0 = 0.31,
              omega1 = 0.22, omega0 = 0.09, rho1 = -0.12, rho0 = 0.18,
              alpha = 0.025)
  fast <- do.call(power_nri, cfg)
  slow <- do.call(naive_power_nri, cfg)

  expect_equal(fast$power, unname(slow["power"]), tolerance = 1e-12)
  expect_equal(fast$type1_error, unname(slow["type1_error"]), tolerance = 1e-12)
})

test_that("the NRI estimand follows the closed-form expression", {
  # delta_NRI = pi1(1-w1) - pi0(1-w0) - rho1 s1 + rho0 s0
  grid <- expand.grid(pi1 = c(0.35, 0.60), pi0 = c(0.20, 0.40),
                      omega1 = c(0.05, 0.25), omega0 = c(0.05, 0.25),
                      rho1 = c(-0.1, 0, 0.1), rho0 = c(-0.1, 0, 0.1))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    s1 <- sqrt(g$pi1 * (1 - g$pi1) * g$omega1 * (1 - g$omega1))
    s0 <- sqrt(g$pi0 * (1 - g$pi0) * g$omega0 * (1 - g$omega0))
    expected <- g$pi1 * (1 - g$omega1) - g$pi0 * (1 - g$omega0) -
      g$rho1 * s1 + g$rho0 * s0
    out <- power_nri(n1 = 20, n0 = 20, pi1 = g$pi1, pi0 = g$pi0,
                     omega1 = g$omega1, omega0 = g$omega0,
                     rho1 = g$rho1, rho0 = g$rho0)
    expect_equal(out$effect_NRI, expected, tolerance = 1e-12)
  }
})

test_that("equal dropout and independence give delta_NRI = (1 - omega) delta_CC", {
  for (omega in c(0.05, 0.10, 0.25)) {
    out <- power_nri(n1 = 30, n0 = 30, pi1 = 0.55, pi0 = 0.35,
                     omega1 = omega, omega0 = omega)
    expect_equal(out$effect_NRI, (1 - omega) * 0.20, tolerance = 1e-12)
  }
})

test_that("unequal dropout makes the NRI estimand non-null under equal latent pi", {
  # This is the setting Reviewer 2 asked about: pi1 = pi0 but omega1 != omega0.
  # The NRI contrast is not zero, so the test is sized against a different
  # null hypothesis rather than exhibiting an inflated type I error rate.
  out <- power_nri(n1 = 60, n0 = 60, pi1 = 0.40, pi0 = 0.40,
                   omega1 = 0.05, omega0 = 0.20)
  expect_equal(out$effect_full, 0, tolerance = 1e-12)
  expect_equal(out$effect_NRI, 0.40 * 0.95 - 0.40 * 0.80, tolerance = 1e-12)
  expect_gt(out$effect_NRI, 0)
  # The reported type I error rate is computed at the pooled NRI probability,
  # so it stays close to the nominal level
  expect_lt(abs(out$type1_error - 0.025), 0.01)
})

test_that("power_nri rejects invalid input", {
  expect_error(power_nri(0, 40, 0.6, 0.4, 0.1, 0.1), "positive")
  expect_error(power_nri(40, 40, 0.6, 1.4, 0.1, 0.1), "pi0 must be")
  expect_error(power_nri(40, 40, 0.6, 0.4, 0.1, 0.1, rho0 = -5), "rho0 must be")
})

test_that("pi_null defaults to the pooled value and can be set explicitly", {
  cfg <- list(n1 = 80, n0 = 80, pi1 = 0.45, pi0 = 0.30,
              omega1 = 0.10, omega0 = 0.10)
  default <- do.call(power_nri, cfg)
  pooled <- (0.45 * 0.90 + 0.30 * 0.90) / 2

  expect_equal(default$pi_null, pooled, tolerance = 1e-12)

  explicit <- do.call(power_nri, c(cfg, list(pi_null = pooled)))
  expect_equal(explicit$type1_error, default$type1_error, tolerance = 1e-14)
  expect_equal(explicit$power, default$power, tolerance = 1e-14)
})

test_that("pi_null changes the size but never the power", {
  # Reference values from the independent Python enumeration, transcribed to
  # ten decimal places. They therefore carry a rounding error of up to 5e-11,
  # which at values near 0.025 is about 2e-9 in relative terms. expect_equal()
  # compares numerics relatively, so the comparison below is made absolute at
  # the precision the constants actually have. An absolute 1e-9 is still about
  # 4e-8 relative here, so nothing is being waved through.
  cfg <- list(n1 = 80, n0 = 80, pi1 = 0.45, pi0 = 0.30,
              omega1 = 0.10, omega0 = 0.10)
  expected <- c(
    "0.3375" = 0.0255719673,   # pooled
    "0.27"   = 0.0258363803,   # pi0 on the NRI scale
    "0.20"   = 0.0251387765,
    "0.50"   = 0.0238911181,
    "0.70"   = 0.0256583857
  )
  base_power <- do.call(power_nri, cfg)$power
  for (nm in names(expected)) {
    out <- do.call(power_nri, c(cfg, list(pi_null = as.numeric(nm))))
    expect_lt(abs(out$type1_error - unname(expected[nm])), 1e-9,
              label = paste("|type1_error - reference| at pi_null =", nm))
    expect_equal(out$pi_null, as.numeric(nm), tolerance = 1e-12)
    # The alternative is untouched by the null nuisance parameter
    expect_equal(out$power, base_power, tolerance = 1e-14, info = nm)
  }
  # The point of the argument: the size is not constant in the nuisance
  # parameter, so a single value does not establish size control
  sizes <- unname(expected)
  expect_gt(max(sizes) - min(sizes), 1e-3)
})

test_that("power_nri rejects an invalid pi_null", {
  cfg <- list(n1 = 40, n0 = 40, pi1 = 0.6, pi0 = 0.4,
              omega1 = 0.1, omega0 = 0.1)
  expect_error(do.call(power_nri, c(cfg, list(pi_null = 0))), "pi_null must be")
  expect_error(do.call(power_nri, c(cfg, list(pi_null = 1))), "pi_null must be")
  expect_error(do.call(power_nri, c(cfg, list(pi_null = c(0.3, 0.4)))),
               "pi_null must be")
  expect_error(do.call(power_nri, c(cfg, list(pi_null = NA_real_))),
               "pi_null must be")
})

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
      n1 = r$n1, n2 = r$n2, p1 = r$p1, p2 = r$p2,
      omega1 = r$omega1, omega2 = r$omega2,
      rho1 = r$rho1, rho2 = r$rho2, alpha = r$alpha
    )
    info <- sprintf("row %d: n = (%d, %d), rho = (%g, %g), omega = (%g, %g)",
                    i, r$n1, r$n2, r$rho1, r$rho2, r$omega1, r$omega2)
    expect_equal(out$p1_NRI, r$p1_NRI, tolerance = 1e-12, info = info)
    expect_equal(out$p2_NRI, r$p2_NRI, tolerance = 1e-12, info = info)
    expect_equal(out$power, r$power, tolerance = 1e-10, info = info)
    expect_equal(out$type1_error, r$type1_error, tolerance = 1e-10, info = info)
  }
})

test_that("power_nri agrees with a direct matrix enumeration", {
  naive_power_nri <- function(n1, n2, p1, p2, omega1, omega2,
                              rho1 = 0, rho2 = 0, alpha = 0.025) {
    q1 <- p1 * (1 - omega1) -
      rho1 * sqrt(p1 * (1 - p1) * omega1 * (1 - omega1))
    q2 <- p2 * (1 - omega2) -
      rho2 * sqrt(p2 * (1 - p2) * omega2 * (1 - omega2))
    p0 <- (n1 * q1 + n2 * q2) / (n1 + n2)
    z <- stats::qnorm(1 - alpha)
    x1 <- 0:n1
    x2 <- 0:n2
    pp <- outer(x1, x2, function(a, b) (a + b) / (n1 + n2))
    se <- sqrt(pp * (1 - pp) * (1 / n1 + 1 / n2))
    se[se < 1e-10] <- 1e-10
    rej <- (outer(x1 / n1, rep(1, n2 + 1)) -
              outer(rep(1, n1 + 1), x2 / n2)) / se > z
    c(power = sum(outer(stats::dbinom(x1, n1, q1),
                        stats::dbinom(x2, n2, q2))[rej]),
      type1_error = sum(outer(stats::dbinom(x1, n1, p0),
                              stats::dbinom(x2, n2, p0))[rej]))
  }

  cfg <- list(n1 = 63, n2 = 47, p1 = 0.52, p2 = 0.31,
              omega1 = 0.22, omega2 = 0.09, rho1 = -0.12, rho2 = 0.18,
              alpha = 0.025)
  fast <- do.call(power_nri, cfg)
  slow <- do.call(naive_power_nri, cfg)

  expect_equal(fast$power, unname(slow["power"]), tolerance = 1e-12)
  expect_equal(fast$type1_error, unname(slow["type1_error"]), tolerance = 1e-12)
})

test_that("the NRI estimand follows the closed-form expression", {
  # delta_NRI = p1(1-w1) - p2(1-w2) - rho1 s1 + rho2 s2
  grid <- expand.grid(p1 = c(0.35, 0.60), p2 = c(0.20, 0.40),
                      omega1 = c(0.05, 0.25), omega2 = c(0.05, 0.25),
                      rho1 = c(-0.1, 0, 0.1), rho2 = c(-0.1, 0, 0.1))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    s1 <- sqrt(g$p1 * (1 - g$p1) * g$omega1 * (1 - g$omega1))
    s2 <- sqrt(g$p2 * (1 - g$p2) * g$omega2 * (1 - g$omega2))
    expected <- g$p1 * (1 - g$omega1) - g$p2 * (1 - g$omega2) -
      g$rho1 * s1 + g$rho2 * s2
    out <- power_nri(n1 = 20, n2 = 20, p1 = g$p1, p2 = g$p2,
                     omega1 = g$omega1, omega2 = g$omega2,
                     rho1 = g$rho1, rho2 = g$rho2)
    expect_equal(out$effect_NRI, expected, tolerance = 1e-12)
  }
})

test_that("equal dropout and independence give delta_NRI = (1 - omega) delta_CC", {
  for (omega in c(0.05, 0.10, 0.25)) {
    out <- power_nri(n1 = 30, n2 = 30, p1 = 0.55, p2 = 0.35,
                     omega1 = omega, omega2 = omega)
    expect_equal(out$effect_NRI, (1 - omega) * 0.20, tolerance = 1e-12)
  }
})

test_that("unequal dropout makes the NRI estimand non-null under equal latent p", {
  # This is the setting Reviewer 2 asked about: p1 = p2 but omega1 != omega2.
  # The NRI contrast is not zero, so the test is sized against a different
  # null hypothesis rather than exhibiting an inflated type I error rate.
  out <- power_nri(n1 = 60, n2 = 60, p1 = 0.40, p2 = 0.40,
                   omega1 = 0.05, omega2 = 0.20)
  expect_equal(out$effect_latent, 0, tolerance = 1e-12)
  expect_equal(out$effect_NRI, 0.40 * 0.95 - 0.40 * 0.80, tolerance = 1e-12)
  expect_gt(out$effect_NRI, 0)
  # The reported type I error rate is computed at the pooled NRI probability,
  # so it stays close to the nominal level
  expect_lt(abs(out$type1_error - 0.025), 0.01)
})

test_that("power_nri rejects invalid input", {
  expect_error(power_nri(0, 40, 0.6, 0.4, 0.1, 0.1), "positive")
  expect_error(power_nri(40, 40, 0.6, 1.4, 0.1, 0.1), "p2 must be")
  expect_error(power_nri(40, 40, 0.6, 0.4, 0.1, 0.1, rho2 = -5), "rho2 must be")
})

test_that("p_null defaults to the pooled value and can be set explicitly", {
  cfg <- list(n1 = 80, n2 = 80, p1 = 0.45, p2 = 0.30,
              omega1 = 0.10, omega2 = 0.10)
  default <- do.call(power_nri, cfg)
  pooled <- (0.45 * 0.90 + 0.30 * 0.90) / 2

  expect_equal(default$p_null, pooled, tolerance = 1e-12)

  explicit <- do.call(power_nri, c(cfg, list(p_null = pooled)))
  expect_equal(explicit$type1_error, default$type1_error, tolerance = 1e-14)
  expect_equal(explicit$power, default$power, tolerance = 1e-14)
})

test_that("p_null changes the size but never the power", {
  # Reference values from the independent Python enumeration, transcribed to
  # ten decimal places. They therefore carry a rounding error of up to 5e-11,
  # which at values near 0.025 is about 2e-9 in relative terms. expect_equal()
  # compares numerics relatively, so the comparison below is made absolute at
  # the precision the constants actually have. An absolute 1e-9 is still about
  # 4e-8 relative here, so nothing is being waved through.
  cfg <- list(n1 = 80, n2 = 80, p1 = 0.45, p2 = 0.30,
              omega1 = 0.10, omega2 = 0.10)
  expected <- c(
    "0.3375" = 0.0255719673,   # pooled
    "0.27"   = 0.0258363803,   # p2 on the NRI scale
    "0.20"   = 0.0251387765,
    "0.50"   = 0.0238911181,
    "0.70"   = 0.0256583857
  )
  base_power <- do.call(power_nri, cfg)$power
  for (nm in names(expected)) {
    out <- do.call(power_nri, c(cfg, list(p_null = as.numeric(nm))))
    expect_lt(abs(out$type1_error - unname(expected[nm])), 1e-9,
              label = paste("|type1_error - reference| at p_null =", nm))
    expect_equal(out$p_null, as.numeric(nm), tolerance = 1e-12)
    # The alternative is untouched by the null nuisance parameter
    expect_equal(out$power, base_power, tolerance = 1e-14, info = nm)
  }
  # The point of the argument: the size is not constant in the nuisance
  # parameter, so a single value does not establish size control
  sizes <- unname(expected)
  expect_gt(max(sizes) - min(sizes), 1e-3)
})

test_that("power_nri rejects an invalid p_null", {
  cfg <- list(n1 = 40, n2 = 40, p1 = 0.6, p2 = 0.4,
              omega1 = 0.1, omega2 = 0.1)
  expect_error(do.call(power_nri, c(cfg, list(p_null = 0))), "p_null must be")
  expect_error(do.call(power_nri, c(cfg, list(p_null = 1))), "p_null must be")
  expect_error(do.call(power_nri, c(cfg, list(p_null = c(0.3, 0.4)))),
               "p_null must be")
  expect_error(do.call(power_nri, c(cfg, list(p_null = NA_real_))),
               "p_null must be")
})

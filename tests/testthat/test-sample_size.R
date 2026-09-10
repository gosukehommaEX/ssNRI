test_that("sample_size_cc reproduces the textbook formula without dropout", {
  # Two-group comparison of proportions, one-sided, pooled variance under H0
  manual_n2 <- function(p1, p2, r, alpha, power) {
    rp <- (r * p1 + p2) / (1 + r)
    v0 <- rp * (1 - rp)
    v1 <- (p1 * (1 - p1) / r + p2 * (1 - p2)) / (1 + 1 / r)
    ceiling((1 + 1 / r) / ((p1 - p2) ^ 2) *
              (stats::qnorm(1 - alpha) * sqrt(v0) +
                 stats::qnorm(power) * sqrt(v1)) ^ 2)
  }
  for (r in c(0.5, 1, 2)) {
    ss <- sample_size_cc(p1 = 0.6, p2 = 0.4, r = r)
    expect_equal(ss$n2, manual_n2(0.6, 0.4, r, 0.025, 0.8))
    expect_equal(ss$n1, ceiling(r * ss$n2))
    expect_equal(ss$n_total, ss$n1 + ss$n2)
    # No dropout, so no inflation
    expect_equal(ss$n_total_adjust, ss$n_total)
  }
})

test_that("the simple inflation step is exactly n / (1 - omega)", {
  base <- sample_size_cc(p1 = 0.6, p2 = 0.4)
  infl <- sample_size_cc(p1 = 0.6, p2 = 0.4, omega1 = 0.2, omega2 = 0.3)
  expect_equal(infl$n1, base$n1)
  expect_equal(infl$n2, base$n2)
  expect_equal(infl$n1_adjust, ceiling(base$n1 / 0.8))
  expect_equal(infl$n2_adjust, ceiling(base$n2 / 0.7))
})

test_that("sample_size_nri returns a size that attains the target power", {
  # The closed form uses a normal approximation, so exact power lands near the
  # target rather than on it; check that it is close and not below by much
  settings <- list(
    list(p1 = 0.60, p2 = 0.40, omega1 = 0.20, omega2 = 0.20, rho1 = 0, rho2 = 0),
    list(p1 = 0.35, p2 = 0.20, omega1 = 0.10, omega2 = 0.10, rho1 = 0, rho2 = 0),
    list(p1 = 0.55, p2 = 0.35, omega1 = 0.15, omega2 = 0.05, rho1 = 0.1, rho2 = -0.1)
  )
  for (s in settings) {
    ss <- do.call(sample_size_nri, c(s, list(target_power = 0.8)))
    pw <- power_nri(n1 = ss$n1, n2 = ss$n2, p1 = s$p1, p2 = s$p2,
                    omega1 = s$omega1, omega2 = s$omega2,
                    rho1 = s$rho1, rho2 = s$rho2)$power
    expect_gt(pw, 0.78)
    expect_lt(pw, 0.83)
  }
})

test_that("the simple inflation method falls short under equal dropout", {
  # delta_NRI = (1 - omega) delta_CC < delta_CC, so a larger sample is needed
  for (omega in c(0.05, 0.10, 0.20, 0.30)) {
    cc <- sample_size_cc(p1 = 0.6, p2 = 0.4, omega1 = omega, omega2 = omega)
    nri <- sample_size_nri(p1 = 0.6, p2 = 0.4, omega1 = omega, omega2 = omega)
    expect_lt(cc$n_total_adjust, nri$n_total)
  }
})

test_that("sample_size_nri with no dropout matches the dropout-free design", {
  cc <- sample_size_cc(p1 = 0.6, p2 = 0.4)
  nri <- sample_size_nri(p1 = 0.6, p2 = 0.4)
  expect_equal(nri$n_total, cc$n_total)
  expect_equal(nri$effect_NRI, nri$effect_latent, tolerance = 1e-12)
})

test_that("sample size functions reject invalid input", {
  expect_error(sample_size_cc(p1 = 1.2, p2 = 0.4), "p1 must be")
  expect_error(sample_size_cc(p1 = 0.6, p2 = 0.4, r = 0), "r must be positive")
  expect_error(sample_size_nri(p1 = 0.6, p2 = 0.4, omega1 = 0.2, rho1 = 9),
               "rho1 must be")
  # p1 <= p2 warns but still computes, because the formula uses the squared effect
  expect_warning(sample_size_nri(p1 = 0.40, p2 = 0.45), "superiority")
  expect_warning(sample_size_cc(p1 = 0.40, p2 = 0.45), "superiority")
})

# Direct tests of the internal helpers. Both of the bugs found so far lived
# here and were only detected through the exported functions, where the symptom
# was far from the cause, so the helpers are now pinned on their own.

test_that(".joint_cells keeps its names, including after clamping", {
  expected <- c("pi_10", "pi_11", "pi_00", "pi_01")

  # Interior of the feasible region
  cells <- ssNRI:::.joint_cells(pi = 0.4, omega = 0.15, rho = 0.1)
  expect_named(cells, expected)
  expect_false(anyNA(cells))

  # At the feasible boundary, where the clamp actually fires
  for (pi in c(0.2, 0.5, 0.8)) {
    for (omega in c(0.05, 0.2, 0.45)) {
      b <- rho_bounds(pi = pi, omega = omega)
      for (rho in c(b$rho_lower, 0, b$rho_upper)) {
        cells <- ssNRI:::.joint_cells(pi = pi, omega = omega, rho = rho)
        expect_named(cells, expected,
                     info = sprintf("pi = %g, omega = %g, rho = %g",
                                    pi, omega, rho))
        expect_true(all(cells >= 0))
        expect_true(all(cells <= 1))
      }
    }
  }
})

test_that(".joint_cells reproduces the marginals and the correlation", {
  grid <- expand.grid(pi = c(0.25, 0.50, 0.75), omega = c(0.05, 0.20, 0.40),
                      rho = c(-0.15, 0, 0.15))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    b <- rho_bounds(pi = g$pi, omega = g$omega)
    if (g$rho < b$rho_lower || g$rho > b$rho_upper) next
    cells <- ssNRI:::.joint_cells(pi = g$pi, omega = g$omega, rho = g$rho)
    info <- sprintf("pi = %g, omega = %g, rho = %g", g$pi, g$omega, g$rho)
    expect_equal(sum(cells), 1, tolerance = 1e-12, info = info)
    expect_equal(cells[["pi_10"]] + cells[["pi_11"]], g$pi,
                 tolerance = 1e-12, info = info)
    expect_equal(cells[["pi_11"]] + cells[["pi_01"]], g$omega,
                 tolerance = 1e-12, info = info)
    back <- (cells[["pi_11"]] - g$pi * g$omega) /
      sqrt(g$pi * (1 - g$pi) * g$omega * (1 - g$omega))
    expect_equal(back, g$rho, tolerance = 1e-10, info = info)
  }
})

# Deliberately naive threshold: scan x1 upwards until the test rejects
brute_threshold <- function(m1, m0, z_alpha) {
  vapply(0:m0, function(j) {
    for (u in 0:m1) {
      pp <- (u + j) / (m1 + m0)
      se <- max(sqrt(pp * (1 - pp) * (1 / m1 + 1 / m0)), 1e-10)
      if ((u / m1 - j / m0) / se > z_alpha) return(u)
    }
    m1 + 1
  }, numeric(1))
}

test_that(".reject_threshold matches a direct scan", {
  z <- stats::qnorm(0.975)
  for (m1 in c(1, 2, 5, 17, 40, 63)) {
    for (m0 in c(1, 3, 8, 23, 55)) {
      fast <- ssNRI:::.reject_threshold(m1, m0, z)
      expect_equal(length(fast), m0 + 1L,
                   info = sprintf("m1 = %d, m0 = %d", m1, m0))
      expect_equal(fast, brute_threshold(m1, m0, z),
                   info = sprintf("m1 = %d, m0 = %d", m1, m0))
    }
  }
})

test_that(".reject_threshold survives group sizes past the integer limit", {
  # The quadratic coefficients involve a product of five counts. In integer
  # arithmetic that overflows from m1 = m0 = 64 upwards, which used to produce
  # NA and an unrelated error further downstream.
  z <- stats::qnorm(0.975)
  for (m in c(64, 65, 90, 150)) {
    fast <- ssNRI:::.reject_threshold(m, m, z)
    expect_false(anyNA(fast), info = sprintf("m = %d", m))
    expect_equal(fast, brute_threshold(m, m, z), info = sprintf("m = %d", m))
  }
  # Unequal sizes, one of them large
  expect_equal(ssNRI:::.reject_threshold(200, 30, z),
               brute_threshold(200, 30, z))
  expect_equal(ssNRI:::.reject_threshold(30, 200, z),
               brute_threshold(30, 200, z))
})

test_that(".reject_threshold is non-decreasing in the group 0 count", {
  z <- stats::qnorm(0.975)
  for (m1 in c(7, 31, 80)) {
    for (m0 in c(5, 29, 77)) {
      cth <- ssNRI:::.reject_threshold(m1, m0, z)
      expect_true(all(diff(cth) >= 0),
                  info = sprintf("m1 = %d, m0 = %d", m1, m0))
      expect_true(all(cth >= 0 & cth <= m1 + 1),
                  info = sprintf("m1 = %d, m0 = %d", m1, m0))
    }
  }
})

test_that(".reject_threshold responds to alpha in the right direction", {
  cth_strict <- ssNRI:::.reject_threshold(50, 50, stats::qnorm(0.995))
  cth_loose <- ssNRI:::.reject_threshold(50, 50, stats::qnorm(0.95))
  expect_true(all(cth_strict >= cth_loose))
})

test_that(".common_rho_range intersects the two feasible ranges", {
  rng <- ssNRI:::.common_rho_range(pi1 = 0.6, pi0 = 0.4,
                                   omega1 = 0.2, omega0 = 0.2)
  b1 <- rho_bounds(pi = 0.6, omega = 0.2)
  b0 <- rho_bounds(pi = 0.4, omega = 0.2)
  expect_equal(rng[1], max(b1$rho_lower, b0$rho_lower))
  expect_equal(rng[2], min(b1$rho_upper, b0$rho_upper))
})

test_that(".clamp_rho moves a value into the feasible range and leaves others", {
  b <- rho_bounds(pi = 0.3, omega = 0.1)
  expect_equal(ssNRI:::.clamp_rho(0, 0.3, 0.1), 0)
  expect_equal(ssNRI:::.clamp_rho(99, 0.3, 0.1), b$rho_upper)
  expect_equal(ssNRI:::.clamp_rho(-99, 0.3, 0.1), b$rho_lower)
})

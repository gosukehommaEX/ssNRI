test_that("rho_bounds returns the Prentice bounds", {
  b <- rho_bounds(p = 0.35, omega = 0.10)
  expect_equal(b$rho_lower,
               max(-sqrt(0.35 * 0.10 / (0.65 * 0.90)),
                   -sqrt(0.65 * 0.90 / (0.35 * 0.10))),
               tolerance = 1e-12)
  expect_equal(b$rho_upper,
               min(sqrt(0.35 * 0.90 / (0.10 * 0.65)),
                   sqrt(0.10 * 0.65 / (0.35 * 0.90))),
               tolerance = 1e-12)
})

test_that("the bounds are attainable and nothing outside them is", {
  # At rho_upper the joint cell probability pi^(1,1) or pi^(0,0) hits a boundary
  for (p in c(0.2, 0.5, 0.8)) {
    for (omega in c(0.05, 0.2, 0.4)) {
      b <- rho_bounds(p = p, omega = omega)
      s <- sqrt(p * (1 - p) * omega * (1 - omega))
      for (rho in c(b$rho_lower, 0, b$rho_upper)) {
        cells <- c(
          p * (1 - omega) - rho * s,
          p * omega + rho * s,
          (1 - p) * (1 - omega) + rho * s,
          (1 - p) * omega - rho * s
        )
        expect_true(all(cells >= -1e-12))
        expect_equal(sum(cells), 1, tolerance = 1e-12)
      }
    }
  }
})

test_that("omega = 0 gives a degenerate range", {
  b <- rho_bounds(p = 0.4, omega = 0)
  expect_equal(b$rho_lower, 0)
  expect_equal(b$rho_upper, 0)
})

test_that("rho_bounds rejects invalid input", {
  expect_error(rho_bounds(p = 0, omega = 0.1), "strictly between")
  expect_error(rho_bounds(p = 0.4, omega = 1), "between 0")
})

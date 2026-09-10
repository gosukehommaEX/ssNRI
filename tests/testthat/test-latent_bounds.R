test_that("the identified range has width exactly omega under NRI", {
  for (nd in c(10, 30, 60)) {
    lb <- latent_bounds(p_obs = 0.40, N_randomized = 200, N_dropout = nd)
    expect_equal(lb$omega, nd / 200)
    expect_equal(lb$p_lower, 0.40)
    expect_equal(lb$p_upper, 0.40 + nd / 200)
    expect_equal(lb$p_upper - lb$p_lower, lb$omega, tolerance = 1e-12)
    expect_equal(lb$p_midpoint, (lb$p_lower + lb$p_upper) / 2)
  }
})

test_that("a complete case rate is converted to the NRI scale first", {
  lb <- latent_bounds(p_obs = 0.50, N_randomized = 100, N_dropout = 20,
                      method = "cc")
  expect_equal(lb$p_lower, 0.50 * 0.80, tolerance = 1e-12)
  expect_equal(lb$p_upper, 0.50 * 0.80 + 0.20, tolerance = 1e-12)
})

test_that("both endpoints of the range are attainable", {
  # At p_lower no dropout responded; at p_upper every dropout responded.
  # Both must give back the reported NRI rate.
  lb <- latent_bounds(p_obs = 0.40, N_randomized = 200, N_dropout = 30)
  omega <- lb$omega
  for (p in c(lb$p_lower, lb$p_upper)) {
    rho <- infer_rho(p_obs = 0.40, omega = omega, p = p, method = "nri")$rho
    p_nri <- p * (1 - omega) - rho * sqrt(p * (1 - p) * omega * (1 - omega))
    expect_equal(p_nri, 0.40, tolerance = 1e-10)
  }
})

test_that("no dropout identifies the latent probability", {
  expect_message(
    lb <- latent_bounds(p_obs = 0.40, N_randomized = 100, N_dropout = 0),
    "with certainty"
  )
  expect_equal(lb$p_lower, lb$p_upper)
  expect_equal(lb$p_lower, 0.40)
})

test_that("an incompatible reported rate is flagged", {
  expect_warning(
    latent_bounds(p_obs = 0.95, N_randomized = 100, N_dropout = 30),
    "exceeds 1"
  )
})

test_that("latent_bounds rejects invalid input", {
  expect_error(latent_bounds(p_obs = 1.5, N_randomized = 100, N_dropout = 10),
               "p_obs must be")
  expect_error(latent_bounds(p_obs = 0.4, N_randomized = 0, N_dropout = 0),
               "positive integer")
  expect_error(latent_bounds(p_obs = 0.4, N_randomized = 100, N_dropout = 200),
               "not exceeding")
})

test_that("the best and worst case rules shift the interval by omega", {
  # Counting dropouts as responders: reported = Pr(R = 1, D = 0) + omega
  bc_t <- latent_bounds(p_obs = 0.50, N_randomized = 100, N_dropout = 20,
                        method = "bc", group = "treatment")
  wc_c <- latent_bounds(p_obs = 0.50, N_randomized = 100, N_dropout = 20,
                        method = "wc", group = "control")
  expect_equal(bc_t$p_nri, 0.30, tolerance = 1e-12)
  expect_equal(bc_t$p_lower, 0.30, tolerance = 1e-12)
  expect_equal(bc_t$p_upper, 0.50, tolerance = 1e-12)
  expect_equal(wc_c$p_lower, bc_t$p_lower, tolerance = 1e-14)
  expect_equal(wc_c$p_upper, bc_t$p_upper, tolerance = 1e-14)

  # Counting dropouts as non-responders: identical to NRI
  nri  <- latent_bounds(p_obs = 0.50, N_randomized = 100, N_dropout = 20)
  bc_c <- latent_bounds(p_obs = 0.50, N_randomized = 100, N_dropout = 20,
                        method = "bc", group = "control")
  wc_t <- latent_bounds(p_obs = 0.50, N_randomized = 100, N_dropout = 20,
                        method = "wc", group = "treatment")
  expect_equal(bc_c$p_lower, nri$p_lower, tolerance = 1e-14)
  expect_equal(bc_c$p_upper, nri$p_upper, tolerance = 1e-14)
  expect_equal(wc_t$p_lower, nri$p_lower, tolerance = 1e-14)
  expect_equal(wc_t$p_upper, nri$p_upper, tolerance = 1e-14)

  # Every rule leaves an interval of width exactly omega
  for (lb in list(nri, bc_t, bc_c, wc_t, wc_c)) {
    expect_equal(lb$p_upper - lb$p_lower, lb$omega, tolerance = 1e-12)
  }
})

test_that("every endpoint reproduces the rate the historical trial reported", {
  omega <- 0.20
  for (m in c("nri", "cc", "bc", "wc")) {
    for (g in c("treatment", "control")) {
      lb <- latent_bounds(p_obs = 0.50, N_randomized = 100, N_dropout = 20,
                          method = m, group = g)
      for (pc in c(lb$p_lower, lb$p_midpoint, lb$p_upper)) {
        if (pc <= 0 || pc >= 1) next
        rho <- infer_rho(p_obs = 0.50, omega = omega, p = pc,
                         method = m, group = g)$rho
        p_nri <- pc * (1 - omega) -
          rho * sqrt(pc * (1 - pc) * omega * (1 - omega))
        expect_equal(p_nri, lb$p_nri, tolerance = 1e-10,
                     info = sprintf("%s / %s at p = %g", m, g, pc))
      }
    }
  }
})

test_that("a rate below the dropout probability is flagged", {
  expect_warning(
    latent_bounds(p_obs = 0.05, N_randomized = 100, N_dropout = 30,
                  method = "bc", group = "treatment"),
    "below 0"
  )
})

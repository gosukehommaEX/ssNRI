test_that("the implied correlation inverts the joint cell probability", {
  grid <- expand.grid(p = c(0.30, 0.50, 0.70), omega = c(0.05, 0.15, 0.30),
                      rho = c(-0.2, 0, 0.2))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    b <- rho_bounds(p = g$p, omega = g$omega)
    if (g$rho < b$rho_lower || g$rho > b$rho_upper) next
    p_nri <- g$p * (1 - g$omega) -
      g$rho * sqrt(g$p * (1 - g$p) * g$omega * (1 - g$omega))
    back <- infer_rho(p_obs = p_nri, omega = g$omega, p = g$p,
                      method = "nri")$rho
    expect_equal(back, g$rho, tolerance = 1e-10)
  }
})

test_that("a complete case rate gives the same correlation as its NRI form", {
  p <- 0.45
  omega <- 0.20
  rho <- 0.15
  p_nri <- p * (1 - omega) - rho * sqrt(p * (1 - p) * omega * (1 - omega))
  p_cc <- p_nri / (1 - omega)
  from_nri <- infer_rho(p_obs = p_nri, omega = omega, p = p, method = "nri")$rho
  from_cc <- infer_rho(p_obs = p_cc, omega = omega, p = p, method = "cc")$rho
  expect_equal(from_nri, from_cc, tolerance = 1e-12)
  expect_equal(from_cc, rho, tolerance = 1e-10)
})

test_that("a best or worst case rate is a joint cell, not an identified p", {
  # A best case treatment rate counts every dropout as a responder, so the
  # reported value is Pr(R = 1, D = 0) + omega. The value that setting the
  # remaining joint cell to zero used to give is the upper endpoint of the
  # interval, recovered by passing p = p_obs.
  res <- infer_rho(p_obs = 0.50, omega = 0.20, p = 0.50, method = "bc",
                   group = "treatment")
  expect_equal(res$rho, sqrt(0.20 * 0.50 / (0.50 * 0.80)), tolerance = 1e-12)
  b <- rho_bounds(p = 0.50, omega = 0.20)
  expect_equal(res$rho, b$rho_upper, tolerance = 1e-12)

  # Every point of the interval reproduces the reported rate
  lb <- latent_bounds(p_obs = 0.50, N_randomized = 100, N_dropout = 20,
                      method = "bc", group = "treatment")
  expect_equal(lb$p_lower, 0.30, tolerance = 1e-12)
  expect_equal(lb$p_upper, 0.50, tolerance = 1e-12)
  for (pc in c(0.32, 0.40, 0.50)) {
    r <- infer_rho(p_obs = 0.50, omega = 0.20, p = pc, method = "bc",
                   group = "treatment")$rho
    reported <- pc * (1 - 0.20) -
      r * sqrt(pc * (1 - pc) * 0.20 * 0.80) + 0.20
    expect_equal(reported, 0.50, tolerance = 1e-10)
  }
})

test_that("the two rules that count dropouts as responders agree", {
  # Each rule places the reported rate on a different part of the NRI scale, so
  # each has its own identified interval and the assumed p has to be taken from
  # the right one. Counting dropouts as responders gives [0.30, 0.50], counting
  # them as non-responders gives [0.50, 0.70]. Passing a p from the wrong
  # interval is exactly what produces an unattainable correlation, so these
  # calls also assert that no warning is raised.
  bc_t <- expect_no_warning(
    infer_rho(p_obs = 0.50, omega = 0.20, p = 0.40, method = "bc",
              group = "treatment")$rho
  )
  wc_c <- expect_no_warning(
    infer_rho(p_obs = 0.50, omega = 0.20, p = 0.40, method = "wc",
              group = "control")$rho
  )
  expect_equal(bc_t, wc_c, tolerance = 1e-14)

  # The remaining two count dropouts as non-responders, as NRI does
  nri <- expect_no_warning(
    infer_rho(p_obs = 0.50, omega = 0.20, p = 0.60, method = "nri")$rho
  )
  bc_c <- expect_no_warning(
    infer_rho(p_obs = 0.50, omega = 0.20, p = 0.60, method = "bc",
              group = "control")$rho
  )
  wc_t <- expect_no_warning(
    infer_rho(p_obs = 0.50, omega = 0.20, p = 0.60, method = "wc",
              group = "treatment")$rho
  )
  expect_equal(bc_c, nri, tolerance = 1e-14)
  expect_equal(wc_t, nri, tolerance = 1e-14)
})

test_that("no candidate inside an identified interval is unattainable", {
  # A warning here would mean the interval and the correlation map disagree
  for (m in c("nri", "cc", "bc", "wc")) {
    for (g in c("treatment", "control")) {
      lb <- latent_bounds(p_obs = 0.50, N_randomized = 100, N_dropout = 20,
                          method = m, group = g)
      grid <- seq(lb$p_lower, lb$p_upper, length.out = 9)
      grid <- grid[grid > 0 & grid < 1]
      expect_no_warning(
        for (pc in grid) {
          infer_rho(p_obs = 0.50, omega = lb$omega, p = pc,
                    method = m, group = g)
        }
      )
    }
  }
})

test_that("an assumed p is required under every rule", {
  for (m in c("nri", "cc", "bc", "wc")) {
    for (g in c("treatment", "control")) {
      expect_error(infer_rho(p_obs = 0.40, omega = 0.15, method = m,
                             group = g),
                   "latent_bounds", info = paste(m, g))
    }
  }
})

test_that("an unattainable combination is flagged", {
  expect_warning(
    infer_rho(p_obs = 0.05, omega = 0.10, p = 0.90, method = "nri"),
    "outside the feasible range"
  )
})

test_that("no dropout gives rho = 0", {
  res <- infer_rho(p_obs = 0.40, omega = 0, method = "nri")
  expect_equal(res$rho, 0)
  expect_equal(res$p, 0.40)
  expect_error(infer_rho(p_obs = 0.40, omega = 0, p = 0.50), "must equal p")
})

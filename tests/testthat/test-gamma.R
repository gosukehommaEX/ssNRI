# The gamma parameterization is an exact reexpression of the correlation
# parameterization, so every test here checks an identity rather than a
# tolerance-limited approximation. The reference relations were also verified
# independently in Python over a 98 x 99 grid of (p, omega).
#
# The grid below is swept once and the results are reduced to worst-case
# discrepancies, which are then asserted. Asserting cell by cell would multiply
# the same evidence into thousands of expectations and dominate the run time of
# the whole test suite; the informative part of a failure is the worst cell, so
# each label carries it.

GRID <- expand.grid(p = seq(0.05, 0.95, 0.05), omega = seq(0.05, 0.95, 0.05))

SWEEP <- local({
  n <- nrow(GRID)
  round_trip <- numeric(n)
  bound_lower <- numeric(n)
  bound_upper <- numeric(n)
  binding_cell <- numeric(n)
  negative_cell <- numeric(n)
  ordered <- logical(n)

  for (i in seq_len(n)) {
    p <- GRID$p[i]
    w <- GRID$omega[i]
    rb <- rho_bounds(p = p, omega = w)
    gb <- gamma_bounds(p = p, omega = w)

    # rho -> gamma -> rho across the whole feasible range
    rho <- seq(rb$rho_lower, rb$rho_upper, length.out = 5)
    g <- pmin(pmax(rho_to_gamma(p = p, omega = w, rho = rho), 0), 1)
    round_trip[i] <- max(abs(gamma_to_rho(p = p, omega = w, gamma = g) - rho))

    # the gamma bounds carry the same information as the Prentice bounds
    bound_lower[i] <- abs(gb$rho_lower - rb$rho_lower)
    bound_upper[i] <- abs(gb$rho_upper - rb$rho_upper)

    # at either endpoint some joint cell is exhausted, and none is negative
    cells <- vapply(c(gb$gamma_lower, gb$gamma_upper), function(gg) {
      c(gg * w, (1 - gg) * w, p - gg * w, 1 - p - w + gg * w)
    }, numeric(4))
    binding_cell[i] <- max(abs(apply(cells, 2, min)))
    negative_cell[i] <- -min(cells)
    ordered[i] <- gb$gamma_lower <= gb$gamma_upper
  }

  worst <- function(v) {
    k <- which.max(v)
    sprintf("worst at p = %g, omega = %g (%.3e)",
            GRID$p[k], GRID$omega[k], v[k])
  }
  list(round_trip = round_trip, bound_lower = bound_lower,
       bound_upper = bound_upper, binding_cell = binding_cell,
       negative_cell = negative_cell, ordered = ordered, worst = worst)
})

test_that("rho and gamma are inverse maps of one another", {
  expect_lt(max(SWEEP$round_trip), 1e-12,
            label = paste("round trip error,", SWEEP$worst(SWEEP$round_trip)))
})

test_that("gamma equals the latent response probability under independence", {
  expect_equal(rho_to_gamma(p = 0.45, omega = 0.10, rho = 0), 0.45)
  expect_equal(rho_to_gamma(p = 0.20, omega = 0.40, rho = 0), 0.20)
  expect_equal(gamma_to_rho(p = 0.45, omega = 0.10, gamma = 0.45), 0)
})

test_that("gamma_bounds reproduces the Prentice correlation bounds", {
  expect_lt(max(SWEEP$bound_lower), 1e-12,
            label = paste("lower bound,", SWEEP$worst(SWEEP$bound_lower)))
  expect_lt(max(SWEEP$bound_upper), 1e-12,
            label = paste("upper bound,", SWEEP$worst(SWEEP$bound_upper)))
})

test_that("the gamma bounds are exactly where a joint cell hits zero", {
  # An endpoint either exhausts one of the four cells or is pinned at 0 or 1,
  # and gamma * omega and (1 - gamma) * omega are themselves cells, so the
  # smallest cell is zero at both endpoints in every case
  expect_lt(max(SWEEP$binding_cell), 1e-12,
            label = paste("binding cell,", SWEEP$worst(SWEEP$binding_cell)))
  expect_lt(max(SWEEP$negative_cell), 1e-12,
            label = paste("negative cell,", SWEEP$worst(SWEEP$negative_cell)))
  expect_true(all(SWEEP$ordered))
})

test_that("the branch of each gamma bound switches where expected", {
  # Upper bound: p / omega binds once omega exceeds p
  expect_equal(gamma_bounds(p = 0.45, omega = 0.40)$gamma_upper, 1)
  expect_equal(gamma_bounds(p = 0.45, omega = 0.50)$gamma_upper, 0.45 / 0.50,
               tolerance = 1e-14)
  # Lower bound: (p + omega - 1) / omega binds once p + omega exceeds 1
  expect_equal(gamma_bounds(p = 0.45, omega = 0.50)$gamma_lower, 0)
  expect_equal(gamma_bounds(p = 0.45, omega = 0.60)$gamma_lower,
               0.05 / 0.60, tolerance = 1e-14)
})

test_that("the NRI estimand written in gamma matches the correlation form", {
  cfg <- expand.grid(p = c(0.30, 0.45, 0.70),
                     omega = c(0.05, 0.15, 0.30),
                     rho = c(-0.2, 0, 0.35))
  s <- sqrt(cfg$p * (1 - cfg$p) * cfg$omega * (1 - cfg$omega))
  g <- rho_to_gamma(p = cfg$p, omega = cfg$omega, rho = cfg$rho)
  expect_lt(max(abs((cfg$p - g * cfg$omega) -
                      (cfg$p * (1 - cfg$omega) - cfg$rho * s))), 1e-14)
})

test_that("gamma reproduces the NRI response probabilities of power_nri", {
  out <- power_nri(n1 = 60, n2 = 60, p1 = 0.45, p2 = 0.30,
                   omega1 = 0.12, omega2 = 0.20,
                   rho1 = 0.25, rho2 = -0.15)
  g1 <- rho_to_gamma(p = 0.45, omega = 0.12, rho = 0.25)
  g2 <- rho_to_gamma(p = 0.30, omega = 0.20, rho = -0.15)
  expect_equal(out$p1_NRI, 0.45 - g1 * 0.12, tolerance = 1e-14)
  expect_equal(out$p2_NRI, 0.30 - g2 * 0.20, tolerance = 1e-14)
  expect_equal(out$effect_NRI, (0.45 - g1 * 0.12) - (0.30 - g2 * 0.20),
               tolerance = 1e-14)
})

test_that("the conversions are vectorized and recycle", {
  g <- rho_to_gamma(p = 0.45, omega = 0.10, rho = c(-0.3, 0, 0.3))
  expect_length(g, 3)
  expect_equal(g[2], 0.45)
  r <- gamma_to_rho(p = c(0.4, 0.6), omega = c(0.1, 0.2),
                    gamma = c(0.4, 0.6))
  expect_equal(r, c(0, 0))
})

test_that("the gamma functions reject invalid input", {
  expect_error(rho_to_gamma(p = 0, omega = 0.1, rho = 0), "p must be")
  expect_error(rho_to_gamma(p = 0.4, omega = 0, rho = 0), "omega must be")
  expect_error(rho_to_gamma(p = 0.4, omega = 0.1, rho = NA), "finite")
  expect_error(gamma_to_rho(p = 0.4, omega = 0.1, gamma = 1.2), "gamma must be")
  expect_error(gamma_bounds(p = 0.4, omega = 1), "omega must be")
  expect_error(gamma_bounds(p = c(0.4, 0.5), omega = 0.1), "length one")
})

test_that("gamma_bounds has working print and plot methods", {
  gb <- gamma_bounds(p = 0.45, omega = 0.10)
  expect_output(print(gb), "Feasible Range of the Dropout Response")
  expect_s3_class(gb, "gamma_bounds")
  skip_if_not_installed("ggplot2")
  pl <- plot(gb, n_points = 11)
  expect_s3_class(pl, "ggplot")
})

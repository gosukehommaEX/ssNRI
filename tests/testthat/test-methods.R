# Smoke tests for the print and plot methods.

test_that("ggplot2 is available, so the plot tests below are not skipped", {
  # ggplot2 is a hard dependency. Asserting it here turns a missing ggplot2
  # into a failure instead of a silent skip, which is the reason the plot
  # tests can safely use skip_if_not_installed() without hiding anything.
  expect_true(requireNamespace("ggplot2", quietly = TRUE))
})

pcc <- power_cc(n1 = 30, n0 = 30, pi1 = 0.6, pi0 = 0.4,
                omega1 = 0.2, omega0 = 0.2, rho1 = 0.1, rho0 = 0.1)
pnri <- power_nri(n1 = 30, n0 = 30, pi1 = 0.6, pi0 = 0.4,
                  omega1 = 0.2, omega0 = 0.2)
sscc <- sample_size_cc(pi1 = 0.6, pi0 = 0.4, omega1 = 0.2, omega0 = 0.2)
ssnri <- sample_size_nri(pi1 = 0.6, pi0 = 0.4, omega1 = 0.2, omega0 = 0.2)
rb <- rho_bounds(pi = 0.4, omega = 0.15)
lb <- latent_bounds(pi_obs = 0.40, N_randomized = 200, N_dropout = 30)
ir <- infer_rho(pi_obs = 0.40, omega = 0.15, pi = 0.45, method = "nri")
gb <- gamma_bounds(pi = 0.4, omega = 0.15)

objects <- list(power_cc = pcc, power_nri = pnri,
                sample_size_cc = sscc, sample_size_nri = ssnri,
                rho_bounds = rb, latent_bounds = lb, infer_rho = ir,
                gamma_bounds = gb)

test_that("every constructor returns its own S3 class", {
  for (nm in names(objects)) {
    expect_s3_class(objects[[nm]], nm)
    expect_s3_class(objects[[nm]], "data.frame")
    expect_equal(nrow(objects[[nm]]), 1L, info = nm)
  }
})

test_that("print methods emit output and return their argument invisibly", {
  for (nm in names(objects)) {
    out <- capture.output(res <- print(objects[[nm]]))
    expect_gt(length(out), 3)
    expect_identical(res, objects[[nm]], info = nm)
    expect_false(any(grepl("NA", out, fixed = TRUE)), info = nm)
  }
})

test_that("print is what autoprinting dispatches to", {
  for (nm in names(objects)) {
    direct <- capture.output(print(objects[[nm]]))
    auto <- capture.output(objects[[nm]])
    expect_identical(direct, auto, info = nm)
  }
})

test_that("plot methods build a ggplot without error", {
  skip_if_not_installed("ggplot2")
  plots <- list(
    "power_cc, rho"        = plot(pcc, vary = "rho", n_points = 5),
    "power_cc, omega"      = plot(pcc, vary = "omega", n_points = 5),
    "power_nri, rho"       = plot(pnri, vary = "rho", n_points = 5),
    "power_nri, omega"     = plot(pnri, vary = "omega", n_points = 5),
    "sample_size_nri, rho" = plot(ssnri, vary = "rho", n_points = 5),
    "sample_size_nri, om"  = plot(ssnri, vary = "omega", n_points = 5),
    "rho_bounds"           = plot(rb, n_points = 5),
    "latent_bounds"        = plot(lb, n_points = 5),
    "gamma_bounds"         = plot(gb, n_points = 5)
  )
  for (nm in names(plots)) {
    expect_s3_class(plots[[nm]], "ggplot")
    # Building is what catches a bad aesthetic or a missing column
    built <- ggplot2::ggplot_build(plots[[nm]])
    expect_true(nrow(built$data[[1]]) > 0, info = nm)
  }
})

test_that("plot methods reject a degenerate sweep", {
  skip_if_not_installed("ggplot2")
  no_dropout <- power_cc(n1 = 20, n0 = 20, pi1 = 0.6, pi0 = 0.4,
                         omega1 = 0, omega0 = 0)
  expect_error(plot(no_dropout, vary = "rho"), "single point")

  expect_message(
    identified <- latent_bounds(pi_obs = 0.4, N_randomized = 100,
                                N_dropout = 0)
  )
  expect_error(plot(identified), "nothing to plot")

  expect_error(plot(pcc, n_points = 2), "at least 3")
  expect_error(plot(rb, omega_range = c(0.5, 0.1)), "increasing")
})

test_that("the sample size plot puts the simple inflation curve below the NRI one", {
  skip_if_not_installed("ggplot2")
  built <- ggplot2::ggplot_build(plot(ssnri, vary = "omega", n_points = 7))
  dat <- built$plot$data
  nri <- dat$n_total[dat$method == "NRI formula"]
  cc <- dat$n_total[dat$method == "Simple inflation"]
  expect_equal(length(nri), length(cc))
  expect_true(all(cc < nri))
})

# Internal helpers shared by the exact evaluation and plotting functions.

#' Smallest group 1 count that rejects, for every group 0 count
#'
#' For fixed group sizes m1 and m0 and a fixed group 0 count x0, the one-sided
#' pooled-variance Z statistic is strictly increasing in the group 1 count x1
#' wherever it is positive. Writing the statistic as
#' \code{(x1 - x0 m1 / m0) / sqrt((x1 + x0)(M - x1 - x0))} up to a positive
#' factor, with \code{M = m1 + m0}, the derivative of its logarithm is
#' \code{1 / (x1 - x0 m1 / m0) - (1 / 2)(1 / (x1 + x0) - 1 / (M - x1 - x0))},
#' and the first term dominates because \code{x1 - x0 m1 / m0 < x1 + x0}.
#' Rejection requires a positive statistic, so the rejection region is an upper
#' set in x1 and is described by a single threshold for each x0.
#'
#' The threshold is obtained as the larger root of the quadratic inequality
#' equivalent to \code{Z > z_alpha} and is then verified against the defining
#' inequality, so it is exact rather than approximate.
#'
#' All arithmetic is carried out in double precision. The quadratic
#' coefficients involve products of up to five counts, which overflows R's
#' 32-bit integer type from group sizes of about 64 upwards.
#'
#' @param m1 Numeric scalar. Group 1 size, m1 >= 1.
#' @param m0 Numeric scalar. Group 0 size, m0 >= 1.
#' @param z_alpha Numeric scalar. Critical value of the one-sided test.
#'
#' @return A numeric vector of length \code{m0 + 1}. Element \code{j + 1} is the
#'   smallest x1 for which the test rejects when x0 equals j, or \code{m1 + 1}
#'   when no x1 rejects.
#'
#' @noRd
.reject_threshold <- function(m1, m0, z_alpha) {

  m1 <- as.numeric(m1)
  m0 <- as.numeric(m0)
  M <- m1 + m0
  inv <- 1 / m1 + 1 / m0
  x0 <- as.numeric(0:m0)
  z2 <- z_alpha * z_alpha

  # Quadratic form of Z > z_alpha in x1. The leading coefficient is positive,
  # so the rejection region is the upper branch beyond the larger root
  qa <- M * m0 * m0 + z2 * m1 * m0
  qb <- -2 * x0 * M * m1 * m0 - z2 * m1 * m0 * (M - 2 * x0)
  qc <- M * x0 * x0 * m1 * m1 - z2 * m1 * m0 * x0 * (M - x0)
  disc <- qb * qb - 4 * qa * qc
  root <- ifelse(disc < 0, -qb / (2 * qa),
                 (-qb + sqrt(pmax(disc, 0))) / (2 * qa))

  cth <- pmax(ceiling(root - 1e-9), 0)
  cth <- pmax(cth, floor(x0 * m1 / m0))

  # Verify against the defining inequality. The root is accurate to well within
  # one unit, so this settles in a single step, but it is kept so that the
  # thresholds are exact by construction rather than by numerical luck
  z_at <- function(u) {
    u <- pmin(pmax(u, 0), m1)
    pp <- (u + x0) / M
    se <- sqrt(pp * (1 - pp) * inv)
    se[se < 1e-10] <- 1e-10
    (u / m1 - x0 / m0) / se
  }
  repeat {
    move <- (cth <= m1) & (z_at(cth) <= z_alpha)
    if (!any(move)) break
    cth[move] <- cth[move] + 1
  }
  repeat {
    move <- (cth >= 1) & (z_at(cth - 1) > z_alpha)
    if (!any(move)) break
    cth[move] <- cth[move] - 1
  }
  cth[cth > m1] <- m1 + 1

  return(cth)
}

#' Bivariate Bernoulli joint cell probabilities
#'
#' @param pi Numeric. Marginal probability that the latent response equals one.
#' @param omega Numeric. Marginal dropout probability.
#' @param rho Numeric. Correlation between the two indicators.
#'
#' @return A named numeric vector with elements \code{pi_10}, \code{pi_11},
#'   \code{pi_00} and \code{pi_01}, clamped to the unit interval.
#'
#' @noRd
.joint_cells <- function(pi, omega, rho) {
  s <- sqrt(pi * (1 - pi) * omega * (1 - omega))
  cells <- c(
    pi_10 = pi * (1 - omega) - rho * s,
    pi_11 = pi * omega + rho * s,
    pi_00 = (1 - pi) * (1 - omega) + rho * s,
    pi_01 = (1 - pi) * omega - rho * s
  )
  # Clamp away tiny negative values at the edge of the feasible region. pmax()
  # and pmin() copy attributes from their FIRST argument, so cells has to come
  # first or the names are silently dropped; they are restored explicitly as
  # well, because callers select cells by name
  out <- pmin(pmax(cells, 0), 1)
  names(out) <- names(cells)
  return(out)
}

#' Correlation range feasible in both groups simultaneously
#'
#' @param pi1,pi0 Numeric. Latent response probabilities.
#' @param omega1,omega0 Numeric. Dropout probabilities.
#'
#' @return A numeric vector of length two, the intersection of the two feasible
#'   ranges.
#'
#' @noRd
.common_rho_range <- function(pi1, pi0, omega1, omega0) {
  b1 <- rho_bounds(pi = pi1, omega = omega1)
  b0 <- rho_bounds(pi = pi0, omega = omega0)
  c(max(b1$rho_lower, b0$rho_lower), min(b1$rho_upper, b0$rho_upper))
}

#' Move a correlation into its feasible range
#'
#' @param rho Numeric. Correlation to clamp.
#' @param pi Numeric. Latent response probability.
#' @param omega Numeric. Dropout probability.
#'
#' @return The clamped correlation.
#'
#' @noRd
.clamp_rho <- function(rho, pi, omega) {
  b <- rho_bounds(pi = pi, omega = omega)
  min(max(rho, b$rho_lower), b$rho_upper)
}

#' Stop when ggplot2 is not installed
#'
#' @param what Character. Name of the calling method, used in the message.
#'
#' @return Invisibly \code{TRUE}; called for the side effect of stopping.
#'
#' @noRd
.require_ggplot2 <- function(what) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(what, " requires the ggplot2 package. ",
         "Install it with install.packages(\"ggplot2\").", call. = FALSE)
  }
  invisible(TRUE)
}

#' Sensitivity plot shared by the power classes
#'
#' @param x An object of class \code{power_cc} or \code{power_nri}.
#' @param vary Character. Either \code{"rho"} or \code{"omega"}.
#' @param n_points Integer. Number of grid points.
#' @param fun The evaluation function to call at each grid point.
#' @param label Character. Title prefix.
#'
#' @return A \code{ggplot} object.
#'
#' @noRd
.plot_power <- function(x, vary, n_points, fun, label) {

  .require_ggplot2("The plot method")
  if (n_points < 3) stop("n_points must be at least 3")

  if (vary == "rho") {
    rng <- .common_rho_range(x$pi1, x$pi0, x$omega1, x$omega0)
    if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[2] - rng[1] < 1e-8) {
      stop("The correlation range feasible in both groups is a single point; ",
           "use vary = \"omega\" instead", call. = FALSE)
    }
    grid <- seq(rng[1], rng[2], length.out = n_points)
    res <- lapply(grid, function(g) {
      fun(n1 = x$n1, n0 = x$n0, pi1 = x$pi1, pi0 = x$pi0,
          omega1 = x$omega1, omega0 = x$omega0,
          rho1 = g, rho0 = g, alpha = x$alpha)
    })
    x_lab <- expression(rho)
    marker <- if (isTRUE(all.equal(x$rho1, x$rho0))) x$rho1 else NA_real_
    sub <- parse(text = sprintf(
      'paste("(", omega[1], ", ", omega[0], ") = (%.2f, %.2f)")',
      x$omega1, x$omega0
    ))
  } else {
    grid <- seq(0.01, 0.50, length.out = n_points)
    res <- lapply(grid, function(g) {
      fun(n1 = x$n1, n0 = x$n0, pi1 = x$pi1, pi0 = x$pi0,
          omega1 = g, omega0 = g,
          rho1 = .clamp_rho(x$rho1, x$pi1, g),
          rho0 = .clamp_rho(x$rho0, x$pi0, g),
          alpha = x$alpha)
    })
    x_lab <- expression(omega)
    marker <- if (isTRUE(all.equal(x$omega1, x$omega0))) x$omega1 else NA_real_
    sub <- parse(text = sprintf(
      paste0('paste("(", rho[1], ", ", rho[0], ") = (%.2f, %.2f), ',
             'moved into the feasible range where needed")'),
      x$rho1, x$rho0
    ))
  }

  dat <- data.frame(
    xval = rep(grid, 2),
    value = c(vapply(res, function(z) z$power, numeric(1)),
              vapply(res, function(z) z$type1_error, numeric(1))),
    quantity = factor(rep(c("Power", "Type I error rate"),
                          each = length(grid)),
                      levels = c("Power", "Type I error rate"))
  )

  out <- ggplot2::ggplot(
    dat, ggplot2::aes(x = .data$xval, y = .data$value)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::facet_wrap(ggplot2::vars(.data$quantity), scales = "free_y") +
    ggplot2::labs(
      x = x_lab, y = NULL,
      title = parse(text = sprintf(
        paste0('paste("%s, (", n[1], ", ", n[0], ") = (%d, %d), ',
               '(", pi[1], ", ", pi[0], ") = (%.2f, %.2f)")'),
        label, x$n1, x$n0, x$pi1, x$pi0
      )),
      subtitle = sub
    ) +
    ggplot2::theme_bw(base_size = 12)

  if (is.finite(marker)) {
    out <- out + ggplot2::geom_vline(xintercept = marker, linetype = "dashed")
  }

  out
}

#' Sample Size for a Trial Analyzed under Non-Responder Imputation
#'
#' Computes the required sample size for a two-group superiority trial with a
#' binary endpoint whose primary analysis applies non-responder imputation
#' (NRI), so that every dropout counts as a non-responder. The calculation is
#' carried out on the NRI scale, which is the scale the analysis actually tests.
#'
#' @param p1 Numeric. Latent response probability for group 1, 0 < p1 < 1.
#' @param p2 Numeric. Latent response probability for group 2, 0 < p2 < 1.
#' @param omega1 Numeric. Dropout probability for group 1, 0 <= omega1 < 1
#'   (default = 0).
#' @param omega2 Numeric. Dropout probability for group 2, 0 <= omega2 < 1
#'   (default = 0).
#' @param rho1 Numeric. Correlation between the latent response and the dropout
#'   indicator in group 1 (default = 0).
#' @param rho2 Numeric. Correlation between the latent response and the dropout
#'   indicator in group 2 (default = 0).
#' @param r Numeric. Allocation ratio n1 / n2, r > 0 (default = 1).
#' @param alpha Numeric. One-sided significance level (default = 0.025).
#' @param target_power Numeric. Target power (default = 0.8).
#'
#' @return A \code{data.frame} of class \code{sample_size_nri} with one row
#' containing the inputs together with \code{p1_NRI}, \code{p2_NRI},
#' \code{effect_latent}, \code{effect_NRI}, \code{n1}, \code{n2} and
#' \code{n_total}.
#'
#' @details
#' The response probability seen by an NRI analysis is the joint probability of
#' responding and completing,
#'
#'   pj_NRI = pj (1 - omegaj) - rhoj sqrt(pj (1 - pj) omegaj (1 - omegaj)),
#'
#' and the effect the trial is powered against is delta_NRI = p1_NRI - p2_NRI.
#' Under a common dropout probability and within-group independence this is
#' (1 - omega) times the complete case effect, so the simple inflation method,
#' which powers against the larger complete case effect, falls short. Unequal
#' dropout probabilities or unequal correlations can move delta_NRI in either
#' direction relative to delta_CC.
#'
#' No inflation step is applied: the returned sample size is already the number
#' of patients to randomize, because dropouts remain in the analysis as
#' non-responders rather than being removed from it.
#'
#' @examples
#' sample_size_nri(p1 = 0.6, p2 = 0.4, omega1 = 0.2, omega2 = 0.2)
#'
#' # Compare with the simple inflation method under the same assumptions
#' sample_size_cc(p1 = 0.6, p2 = 0.4, omega1 = 0.2, omega2 = 0.2)$n_total_adjust
#'
#' @seealso \code{\link{sample_size_cc}}, \code{\link{power_nri}}
#'
#' @importFrom stats qnorm
#' @export
sample_size_nri <- function(p1, p2, omega1 = 0, omega2 = 0,
                            rho1 = 0, rho2 = 0,
                            r = 1, alpha = 0.025, target_power = 0.8) {

  # Input validation
  if (p1 <= 0 || p1 >= 1) stop("p1 must be between 0 and 1")
  if (p2 <= 0 || p2 >= 1) stop("p2 must be between 0 and 1")
  if (omega1 < 0 || omega1 >= 1) stop("omega1 must be between 0 and 1")
  if (omega2 < 0 || omega2 >= 1) stop("omega2 must be between 0 and 1")
  if (r <= 0) stop("r must be positive")
  if (alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1")
  if (target_power <= 0 || target_power >= 1) {
    stop("target_power must be between 0 and 1")
  }
  if (p1 <= p2) warning("p1 should be greater than p2 for superiority trial")

  # Validate correlation bounds
  ctol <- 1e-6
  bounds1 <- rho_bounds(p = p1, omega = omega1)
  bounds2 <- rho_bounds(p = p2, omega = omega2)

  if (rho1 < bounds1$rho_lower - ctol || rho1 > bounds1$rho_upper + ctol) {
    stop(paste0("rho1 must be between ", round(bounds1$rho_lower, 4),
                " and ", round(bounds1$rho_upper, 4)))
  }
  if (rho2 < bounds2$rho_lower - ctol || rho2 > bounds2$rho_upper + ctol) {
    stop(paste0("rho2 must be between ", round(bounds2$rho_lower, 4),
                " and ", round(bounds2$rho_upper, 4)))
  }

  # Response probabilities on the NRI scale
  p1_NRI <- p1 * (1 - omega1) -
    rho1 * sqrt(p1 * (1 - p1) * omega1 * (1 - omega1))
  p2_NRI <- p2 * (1 - omega2) -
    rho2 * sqrt(p2 * (1 - p2) * omega2 * (1 - omega2))

  prob_tol <- 1e-10
  if (p1_NRI < -prob_tol || p1_NRI > 1 + prob_tol) {
    stop("Calculated p1_NRI is outside valid range. Check input parameters.")
  }
  if (p2_NRI < -prob_tol || p2_NRI > 1 + prob_tol) {
    stop("Calculated p2_NRI is outside valid range. Check input parameters.")
  }

  effect_NRI <- p1_NRI - p2_NRI
  if (abs(effect_NRI) < 1e-12) {
    stop("The NRI treatment effect is zero, so no finite sample size attains ",
         "the target power. Check p1, p2, omega1, omega2, rho1 and rho2.")
  }

  # Closed-form sample size on the NRI scale
  p_pooled <- (r * p1_NRI + p2_NRI) / (1 + r)
  var_h0 <- p_pooled * (1 - p_pooled)
  var_h1 <- (p1_NRI * (1 - p1_NRI) / r + p2_NRI * (1 - p2_NRI)) / (1 + 1 / r)

  n2 <- (1 + 1 / r) / (effect_NRI ^ 2) *
    (qnorm(1 - alpha) * sqrt(var_h0) + qnorm(target_power) * sqrt(var_h1)) ^ 2

  n2 <- ceiling(n2)
  n1 <- ceiling(r * n2)
  n_total <- n1 + n2

  result <- data.frame(
    p1 = p1,
    p2 = p2,
    omega1 = omega1,
    omega2 = omega2,
    rho1 = rho1,
    rho2 = rho2,
    r = r,
    alpha = alpha,
    target_power = target_power,
    p1_NRI = p1_NRI,
    p2_NRI = p2_NRI,
    effect_latent = p1 - p2,
    effect_NRI = effect_NRI,
    n1 = n1,
    n2 = n2,
    n_total = n_total
  )

  class(result) <- c("sample_size_nri", "data.frame")

  return(result)
}

#' Print Method for sample_size_nri Objects
#'
#' @param x An object of class \code{sample_size_nri}.
#' @param ... Additional arguments (not used).
#'
#' @return \code{x}, invisibly.
#'
#' @export
print.sample_size_nri <- function(x, ...) {
  cat("Sample Size for an NRI Analysis\n")
  cat("===============================\n")
  cat("Design parameters:\n")
  cat(sprintf("             p1 = %.3f\n", x$p1))
  cat(sprintf("             p2 = %.3f\n", x$p2))
  cat(sprintf("         omega1 = %.3f\n", x$omega1))
  cat(sprintf("         omega2 = %.3f\n", x$omega2))
  cat(sprintf("           rho1 = %.3f\n", x$rho1))
  cat(sprintf("           rho2 = %.3f\n", x$rho2))
  cat(sprintf("              r = %.4g\n", x$r))
  cat(sprintf("          alpha = %.4f\n", x$alpha))
  cat(sprintf("   target_power = %.2f\n", x$target_power))

  cat("\nNRI response probabilities:\n")
  cat(sprintf("         p1_NRI = %.4f\n", x$p1_NRI))
  cat(sprintf("         p2_NRI = %.4f\n", x$p2_NRI))
  cat(sprintf("  effect_latent = %.4f\n", x$effect_latent))
  cat(sprintf("     effect_NRI = %.4f\n", x$effect_NRI))

  cat("\nRequired sample size (patients to randomize):\n")
  cat(sprintf("             n1 = %d\n", x$n1))
  cat(sprintf("             n2 = %d\n", x$n2))
  cat(sprintf("        n_total = %d\n", x$n_total))

  invisible(x)
}

#' Plot Method for sample_size_nri Objects
#'
#' Required total sample size under the NRI formula and under the simple
#' inflation method, as a function of one design input, with the remaining
#' inputs held at the values stored in the object. The vertical distance between
#' the two curves is the shortfall the simple inflation method incurs. A dashed
#' vertical line marks the value the object itself was evaluated at.
#'
#' @param x An object of class \code{sample_size_nri}.
#' @param vary Character. Design input to vary: \code{"rho"} (default) sets a
#'   common correlation in both groups and sweeps the range feasible in both, or
#'   \code{"omega"} sets a common dropout probability and sweeps 0.01 to 0.50,
#'   moving the correlations into the feasible range where a value would
#'   otherwise be infeasible.
#' @param n_points Integer. Number of grid points (default = 41).
#' @param ... Additional arguments (not used).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' ss <- sample_size_nri(p1 = 0.6, p2 = 0.4, omega1 = 0.2, omega2 = 0.2)
#' plot(ss, vary = "omega", n_points = 11)
#'
#' @importFrom rlang .data
#' @export
plot.sample_size_nri <- function(x, vary = c("rho", "omega"),
                                 n_points = 41, ...) {

  vary <- match.arg(vary)
  .require_ggplot2("The plot method")
  if (n_points < 3) stop("n_points must be at least 3")

  if (vary == "rho") {
    rng <- .common_rho_range(x$p1, x$p2, x$omega1, x$omega2)
    if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[2] - rng[1] < 1e-8) {
      stop("The correlation range feasible in both groups is a single point; ",
           "use vary = \"omega\" instead", call. = FALSE)
    }
    grid <- seq(rng[1], rng[2], length.out = n_points)
    omega_of <- function(g) c(x$omega1, x$omega2)
    rho_of <- function(g) c(g, g)
    x_lab <- expression(rho)
    marker <- if (isTRUE(all.equal(x$rho1, x$rho2))) x$rho1 else NA_real_
  } else {
    grid <- seq(0.01, 0.50, length.out = n_points)
    omega_of <- function(g) c(g, g)
    rho_of <- function(g) c(.clamp_rho(x$rho1, x$p1, g),
                            .clamp_rho(x$rho2, x$p2, g))
    x_lab <- expression(omega)
    marker <- if (isTRUE(all.equal(x$omega1, x$omega2))) x$omega1 else NA_real_
  }

  n_nri <- numeric(length(grid))
  n_cc <- numeric(length(grid))
  for (i in seq_along(grid)) {
    om <- omega_of(grid[i])
    rh <- rho_of(grid[i])
    # The NRI effect can vanish at the edge of the feasible correlation range,
    # in which case no finite sample size exists and the point is left blank
    n_nri[i] <- tryCatch(
      sample_size_nri(
        p1 = x$p1, p2 = x$p2, omega1 = om[1], omega2 = om[2],
        rho1 = rh[1], rho2 = rh[2], r = x$r, alpha = x$alpha,
        target_power = x$target_power
      )$n_total,
      error = function(e) NA_real_
    )
    n_cc[i] <- tryCatch(
      sample_size_cc(
        p1 = x$p1, p2 = x$p2, omega1 = om[1], omega2 = om[2],
        r = x$r, alpha = x$alpha, target_power = x$target_power
      )$n_total_adjust,
      error = function(e) NA_real_
    )
  }

  dat <- data.frame(
    xval = rep(grid, 2),
    n_total = c(n_nri, n_cc),
    method = factor(rep(c("NRI formula", "Simple inflation"),
                        each = length(grid)),
                    levels = c("NRI formula", "Simple inflation"))
  )

  out <- ggplot2::ggplot(
    dat, ggplot2::aes(x = .data$xval, y = .data$n_total,
                      colour = .data$method, linetype = .data$method)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(
      x = x_lab, y = "Required total sample size",
      colour = NULL, linetype = NULL,
      title = sprintf("p = (%.2f, %.2f), r = %.4g, target power = %.2f",
                      x$p1, x$p2, x$r, x$target_power)
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "top")

  if (is.finite(marker)) {
    out <- out + ggplot2::geom_vline(xintercept = marker, linetype = "dashed")
  }

  out
}

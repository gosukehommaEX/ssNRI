#' Sample Size for a Trial Analyzed under Non-Responder Imputation
#'
#' Computes the required sample size for a two-group superiority trial with a
#' binary endpoint whose primary analysis applies non-responder imputation
#' (NRI), so that every dropout counts as a non-responder. The calculation is
#' carried out on the NRI scale, which is the scale the analysis actually tests.
#'
#' @param pi1 Numeric. Latent response probability for group 1, 0 < pi1 < 1.
#' @param pi0 Numeric. Latent response probability for group 0, 0 < pi0 < 1.
#' @param omega1 Numeric. Dropout probability for group 1, 0 <= omega1 < 1
#'   (default = 0).
#' @param omega0 Numeric. Dropout probability for group 0, 0 <= omega0 < 1
#'   (default = 0).
#' @param rho1 Numeric. Correlation between the latent response and the dropout
#'   indicator in group 1 (default = 0).
#' @param rho0 Numeric. Correlation between the latent response and the dropout
#'   indicator in group 0 (default = 0).
#' @param r Numeric. Allocation ratio n1 / n0, r > 0 (default = 1).
#' @param alpha Numeric. One-sided significance level (default = 0.025).
#' @param target_power Numeric. Target power (default = 0.8).
#'
#' @return A \code{data.frame} of class \code{sample_size_nri} with one row
#' containing the inputs together with \code{pi1_NRI}, \code{pi0_NRI},
#' \code{effect_full}, \code{effect_NRI}, \code{n1}, \code{n0} and
#' \code{n_total}.
#'
#' @details
#' The response probability seen by an NRI analysis is the joint probability of
#' responding and completing,
#'
#'   pi_j_NRI = pi_j (1 - omega_j) - rho_j sqrt(pi_j (1 - pi_j) omega_j (1 - omega_j)),
#'
#' and the effect the trial is powered against is delta_NRI = pi1_NRI - pi0_NRI.
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
#' sample_size_nri(pi1 = 0.6, pi0 = 0.4, omega1 = 0.2, omega0 = 0.2)
#'
#' # Compare with the simple inflation method under the same assumptions
#' sample_size_cc(pi1 = 0.6, pi0 = 0.4, omega1 = 0.2, omega0 = 0.2)$n_total_adjust
#'
#' @seealso \code{\link{sample_size_cc}}, \code{\link{power_nri}}
#'
#' @importFrom stats qnorm
#' @export
sample_size_nri <- function(pi1, pi0, omega1 = 0, omega0 = 0,
                            rho1 = 0, rho0 = 0,
                            r = 1, alpha = 0.025, target_power = 0.8) {

  # Input validation
  if (pi1 <= 0 || pi1 >= 1) stop("pi1 must be between 0 and 1")
  if (pi0 <= 0 || pi0 >= 1) stop("pi0 must be between 0 and 1")
  if (omega1 < 0 || omega1 >= 1) stop("omega1 must be between 0 and 1")
  if (omega0 < 0 || omega0 >= 1) stop("omega0 must be between 0 and 1")
  if (r <= 0) stop("r must be positive")
  if (alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1")
  if (target_power <= 0 || target_power >= 1) {
    stop("target_power must be between 0 and 1")
  }
  if (pi1 <= pi0) warning("pi1 should be greater than pi0 for superiority trial")

  # Validate correlation bounds
  ctol <- 1e-6
  bounds1 <- rho_bounds(pi = pi1, omega = omega1)
  bounds0 <- rho_bounds(pi = pi0, omega = omega0)

  if (rho1 < bounds1$rho_lower - ctol || rho1 > bounds1$rho_upper + ctol) {
    stop(paste0("rho1 must be between ", round(bounds1$rho_lower, 4),
                " and ", round(bounds1$rho_upper, 4)))
  }
  if (rho0 < bounds0$rho_lower - ctol || rho0 > bounds0$rho_upper + ctol) {
    stop(paste0("rho0 must be between ", round(bounds0$rho_lower, 4),
                " and ", round(bounds0$rho_upper, 4)))
  }

  # Response probabilities on the NRI scale
  pi1_NRI <- pi1 * (1 - omega1) -
    rho1 * sqrt(pi1 * (1 - pi1) * omega1 * (1 - omega1))
  pi0_NRI <- pi0 * (1 - omega0) -
    rho0 * sqrt(pi0 * (1 - pi0) * omega0 * (1 - omega0))

  prob_tol <- 1e-10
  if (pi1_NRI < -prob_tol || pi1_NRI > 1 + prob_tol) {
    stop("Calculated pi1_NRI is outside valid range. Check input parameters.")
  }
  if (pi0_NRI < -prob_tol || pi0_NRI > 1 + prob_tol) {
    stop("Calculated pi0_NRI is outside valid range. Check input parameters.")
  }

  effect_NRI <- pi1_NRI - pi0_NRI
  if (abs(effect_NRI) < 1e-12) {
    stop("The NRI treatment effect is zero, so no finite sample size attains ",
         "the target power. Check pi1, pi0, omega1, omega0, rho1 and rho0.")
  }

  # Closed-form sample size on the NRI scale
  p_pooled <- (r * pi1_NRI + pi0_NRI) / (1 + r)
  var_h0 <- p_pooled * (1 - p_pooled)
  var_h1 <- (pi1_NRI * (1 - pi1_NRI) / r + pi0_NRI * (1 - pi0_NRI)) / (1 + 1 / r)

  n0 <- (1 + 1 / r) / (effect_NRI ^ 2) *
    (qnorm(1 - alpha) * sqrt(var_h0) + qnorm(target_power) * sqrt(var_h1)) ^ 2

  n0 <- ceiling(n0)
  n1 <- ceiling(r * n0)
  n_total <- n1 + n0

  result <- data.frame(
    pi1 = pi1,
    pi0 = pi0,
    omega1 = omega1,
    omega0 = omega0,
    rho1 = rho1,
    rho0 = rho0,
    r = r,
    alpha = alpha,
    target_power = target_power,
    pi1_NRI = pi1_NRI,
    pi0_NRI = pi0_NRI,
    effect_full = pi1 - pi0,
    effect_NRI = effect_NRI,
    n1 = n1,
    n0 = n0,
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
  cat(sprintf("           pi1 = %.3f\n", x$pi1))
  cat(sprintf("           pi0 = %.3f\n", x$pi0))
  cat(sprintf("        omega1 = %.3f\n", x$omega1))
  cat(sprintf("        omega0 = %.3f\n", x$omega0))
  cat(sprintf("          rho1 = %.3f\n", x$rho1))
  cat(sprintf("          rho0 = %.3f\n", x$rho0))
  cat(sprintf("             r = %.4g\n", x$r))
  cat(sprintf("         alpha = %.4f\n", x$alpha))
  cat(sprintf("  target_power = %.2f\n", x$target_power))

  cat("\nNRI response probabilities:\n")
  cat(sprintf("       pi1_NRI = %.4f\n", x$pi1_NRI))
  cat(sprintf("       pi0_NRI = %.4f\n", x$pi0_NRI))
  cat(sprintf("   effect_full = %.4f\n", x$effect_full))
  cat(sprintf("    effect_NRI = %.4f\n", x$effect_NRI))

  cat("\nRequired sample size (patients to randomize):\n")
  cat(sprintf("            n1 = %d\n", x$n1))
  cat(sprintf("            n0 = %d\n", x$n0))
  cat(sprintf("       n_total = %d\n", x$n_total))

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
#' ss <- sample_size_nri(pi1 = 0.6, pi0 = 0.4, omega1 = 0.2, omega0 = 0.2)
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
    rng <- .common_rho_range(x$pi1, x$pi0, x$omega1, x$omega0)
    if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[2] - rng[1] < 1e-8) {
      stop("The correlation range feasible in both groups is a single point; ",
           "use vary = \"omega\" instead", call. = FALSE)
    }
    grid <- seq(rng[1], rng[2], length.out = n_points)
    omega_of <- function(g) c(x$omega1, x$omega0)
    rho_of <- function(g) c(g, g)
    x_lab <- expression(rho)
    marker <- if (isTRUE(all.equal(x$rho1, x$rho0))) x$rho1 else NA_real_
  } else {
    grid <- seq(0.01, 0.50, length.out = n_points)
    omega_of <- function(g) c(g, g)
    rho_of <- function(g) c(.clamp_rho(x$rho1, x$pi1, g),
                            .clamp_rho(x$rho0, x$pi0, g))
    x_lab <- expression(omega)
    marker <- if (isTRUE(all.equal(x$omega1, x$omega0))) x$omega1 else NA_real_
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
        pi1 = x$pi1, pi0 = x$pi0, omega1 = om[1], omega0 = om[2],
        rho1 = rh[1], rho0 = rh[2], r = x$r, alpha = x$alpha,
        target_power = x$target_power
      )$n_total,
      error = function(e) NA_real_
    )
    n_cc[i] <- tryCatch(
      sample_size_cc(
        pi1 = x$pi1, pi0 = x$pi0, omega1 = om[1], omega0 = om[2],
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
      title = parse(text = sprintf(
        paste0('paste("(", pi[1], ", ", pi[0], ") = (%.2f, %.2f), ',
               'r = %.4g, target power = %.2f")'),
        x$pi1, x$pi0, x$r, x$target_power
      ))
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "top")

  if (is.finite(marker)) {
    out <- out + ggplot2::geom_vline(xintercept = marker, linetype = "dashed")
  }

  out
}

#' Exact Power and Type I Error Rate for an NRI Analysis
#'
#' Evaluates the exact power and the exact type I error rate of a one-sided
#' pooled-variance Z-test applied to a non-responder imputation (NRI) analysis,
#' in which every dropout is counted as a non-responder.
#'
#' @param n1 Integer. Randomized sample size for group 1 (treatment), n1 > 0.
#' @param n0 Integer. Randomized sample size for group 0 (control), n0 > 0.
#' @param pi1 Numeric. Latent response probability for group 1, 0 < pi1 < 1.
#' @param pi0 Numeric. Latent response probability for group 0, 0 < pi0 < 1.
#' @param omega1 Numeric. Dropout probability for group 1, 0 <= omega1 < 1.
#' @param omega0 Numeric. Dropout probability for group 0, 0 <= omega0 < 1.
#' @param rho1 Numeric. Correlation between the latent response and the dropout
#'   indicator in group 1 (default = 0).
#' @param rho0 Numeric. Correlation between the latent response and the dropout
#'   indicator in group 0 (default = 0).
#' @param alpha Numeric. One-sided significance level, 0 < alpha < 1
#'   (default = 0.025).
#' @param pi_null Numeric or NULL. Common response probability under the null
#'   hypothesis, used only for the type I error rate. The default NULL uses the
#'   sample-size weighted average of the two NRI response probabilities. Supply
#'   a value to evaluate the size of the test at a different nuisance parameter.
#'
#' @return A \code{data.frame} of class \code{power_nri} with one row containing
#' the inputs together with:
#' \describe{
#'   \item{pi1_NRI}{Observed response probability in group 1 under NRI}
#'   \item{pi0_NRI}{Observed response probability in group 0 under NRI}
#'   \item{effect_full}{Full-data treatment effect, pi1 - pi0}
#'   \item{effect_NRI}{NRI treatment effect, pi1_NRI - pi0_NRI}
#'   \item{power}{Exact power}
#'   \item{type1_error}{Exact type I error rate}
#' }
#'
#' @details
#' Under NRI a patient counts as a responder only when the latent response
#' equals one and the patient completes the trial, so the number of observed
#' responders in group j follows a binomial distribution with size n_j and
#' success probability
#'
#'   pi_j_NRI = P(R_ij = 1, D_ij = 0)
#'          = pi_j (1 - omega_j) - rho_j sqrt(pi_j (1 - pi_j) omega_j (1 - omega_j)).
#'
#' The randomized sample size is fixed, so no completer-count stage is involved
#' and the evaluation is a single enumeration over the observed responder pair.
#' As in \code{\link{power_cc}}, the rejection region is an upper set in the
#' group 1 count for every fixed group 0 count, so the enumeration reduces to
#' upper binomial tail probabilities evaluated at one threshold per group 0
#' count. This is exact; nothing is discarded.
#'
#' The type I error rate is evaluated at a common response probability, which by
#' default is the sample-size weighted average
#' pi_bar = (n1 pi1_NRI + n0 pi0_NRI) / (n1 + n0). That choice is a convention rather
#' than a property of the test: the size of a two-sample binomial test depends on
#' the nuisance parameter, and \code{pi_null} exists so that the dependence can be
#' examined. The returned \code{pi_null} column records the value used.
#'
#' The NRI estimand differs from the complete case estimand. Writing
#' s_j = sqrt(pi_j (1 - pi_j) omega_j (1 - omega_j)),
#'
#'   delta_NRI = pi1 (1 - omega1) - pi0 (1 - omega0) - rho1 s1 + rho0 s0,
#'
#' which reduces to (1 - omega) delta_CC when the two groups share a dropout
#' probability and the latent response is independent of dropout within each
#' group. The two estimands can differ in sign as well as in magnitude, so a
#' non-null NRI contrast under equal latent response probabilities is not an
#' inflated type I error rate but a correctly sized test of a different null
#' hypothesis.
#'
#' @examples
#' power_nri(n1 = 60, n0 = 60, pi1 = 0.6, pi0 = 0.4,
#'           omega1 = 0.2, omega0 = 0.2, rho1 = 0, rho0 = 0)
#'
#' # Unequal dropout probabilities move the NRI estimand away from zero
#' # even when the latent response probabilities are equal
#' power_nri(n1 = 60, n0 = 60, pi1 = 0.4, pi0 = 0.4,
#'           omega1 = 0.05, omega0 = 0.20)
#'
#' @seealso \code{\link{power_cc}}, \code{\link{rho_bounds}}
#'
#' @importFrom stats dbinom qnorm
#' @export
power_nri <- function(n1, n0, pi1, pi0, omega1, omega0,
                      rho1 = 0, rho0 = 0, alpha = 0.025, pi_null = NULL) {

  # Input validation
  if (n1 <= 0 || n0 <= 0) stop("Sample sizes must be positive")
  if (n1 != round(n1) || n0 != round(n0)) stop("Sample sizes must be integers")
  if (pi1 <= 0 || pi1 >= 1) stop("pi1 must be between 0 and 1")
  if (pi0 <= 0 || pi0 >= 1) stop("pi0 must be between 0 and 1")
  if (omega1 < 0 || omega1 >= 1) stop("omega1 must be between 0 and 1")
  if (omega0 < 0 || omega0 >= 1) stop("omega0 must be between 0 and 1")
  if (alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1")
  if (!is.null(pi_null) && (length(pi_null) != 1 || is.na(pi_null) ||
                           pi_null <= 0 || pi_null >= 1)) {
    stop("pi_null must be a single value between 0 and 1, or NULL")
  }

  n1 <- as.integer(round(n1))
  n0 <- as.integer(round(n0))

  # Validate correlation bounds
  ctol <- 1e-6
  if (omega1 > 0) {
    bounds1 <- rho_bounds(pi = pi1, omega = omega1)
    if (rho1 < bounds1$rho_lower - ctol || rho1 > bounds1$rho_upper + ctol) {
      stop(paste0("rho1 must be between ", round(bounds1$rho_lower, 4),
                  " and ", round(bounds1$rho_upper, 4)))
    }
  }
  if (omega0 > 0) {
    bounds0 <- rho_bounds(pi = pi0, omega = omega0)
    if (rho0 < bounds0$rho_lower - ctol || rho0 > bounds0$rho_upper + ctol) {
      stop(paste0("rho0 must be between ", round(bounds0$rho_lower, 4),
                  " and ", round(bounds0$rho_upper, 4)))
    }
  }

  # Observed response probabilities under NRI
  pi1_NRI <- .joint_cells(pi = pi1, omega = omega1, rho = rho1)[["pi_10"]]
  pi0_NRI <- .joint_cells(pi = pi0, omega = omega0, rho = rho0)[["pi_10"]]

  effect_NRI <- pi1_NRI - pi0_NRI

  # Null value: pooled by default, otherwise as supplied
  pi_bar <- if (is.null(pi_null)) {
    (n1 * pi1_NRI + n0 * pi0_NRI) / (n1 + n0)
  } else {
    pi_null
  }

  z_alpha <- qnorm(1 - alpha)

  # Rejection threshold in x1 for every x0. The randomized sample sizes are
  # fixed, so a single set of thresholds serves both hypotheses
  idx <- .reject_threshold(n1, n0, z_alpha) + 1

  # Upper tail probabilities for group 1: surv[k + 1] = P(X1 >= k)
  surv_h1 <- c(rev(cumsum(rev(dbinom(0:n1, n1, pi1_NRI)))), 0)
  surv_h0 <- c(rev(cumsum(rev(dbinom(0:n1, n1, pi_bar)))), 0)

  power <- sum(dbinom(0:n0, n0, pi0_NRI) * surv_h1[idx])
  type1_error <- sum(dbinom(0:n0, n0, pi_bar) * surv_h0[idx])

  result <- data.frame(
    n1 = n1,
    n0 = n0,
    pi1 = pi1,
    pi0 = pi0,
    omega1 = omega1,
    omega0 = omega0,
    rho1 = rho1,
    rho0 = rho0,
    pi1_NRI = pi1_NRI,
    pi0_NRI = pi0_NRI,
    effect_full = pi1 - pi0,
    effect_NRI = effect_NRI,
    alpha = alpha,
    pi_null = pi_bar,
    power = power,
    type1_error = type1_error
  )

  class(result) <- c("power_nri", "data.frame")

  return(result)
}

#' Print Method for power_nri Objects
#'
#' @param x An object of class \code{power_nri}.
#' @param ... Additional arguments (not used).
#'
#' @return \code{x}, invisibly.
#'
#' @export
print.power_nri <- function(x, ...) {
  cat("Exact Power Evaluation (NRI)\n")
  cat("============================\n")
  cat("Design parameters:\n")
  cat(sprintf("           n1 = %d\n", x$n1))
  cat(sprintf("           n0 = %d\n", x$n0))
  cat(sprintf("          pi1 = %.3f\n", x$pi1))
  cat(sprintf("          pi0 = %.3f\n", x$pi0))
  cat(sprintf("       omega1 = %.3f\n", x$omega1))
  cat(sprintf("       omega0 = %.3f\n", x$omega0))
  cat(sprintf("         rho1 = %.3f\n", x$rho1))
  cat(sprintf("         rho0 = %.3f\n", x$rho0))
  cat(sprintf("        alpha = %.4f\n", x$alpha))

  cat("\nNRI response probabilities:\n")
  cat(sprintf("      pi1_NRI = %.4f\n", x$pi1_NRI))
  cat(sprintf("      pi0_NRI = %.4f\n", x$pi0_NRI))
  cat(sprintf("  effect_full = %.4f\n", x$effect_full))
  cat(sprintf("   effect_NRI = %.4f\n", x$effect_NRI))

  cat("\nEvaluation results:\n")
  cat(sprintf("      pi_null = %.4f\n", x$pi_null))
  cat(sprintf("        power = %.4f\n", x$power))
  cat(sprintf("  type1_error = %.4f\n", x$type1_error))

  invisible(x)
}

#' Plot Method for power_nri Objects
#'
#' Exact power and exact type I error rate as a function of one design input,
#' with the remaining inputs held at the values stored in the object. A dashed
#' vertical line marks the value the object itself was evaluated at.
#'
#' @param x An object of class \code{power_nri}.
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
#' res <- power_nri(n1 = 40, n0 = 40, pi1 = 0.6, pi0 = 0.4,
#'                  omega1 = 0.2, omega0 = 0.2)
#' plot(res, vary = "rho", n_points = 11)
#'
#' @export
plot.power_nri <- function(x, vary = c("rho", "omega"), n_points = 41, ...) {
  vary <- match.arg(vary)
  .plot_power(x, vary = vary, n_points = n_points, fun = power_nri,
              label = "NRI analysis")
}

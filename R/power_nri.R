#' Exact Power and Type I Error Rate for an NRI Analysis
#'
#' Evaluates the exact power and the exact type I error rate of a one-sided
#' pooled-variance Z-test applied to a non-responder imputation (NRI) analysis,
#' in which every dropout is counted as a non-responder.
#'
#' @param n1 Integer. Randomized sample size for group 1 (treatment), n1 > 0.
#' @param n2 Integer. Randomized sample size for group 2 (control), n2 > 0.
#' @param p1 Numeric. Latent response probability for group 1, 0 < p1 < 1.
#' @param p2 Numeric. Latent response probability for group 2, 0 < p2 < 1.
#' @param omega1 Numeric. Dropout probability for group 1, 0 <= omega1 < 1.
#' @param omega2 Numeric. Dropout probability for group 2, 0 <= omega2 < 1.
#' @param rho1 Numeric. Correlation between the latent response and the dropout
#'   indicator in group 1 (default = 0).
#' @param rho2 Numeric. Correlation between the latent response and the dropout
#'   indicator in group 2 (default = 0).
#' @param alpha Numeric. One-sided significance level, 0 < alpha < 1
#'   (default = 0.025).
#' @param p_null Numeric or NULL. Common response probability under the null
#'   hypothesis, used only for the type I error rate. The default NULL uses the
#'   sample-size weighted average of the two NRI response probabilities. Supply
#'   a value to evaluate the size of the test at a different nuisance parameter.
#'
#' @return A \code{data.frame} of class \code{power_nri} with one row containing
#' the inputs together with:
#' \describe{
#'   \item{p1_NRI}{Observed response probability in group 1 under NRI}
#'   \item{p2_NRI}{Observed response probability in group 2 under NRI}
#'   \item{effect_latent}{Latent treatment effect, p1 - p2}
#'   \item{effect_NRI}{NRI treatment effect, p1_NRI - p2_NRI}
#'   \item{power}{Exact power}
#'   \item{type1_error}{Exact type I error rate}
#' }
#'
#' @details
#' Under NRI a patient counts as a responder only when the latent response
#' equals one and the patient completes the trial, so the number of observed
#' responders in group j follows a binomial distribution with size nj and
#' success probability
#'
#'   pj_NRI = P(R = 1, D = 0)
#'          = pj (1 - omegaj) - rhoj sqrt(pj (1 - pj) omegaj (1 - omegaj)).
#'
#' The randomized sample size is fixed, so no completer-count stage is involved
#' and the evaluation is a single enumeration over the observed responder pair.
#' As in \code{\link{power_cc}}, the rejection region is an upper set in the
#' group 1 count for every fixed group 2 count, so the enumeration reduces to
#' upper binomial tail probabilities evaluated at one threshold per group 2
#' count. This is exact; nothing is discarded.
#'
#' The type I error rate is evaluated at a common response probability, which by
#' default is the sample-size weighted average
#' p0 = (n1 p1_NRI + n2 p2_NRI) / (n1 + n2). That choice is a convention rather
#' than a property of the test: the size of a two-sample binomial test depends on
#' the nuisance parameter, and \code{p_null} exists so that the dependence can be
#' examined. The returned \code{p_null} column records the value used.
#'
#' The NRI estimand differs from the complete case estimand. Writing
#' sj = sqrt(pj (1 - pj) omegaj (1 - omegaj)),
#'
#'   delta_NRI = p1 (1 - omega1) - p2 (1 - omega2) - rho1 s1 + rho2 s2,
#'
#' which reduces to (1 - omega) delta_CC when the two groups share a dropout
#' probability and the latent response is independent of dropout within each
#' group. The two estimands can differ in sign as well as in magnitude, so a
#' non-null NRI contrast under equal latent response probabilities is not an
#' inflated type I error rate but a correctly sized test of a different null
#' hypothesis.
#'
#' @examples
#' power_nri(n1 = 60, n2 = 60, p1 = 0.6, p2 = 0.4,
#'           omega1 = 0.2, omega2 = 0.2, rho1 = 0, rho2 = 0)
#'
#' # Unequal dropout probabilities move the NRI estimand away from zero
#' # even when the latent response probabilities are equal
#' power_nri(n1 = 60, n2 = 60, p1 = 0.4, p2 = 0.4,
#'           omega1 = 0.05, omega2 = 0.20)
#'
#' @seealso \code{\link{power_cc}}, \code{\link{rho_bounds}}
#'
#' @importFrom stats dbinom qnorm
#' @export
power_nri <- function(n1, n2, p1, p2, omega1, omega2,
                      rho1 = 0, rho2 = 0, alpha = 0.025, p_null = NULL) {

  # Input validation
  if (n1 <= 0 || n2 <= 0) stop("Sample sizes must be positive")
  if (n1 != round(n1) || n2 != round(n2)) stop("Sample sizes must be integers")
  if (p1 <= 0 || p1 >= 1) stop("p1 must be between 0 and 1")
  if (p2 <= 0 || p2 >= 1) stop("p2 must be between 0 and 1")
  if (omega1 < 0 || omega1 >= 1) stop("omega1 must be between 0 and 1")
  if (omega2 < 0 || omega2 >= 1) stop("omega2 must be between 0 and 1")
  if (alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1")
  if (!is.null(p_null) && (length(p_null) != 1 || is.na(p_null) ||
                           p_null <= 0 || p_null >= 1)) {
    stop("p_null must be a single value between 0 and 1, or NULL")
  }

  n1 <- as.integer(round(n1))
  n2 <- as.integer(round(n2))

  # Validate correlation bounds
  ctol <- 1e-6
  if (omega1 > 0) {
    bounds1 <- rho_bounds(p = p1, omega = omega1)
    if (rho1 < bounds1$rho_lower - ctol || rho1 > bounds1$rho_upper + ctol) {
      stop(paste0("rho1 must be between ", round(bounds1$rho_lower, 4),
                  " and ", round(bounds1$rho_upper, 4)))
    }
  }
  if (omega2 > 0) {
    bounds2 <- rho_bounds(p = p2, omega = omega2)
    if (rho2 < bounds2$rho_lower - ctol || rho2 > bounds2$rho_upper + ctol) {
      stop(paste0("rho2 must be between ", round(bounds2$rho_lower, 4),
                  " and ", round(bounds2$rho_upper, 4)))
    }
  }

  # Observed response probabilities under NRI
  p1_NRI <- .joint_cells(p = p1, omega = omega1, rho = rho1)[["pi_10"]]
  p2_NRI <- .joint_cells(p = p2, omega = omega2, rho = rho2)[["pi_10"]]

  effect_NRI <- p1_NRI - p2_NRI

  # Null value: pooled by default, otherwise as supplied
  p0 <- if (is.null(p_null)) {
    (n1 * p1_NRI + n2 * p2_NRI) / (n1 + n2)
  } else {
    p_null
  }

  z_alpha <- qnorm(1 - alpha)

  # Rejection threshold in x1 for every x2. The randomized sample sizes are
  # fixed, so a single set of thresholds serves both hypotheses
  idx <- .reject_threshold(n1, n2, z_alpha) + 1

  # Upper tail probabilities for group 1: surv[k + 1] = P(X1 >= k)
  surv_h1 <- c(rev(cumsum(rev(dbinom(0:n1, n1, p1_NRI)))), 0)
  surv_h0 <- c(rev(cumsum(rev(dbinom(0:n1, n1, p0)))), 0)

  power <- sum(dbinom(0:n2, n2, p2_NRI) * surv_h1[idx])
  type1_error <- sum(dbinom(0:n2, n2, p0) * surv_h0[idx])

  result <- data.frame(
    n1 = n1,
    n2 = n2,
    p1 = p1,
    p2 = p2,
    omega1 = omega1,
    omega2 = omega2,
    rho1 = rho1,
    rho2 = rho2,
    p1_NRI = p1_NRI,
    p2_NRI = p2_NRI,
    effect_latent = p1 - p2,
    effect_NRI = effect_NRI,
    alpha = alpha,
    p_null = p0,
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
  cat(sprintf("             n1 = %d\n", x$n1))
  cat(sprintf("             n2 = %d\n", x$n2))
  cat(sprintf("             p1 = %.3f\n", x$p1))
  cat(sprintf("             p2 = %.3f\n", x$p2))
  cat(sprintf("         omega1 = %.3f\n", x$omega1))
  cat(sprintf("         omega2 = %.3f\n", x$omega2))
  cat(sprintf("           rho1 = %.3f\n", x$rho1))
  cat(sprintf("           rho2 = %.3f\n", x$rho2))
  cat(sprintf("          alpha = %.4f\n", x$alpha))

  cat("\nNRI response probabilities:\n")
  cat(sprintf("         p1_NRI = %.4f\n", x$p1_NRI))
  cat(sprintf("         p2_NRI = %.4f\n", x$p2_NRI))
  cat(sprintf("  effect_latent = %.4f\n", x$effect_latent))
  cat(sprintf("     effect_NRI = %.4f\n", x$effect_NRI))

  cat("\nEvaluation results:\n")
  cat(sprintf("         p_null = %.4f\n", x$p_null))
  cat(sprintf("          power = %.4f\n", x$power))
  cat(sprintf("    type1_error = %.4f\n", x$type1_error))

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
#' res <- power_nri(n1 = 40, n2 = 40, p1 = 0.6, p2 = 0.4,
#'                  omega1 = 0.2, omega2 = 0.2)
#' plot(res, vary = "rho", n_points = 11)
#'
#' @export
plot.power_nri <- function(x, vary = c("rho", "omega"), n_points = 41, ...) {
  vary <- match.arg(vary)
  .plot_power(x, vary = vary, n_points = n_points, fun = power_nri,
              label = "NRI analysis")
}

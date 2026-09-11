#' Sample Size under the Simple Inflation Method
#'
#' Computes the required sample size for a two-group superiority trial with a
#' binary endpoint under the simple inflation method: the sample size for a
#' design without dropout, divided by the probability of remaining in the trial.
#' This is the calculation that a complete case (CC) analysis under within-group
#' independence between response and dropout justifies.
#'
#' @param pi1 Numeric. Latent response probability for group 1, 0 < pi1 < 1.
#' @param pi0 Numeric. Latent response probability for group 0, 0 < pi0 < 1.
#' @param omega1 Numeric. Dropout probability for group 1, 0 <= omega1 < 1
#'   (default = 0).
#' @param omega0 Numeric. Dropout probability for group 0, 0 <= omega0 < 1
#'   (default = 0).
#' @param r Numeric. Allocation ratio n1 / n0, r > 0 (default = 1).
#' @param alpha Numeric. One-sided significance level (default = 0.025).
#' @param target_power Numeric. Target power (default = 0.8).
#'
#' @return A \code{data.frame} of class \code{sample_size_cc} with one row
#' containing the inputs together with:
#' \describe{
#'   \item{n1, n0, n_total}{Sample size for the design without dropout}
#'   \item{n1_adjust, n0_adjust, n_total_adjust}{Sample size after the
#'     n_j / (1 - omega_j) inflation}
#' }
#'
#' @details
#' Setting both dropout probabilities to zero returns the sample size for the
#' design without dropout, which is the starting point of the inflation.
#'
#' The inflation step assumes that the number of completers is what the analysis
#' uses and that the conditional response probability among completers equals
#' the latent response probability. That equality holds when the latent response
#' and dropout are independent within each group. It does not hold under NRI,
#' where a dropout is counted as a non-responder, and it is that mismatch that
#' \code{\link{sample_size_nri}} addresses.
#'
#' @examples
#' # Design without dropout
#' sample_size_cc(pi1 = 0.6, pi0 = 0.4)
#'
#' # With 20 percent dropout in both groups
#' sample_size_cc(pi1 = 0.6, pi0 = 0.4, omega1 = 0.2, omega0 = 0.2)
#'
#' @seealso \code{\link{sample_size_nri}}, \code{\link{power_cc}}
#'
#' @importFrom stats qnorm
#' @export
sample_size_cc <- function(pi1, pi0, omega1 = 0, omega0 = 0,
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

  # Step 1: sample size assuming no dropout
  r_pooled <- (r * pi1 + pi0) / (1 + r)

  # Variance under H0 (pooled) and H1 (separate)
  var_h0 <- r_pooled * (1 - r_pooled)
  var_h1 <- (pi1 * (1 - pi1) / r + pi0 * (1 - pi0)) / (1 + 1 / r)

  n0 <- (1 + 1 / r) / ((pi1 - pi0) ^ 2) *
    (qnorm(1 - alpha) * sqrt(var_h0) + qnorm(target_power) * sqrt(var_h1)) ^ 2

  n0 <- ceiling(n0)
  n1 <- ceiling(r * n0)
  n_total <- n1 + n0

  # Step 2: apply the n / (1 - omega) inflation
  n1_adjust <- if (omega1 > 0) ceiling(n1 / (1 - omega1)) else n1
  n0_adjust <- if (omega0 > 0) ceiling(n0 / (1 - omega0)) else n0
  n_total_adjust <- n1_adjust + n0_adjust

  # Under within-group independence the conditional response probabilities
  # among completers equal the latent response probabilities
  pi1_CC <- pi1
  pi0_CC <- pi0

  result <- data.frame(
    pi1 = pi1,
    pi0 = pi0,
    omega1 = omega1,
    omega0 = omega0,
    r = r,
    alpha = alpha,
    target_power = target_power,
    pi1_CC = pi1_CC,
    pi0_CC = pi0_CC,
    effect_full = pi1 - pi0,
    effect_CC = pi1_CC - pi0_CC,
    n1 = n1,
    n0 = n0,
    n_total = n_total,
    n1_adjust = n1_adjust,
    n0_adjust = n0_adjust,
    n_total_adjust = n_total_adjust
  )

  class(result) <- c("sample_size_cc", "data.frame")

  return(result)
}

#' Print Method for sample_size_cc Objects
#'
#' @param x An object of class \code{sample_size_cc}.
#' @param ... Additional arguments (not used).
#'
#' @return \code{x}, invisibly.
#'
#' @export
print.sample_size_cc <- function(x, ...) {
  cat("Sample Size under the Simple Inflation Method\n")
  cat("=============================================\n")
  cat("Design parameters:\n")
  cat(sprintf("             pi1 = %.3f\n", x$pi1))
  cat(sprintf("             pi0 = %.3f\n", x$pi0))
  cat(sprintf("          omega1 = %.3f\n", x$omega1))
  cat(sprintf("          omega0 = %.3f\n", x$omega0))
  cat(sprintf("               r = %.4g\n", x$r))
  cat(sprintf("           alpha = %.4f\n", x$alpha))
  cat(sprintf("    target_power = %.2f\n", x$target_power))

  cat("\nComplete case probabilities (within-group independence):\n")
  cat(sprintf("          pi1_CC = %.4f\n", x$pi1_CC))
  cat(sprintf("          pi0_CC = %.4f\n", x$pi0_CC))
  cat(sprintf("     effect_full = %.4f\n", x$effect_full))
  cat(sprintf("       effect_CC = %.4f\n", x$effect_CC))

  cat("\nRequired sample size (no dropout):\n")
  cat(sprintf("              n1 = %d\n", x$n1))
  cat(sprintf("              n0 = %d\n", x$n0))
  cat(sprintf("         n_total = %d\n", x$n_total))

  cat("\nRequired sample size (inflated by n / (1 - omega)):\n")
  cat(sprintf("       n1_adjust = %d\n", x$n1_adjust))
  cat(sprintf("       n0_adjust = %d\n", x$n0_adjust))
  cat(sprintf("  n_total_adjust = %d\n", x$n_total_adjust))

  invisible(x)
}

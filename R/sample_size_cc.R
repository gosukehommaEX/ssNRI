#' Sample Size under the Simple Inflation Method
#'
#' Computes the required sample size for a two-group superiority trial with a
#' binary endpoint under the simple inflation method: the sample size for a
#' design without dropout, divided by the probability of remaining in the trial.
#' This is the calculation that a complete case (CC) analysis under within-group
#' independence between response and dropout justifies.
#'
#' @param p1 Numeric. Latent response probability for group 1, 0 < p1 < 1.
#' @param p2 Numeric. Latent response probability for group 2, 0 < p2 < 1.
#' @param omega1 Numeric. Dropout probability for group 1, 0 <= omega1 < 1
#'   (default = 0).
#' @param omega2 Numeric. Dropout probability for group 2, 0 <= omega2 < 1
#'   (default = 0).
#' @param r Numeric. Allocation ratio n1 / n2, r > 0 (default = 1).
#' @param alpha Numeric. One-sided significance level (default = 0.025).
#' @param target_power Numeric. Target power (default = 0.8).
#'
#' @return A \code{data.frame} of class \code{sample_size_cc} with one row
#' containing the inputs together with:
#' \describe{
#'   \item{n1, n2, n_total}{Sample size for the design without dropout}
#'   \item{n1_adjust, n2_adjust, n_total_adjust}{Sample size after the
#'     nj / (1 - omegaj) inflation}
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
#' sample_size_cc(p1 = 0.6, p2 = 0.4)
#'
#' # With 20 percent dropout in both groups
#' sample_size_cc(p1 = 0.6, p2 = 0.4, omega1 = 0.2, omega2 = 0.2)
#'
#' @seealso \code{\link{sample_size_nri}}, \code{\link{power_cc}}
#'
#' @importFrom stats qnorm
#' @export
sample_size_cc <- function(p1, p2, omega1 = 0, omega2 = 0,
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

  # Step 1: sample size assuming no dropout
  r_pooled <- (r * p1 + p2) / (1 + r)

  # Variance under H0 (pooled) and H1 (separate)
  var_h0 <- r_pooled * (1 - r_pooled)
  var_h1 <- (p1 * (1 - p1) / r + p2 * (1 - p2)) / (1 + 1 / r)

  n2 <- (1 + 1 / r) / ((p1 - p2) ^ 2) *
    (qnorm(1 - alpha) * sqrt(var_h0) + qnorm(target_power) * sqrt(var_h1)) ^ 2

  n2 <- ceiling(n2)
  n1 <- ceiling(r * n2)
  n_total <- n1 + n2

  # Step 2: apply the n / (1 - omega) inflation
  n1_adjust <- if (omega1 > 0) ceiling(n1 / (1 - omega1)) else n1
  n2_adjust <- if (omega2 > 0) ceiling(n2 / (1 - omega2)) else n2
  n_total_adjust <- n1_adjust + n2_adjust

  # Under within-group independence the conditional response probabilities
  # among completers equal the latent response probabilities
  p1_CC <- p1
  p2_CC <- p2

  result <- data.frame(
    p1 = p1,
    p2 = p2,
    omega1 = omega1,
    omega2 = omega2,
    r = r,
    alpha = alpha,
    target_power = target_power,
    p1_CC = p1_CC,
    p2_CC = p2_CC,
    effect_latent = p1 - p2,
    effect_CC = p1_CC - p2_CC,
    n1 = n1,
    n2 = n2,
    n_total = n_total,
    n1_adjust = n1_adjust,
    n2_adjust = n2_adjust,
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
  cat(sprintf("             p1 = %.3f\n", x$p1))
  cat(sprintf("             p2 = %.3f\n", x$p2))
  cat(sprintf("         omega1 = %.3f\n", x$omega1))
  cat(sprintf("         omega2 = %.3f\n", x$omega2))
  cat(sprintf("              r = %.4g\n", x$r))
  cat(sprintf("          alpha = %.4f\n", x$alpha))
  cat(sprintf("   target_power = %.2f\n", x$target_power))

  cat("\nComplete case probabilities (within-group independence):\n")
  cat(sprintf("          p1_CC = %.4f\n", x$p1_CC))
  cat(sprintf("          p2_CC = %.4f\n", x$p2_CC))
  cat(sprintf("  effect_latent = %.4f\n", x$effect_latent))
  cat(sprintf("      effect_CC = %.4f\n", x$effect_CC))

  cat("\nRequired sample size (no dropout):\n")
  cat(sprintf("             n1 = %d\n", x$n1))
  cat(sprintf("             n2 = %d\n", x$n2))
  cat(sprintf("        n_total = %d\n", x$n_total))

  cat("\nRequired sample size (inflated by n / (1 - omega)):\n")
  cat(sprintf("      n1_adjust = %d\n", x$n1_adjust))
  cat(sprintf("      n2_adjust = %d\n", x$n2_adjust))
  cat(sprintf(" n_total_adjust = %d\n", x$n_total_adjust))

  invisible(x)
}

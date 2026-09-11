#' Exact Power and Type I Error Rate for a Complete Case Analysis
#'
#' Evaluates the exact power and the exact type I error rate of a one-sided
#' pooled-variance Z-test applied to a complete case (CC) analysis, accounting
#' for the randomness of both the number of completers and the number of
#' responders among completers.
#'
#' @param n1 Integer. Randomized sample size for group 1 (treatment), n1 > 0.
#' @param n0 Integer. Randomized sample size for group 0 (control), n0 > 0.
#' @param pi1 Numeric. Latent response probability for group 1, 0 < pi1 < 1.
#' @param pi0 Numeric. Latent response probability for group 0, 0 < pi0 < 1.
#' @param omega1 Numeric. Dropout probability for group 1, 0 <= omega1 < 1.
#' @param omega0 Numeric. Dropout probability for group 0, 0 <= omega0 < 1.
#' @param rho1 Numeric. Correlation between the latent response and the dropout
#'   indicator in group 1 (default = 0). Must lie within the feasible range
#'   returned by \code{rho_bounds}.
#' @param rho0 Numeric. Correlation between the latent response and the dropout
#'   indicator in group 0 (default = 0).
#' @param alpha Numeric. One-sided significance level, 0 < alpha < 1
#'   (default = 0.025).
#' @param pi_null Numeric or NULL. Common conditional response probability among
#'   completers under the null hypothesis, used only for the type I error rate.
#'   The default NULL uses the sample-size weighted average of pi1_CC and pi0_CC.
#'   Supply a value to evaluate the size of the test at a different nuisance
#'   parameter. The pooled dropout probability used for the completer stage is
#'   unaffected.
#' @param tol Numeric. Optional truncation threshold for the completer-count
#'   distribution, tol >= 0 (default = 0). With the default the calculation is
#'   exact, in the sense that every completer count with positive probability is
#'   enumerated. A positive value discards completer counts whose probability is
#'   at most tol under both hypotheses; see Details for the resulting error bound.
#'
#' @return A \code{data.frame} of class \code{power_cc} with one row containing
#' the inputs together with:
#' \describe{
#'   \item{pi1_CC}{Conditional response probability among completers in group 1}
#'   \item{pi0_CC}{Conditional response probability among completers in group 0}
#'   \item{effect_full}{Full-data treatment effect, pi1 - pi0}
#'   \item{effect_CC}{Complete case treatment effect, pi1_CC - pi0_CC}
#'   \item{power}{Exact power}
#'   \item{type1_error}{Exact type I error rate}
#' }
#'
#' @details
#' Let M1 and M0 denote the numbers of completers and X1 and X0 the numbers of
#' responders among completers. Under the bivariate Bernoulli model,
#' M_j follows a binomial distribution with size n_j and success probability
#' 1 - omega_j, and, conditionally on M_j = m_j, X_j follows a binomial distribution
#' with size m_j and success probability
#'
#'   pi_j_CC = P(R_ij = 1 | D_ij = 0)
#'         = (pi_j (1 - omega_j) - rho_j sqrt(pi_j (1 - pi_j) omega_j (1 - omega_j)))
#'           / (1 - omega_j).
#'
#' Under within-group independence (rho_j = 0) this reduces to pi_j_CC = pi_j, which
#' is what justifies the simple inflation method implemented in
#' \code{sample_size_cc}.
#'
#' Exact evaluation therefore requires the sum
#'
#'   sum over (m1, m0) of P(M1 = m1) P(M0 = m0) times the conditional rejection
#'   probability given (m1, m0),
#'
#' with completer counts of zero excluded because the test statistic is then
#' undefined.
#'
#' Rather than enumerating every responder pair (x1, x0) within each (m1, m0)
#' cell, this implementation uses the following property. For fixed m1, m0 and
#' x0, the pooled-variance Z statistic is strictly increasing in x1 wherever it
#' is positive, so the rejection region is an upper set in x1 and is described by
#' a single threshold c(m1, m0, x0). The conditional rejection probability is
#' then obtained from upper binomial tail probabilities in one pass over x0. The
#' threshold is located as the larger root of the quadratic inequality equivalent
#' to Z > z(1 - alpha) and is then verified against the defining inequality, so
#' the thresholds are exact rather than approximate. This is a regrouping of the
#' same finite sum: the returned values agree with a full double enumeration to
#' within floating point summation error.
#'
#' Setting \code{tol} above zero additionally discards completer counts of small
#' probability. If the discarded completer-count probability mass is eps, the
#' reported power and type I error rate are lower bounds for the exact values and
#' differ from them by at most eps. The default \code{tol = 0} discards nothing.
#'
#' Under the null hypothesis a common conditional response probability
#'
#'   pi_bar_CC = (n1 pi1_CC + n0 pi0_CC) / (n1 + n0)
#'
#' and a pooled dropout probability
#'
#'   omega_bar = (n1 omega1 + n0 omega0) / (n1 + n0)
#'
#' are used. The first of these is a convention rather than a property of the
#' test: the size of a two-sample binomial test depends on the nuisance
#' parameter, and \code{pi_null} exists so that the dependence can be examined.
#' The returned \code{pi_null} column records the value used.
#'
#' The rejection region depends only on (x1, x0, m1, m0) and is therefore
#' shared by the two hypotheses, so power and type I error rate are accumulated
#' in a single pass.
#'
#' @examples
#' # Within-group independence between response and dropout
#' power_cc(n1 = 60, n0 = 60, pi1 = 0.6, pi0 = 0.4,
#'          omega1 = 0.2, omega0 = 0.2, rho1 = 0, rho0 = 0, alpha = 0.025)
#'
#' # Positive correlation: responders are more likely to drop out
#' power_cc(n1 = 60, n0 = 60, pi1 = 0.6, pi0 = 0.4,
#'          omega1 = 0.2, omega0 = 0.2, rho1 = 0.3, rho0 = 0.3)
#'
#' @seealso \code{\link{rho_bounds}}
#'
#' @importFrom stats dbinom qnorm
#' @export
power_cc <- function(n1, n0, pi1, pi0, omega1, omega0,
                     rho1 = 0, rho0 = 0, alpha = 0.025,
                     pi_null = NULL, tol = 0) {

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
  if (tol < 0) stop("tol must be non-negative")

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

  # Bivariate Bernoulli joint cell probabilities
  cells1 <- .joint_cells(pi = pi1, omega = omega1, rho = rho1)
  cells0 <- .joint_cells(pi = pi0, omega = omega0, rho = rho0)

  # Conditional response probabilities among completers
  pi1_CC <- cells1[["pi_10"]] / (cells1[["pi_00"]] + cells1[["pi_10"]])
  pi0_CC <- cells0[["pi_10"]] / (cells0[["pi_00"]] + cells0[["pi_10"]])

  # Null values: the response probability is pooled by default, otherwise as
  # supplied; the dropout probability is always pooled
  pi_bar_CC <- if (is.null(pi_null)) {
    (n1 * pi1_CC + n0 * pi0_CC) / (n1 + n0)
  } else {
    pi_null
  }
  omega_bar <- (n1 * omega1 + n0 * omega0) / (n1 + n0)

  z_alpha <- qnorm(1 - alpha)

  # Completer-count distributions
  prob_m1_h1 <- dbinom(0:n1, n1, 1 - omega1)
  prob_m0_h1 <- dbinom(0:n0, n0, 1 - omega0)
  prob_m1_h0 <- dbinom(0:n1, n1, 1 - omega_bar)
  prob_m0_h0 <- dbinom(0:n0, n0, 1 - omega_bar)

  # Completer counts to enumerate; m = 0 is always excluded because the test
  # statistic is undefined when a group contributes no completers
  # Kept as doubles throughout: the quadratic coefficients below involve
  # products of five counts, which overflows R's 32-bit integer type from
  # m1 = m0 = 64 upwards
  keep_m1 <- as.numeric(which(pmax(prob_m1_h1, prob_m1_h0) > tol) - 1L)
  keep_m0 <- as.numeric(which(pmax(prob_m0_h1, prob_m0_h0) > tol) - 1L)
  keep_m1 <- keep_m1[keep_m1 >= 1]
  keep_m0 <- keep_m0[keep_m0 >= 1]

  if (length(keep_m1) == 0L || length(keep_m0) == 0L) {
    stop("No completer counts left to enumerate; check omega1, omega0 and tol")
  }

  # Responder distributions in group 0 depend on m0 only, so cache them
  f0_h1 <- vector("list", n0 + 1L)
  f0_h0 <- vector("list", n0 + 1L)
  for (m0 in keep_m0) {
    f0_h1[[m0 + 1]] <- dbinom(0:m0, m0, pi0_CC)
    f0_h0[[m0 + 1]] <- dbinom(0:m0, m0, pi_bar_CC)
  }

  power <- 0
  type1_error <- 0

  for (m1 in keep_m1) {

    # Upper tail probabilities for group 1: surv[k + 1] = P(X1 >= k),
    # k = 0, ..., m1 + 1, so that surv[m1 + 2] = 0 marks "never rejects"
    surv_h1 <- c(rev(cumsum(rev(dbinom(0:m1, m1, pi1_CC)))), 0)
    surv_h0 <- c(rev(cumsum(rev(dbinom(0:m1, m1, pi_bar_CC)))), 0)

    for (m0 in keep_m0) {

      idx <- .reject_threshold(m1, m0, z_alpha) + 1
      power <- power +
        prob_m1_h1[m1 + 1] * prob_m0_h1[m0 + 1] *
        sum(f0_h1[[m0 + 1]] * surv_h1[idx])
      type1_error <- type1_error +
        prob_m1_h0[m1 + 1] * prob_m0_h0[m0 + 1] *
        sum(f0_h0[[m0 + 1]] * surv_h0[idx])
    }
  }

  result <- data.frame(
    n1 = n1,
    n0 = n0,
    pi1 = pi1,
    pi0 = pi0,
    omega1 = omega1,
    omega0 = omega0,
    rho1 = rho1,
    rho0 = rho0,
    pi1_CC = pi1_CC,
    pi0_CC = pi0_CC,
    effect_full = pi1 - pi0,
    effect_CC = pi1_CC - pi0_CC,
    alpha = alpha,
    pi_null = pi_bar_CC,
    power = power,
    type1_error = type1_error
  )

  class(result) <- c("power_cc", "data.frame")

  return(result)
}

#' Print Method for power_cc Objects
#'
#' @param x An object of class \code{power_cc}.
#' @param ... Additional arguments (not used).
#'
#' @return \code{x}, invisibly.
#'
#' @export
print.power_cc <- function(x, ...) {
  cat("Exact Power Evaluation (Complete Case)\n")
  cat("=======================================\n")
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

  cat("\nComplete case probabilities:\n")
  cat(sprintf("       pi1_CC = %.4f\n", x$pi1_CC))
  cat(sprintf("       pi0_CC = %.4f\n", x$pi0_CC))
  cat(sprintf("  effect_full = %.4f\n", x$effect_full))
  cat(sprintf("    effect_CC = %.4f\n", x$effect_CC))

  cat("\nEvaluation results:\n")
  cat(sprintf("      pi_null = %.4f\n", x$pi_null))
  cat(sprintf("        power = %.4f\n", x$power))
  cat(sprintf("  type1_error = %.4f\n", x$type1_error))

  invisible(x)
}

#' Plot Method for power_cc Objects
#'
#' Exact power and exact type I error rate as a function of one design input,
#' with the remaining inputs held at the values stored in the object. A dashed
#' vertical line marks the value the object itself was evaluated at.
#'
#' @param x An object of class \code{power_cc}.
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
#' res <- power_cc(n1 = 40, n0 = 40, pi1 = 0.6, pi0 = 0.4,
#'                 omega1 = 0.2, omega0 = 0.2)
#' plot(res, vary = "rho", n_points = 11)
#'
#' @importFrom rlang .data
#' @export
plot.power_cc <- function(x, vary = c("rho", "omega"), n_points = 41, ...) {
  vary <- match.arg(vary)
  .plot_power(x, vary = vary, n_points = n_points, fun = power_cc,
              label = "Complete case analysis")
}

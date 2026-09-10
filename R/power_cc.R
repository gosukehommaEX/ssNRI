#' Exact Power and Type I Error Rate for a Complete Case Analysis
#'
#' Evaluates the exact power and the exact type I error rate of a one-sided
#' pooled-variance Z-test applied to a complete case (CC) analysis, accounting
#' for the randomness of both the number of completers and the number of
#' responders among completers.
#'
#' @param n1 Integer. Randomized sample size for group 1 (treatment), n1 > 0.
#' @param n2 Integer. Randomized sample size for group 2 (control), n2 > 0.
#' @param p1 Numeric. Latent response probability for group 1, 0 < p1 < 1.
#' @param p2 Numeric. Latent response probability for group 2, 0 < p2 < 1.
#' @param omega1 Numeric. Dropout probability for group 1, 0 <= omega1 < 1.
#' @param omega2 Numeric. Dropout probability for group 2, 0 <= omega2 < 1.
#' @param rho1 Numeric. Correlation between the latent response and the dropout
#'   indicator in group 1 (default = 0). Must lie within the feasible range
#'   returned by \code{rho_bounds}.
#' @param rho2 Numeric. Correlation between the latent response and the dropout
#'   indicator in group 2 (default = 0).
#' @param alpha Numeric. One-sided significance level, 0 < alpha < 1
#'   (default = 0.025).
#' @param p_null Numeric or NULL. Common conditional response probability among
#'   completers under the null hypothesis, used only for the type I error rate.
#'   The default NULL uses the sample-size weighted average of p1_CC and p2_CC.
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
#'   \item{p1_CC}{Conditional response probability among completers in group 1}
#'   \item{p2_CC}{Conditional response probability among completers in group 2}
#'   \item{effect_latent}{Latent treatment effect, p1 - p2}
#'   \item{effect_CC}{Complete case treatment effect, p1_CC - p2_CC}
#'   \item{power}{Exact power}
#'   \item{type1_error}{Exact type I error rate}
#' }
#'
#' @details
#' Let M1 and M2 denote the numbers of completers and X1 and X2 the numbers of
#' responders among completers. Under the bivariate Bernoulli model,
#' Mj follows a binomial distribution with size nj and success probability
#' 1 - omegaj, and, conditionally on Mj = mj, Xj follows a binomial distribution
#' with size mj and success probability
#'
#'   pj_CC = P(R = 1 | D = 0)
#'         = (pj (1 - omegaj) - rhoj sqrt(pj (1 - pj) omegaj (1 - omegaj)))
#'           / (1 - omegaj).
#'
#' Under within-group independence (rhoj = 0) this reduces to pj_CC = pj, which
#' is what justifies the simple inflation method implemented in
#' \code{sample_size_cc}.
#'
#' Exact evaluation therefore requires the sum
#'
#'   sum over (m1, m2) of P(M1 = m1) P(M2 = m2) times the conditional rejection
#'   probability given (m1, m2),
#'
#' with completer counts of zero excluded because the test statistic is then
#' undefined.
#'
#' Rather than enumerating every responder pair (x1, x2) within each (m1, m2)
#' cell, this implementation uses the following property. For fixed m1, m2 and
#' x2, the pooled-variance Z statistic is strictly increasing in x1 wherever it
#' is positive, so the rejection region is an upper set in x1 and is described by
#' a single threshold c(m1, m2, x2). The conditional rejection probability is
#' then obtained from upper binomial tail probabilities in one pass over x2. The
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
#'   p0_CC = (n1 p1_CC + n2 p2_CC) / (n1 + n2)
#'
#' and a pooled dropout probability
#'
#'   d0 = (n1 omega1 + n2 omega2) / (n1 + n2)
#'
#' are used. The first of these is a convention rather than a property of the
#' test: the size of a two-sample binomial test depends on the nuisance
#' parameter, and \code{p_null} exists so that the dependence can be examined.
#' The returned \code{p_null} column records the value used.
#'
#' The rejection region depends only on (x1, x2, m1, m2) and is therefore
#' shared by the two hypotheses, so power and type I error rate are accumulated
#' in a single pass.
#'
#' @examples
#' # Within-group independence between response and dropout
#' power_cc(n1 = 60, n2 = 60, p1 = 0.6, p2 = 0.4,
#'          omega1 = 0.2, omega2 = 0.2, rho1 = 0, rho2 = 0, alpha = 0.025)
#'
#' # Positive correlation: responders are more likely to drop out
#' power_cc(n1 = 60, n2 = 60, p1 = 0.6, p2 = 0.4,
#'          omega1 = 0.2, omega2 = 0.2, rho1 = 0.3, rho2 = 0.3)
#'
#' @seealso \code{\link{rho_bounds}}
#'
#' @importFrom stats dbinom qnorm
#' @export
power_cc <- function(n1, n2, p1, p2, omega1, omega2,
                     rho1 = 0, rho2 = 0, alpha = 0.025,
                     p_null = NULL, tol = 0) {

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
  if (tol < 0) stop("tol must be non-negative")

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

  # Bivariate Bernoulli joint cell probabilities
  cells1 <- .joint_cells(p = p1, omega = omega1, rho = rho1)
  cells2 <- .joint_cells(p = p2, omega = omega2, rho = rho2)

  # Conditional response probabilities among completers
  p1_CC <- cells1[["pi_10"]] / (cells1[["pi_00"]] + cells1[["pi_10"]])
  p2_CC <- cells2[["pi_10"]] / (cells2[["pi_00"]] + cells2[["pi_10"]])

  # Null values: the response probability is pooled by default, otherwise as
  # supplied; the dropout probability is always pooled
  p0_CC <- if (is.null(p_null)) {
    (n1 * p1_CC + n2 * p2_CC) / (n1 + n2)
  } else {
    p_null
  }
  d0 <- (n1 * omega1 + n2 * omega2) / (n1 + n2)

  z_alpha <- qnorm(1 - alpha)

  # Completer-count distributions
  prob_m1_h1 <- dbinom(0:n1, n1, 1 - omega1)
  prob_m2_h1 <- dbinom(0:n2, n2, 1 - omega2)
  prob_m1_h0 <- dbinom(0:n1, n1, 1 - d0)
  prob_m2_h0 <- dbinom(0:n2, n2, 1 - d0)

  # Completer counts to enumerate; m = 0 is always excluded because the test
  # statistic is undefined when a group contributes no completers
  # Kept as doubles throughout: the quadratic coefficients below involve
  # products of five counts, which overflows R's 32-bit integer type from
  # m1 = m2 = 64 upwards
  keep_m1 <- as.numeric(which(pmax(prob_m1_h1, prob_m1_h0) > tol) - 1L)
  keep_m2 <- as.numeric(which(pmax(prob_m2_h1, prob_m2_h0) > tol) - 1L)
  keep_m1 <- keep_m1[keep_m1 >= 1]
  keep_m2 <- keep_m2[keep_m2 >= 1]

  if (length(keep_m1) == 0L || length(keep_m2) == 0L) {
    stop("No completer counts left to enumerate; check omega1, omega2 and tol")
  }

  # Responder distributions in group 2 depend on m2 only, so cache them
  f2_h1 <- vector("list", n2 + 1L)
  f2_h0 <- vector("list", n2 + 1L)
  for (m2 in keep_m2) {
    f2_h1[[m2 + 1]] <- dbinom(0:m2, m2, p2_CC)
    f2_h0[[m2 + 1]] <- dbinom(0:m2, m2, p0_CC)
  }

  power <- 0
  type1_error <- 0

  for (m1 in keep_m1) {

    # Upper tail probabilities for group 1: surv[k + 1] = P(X1 >= k),
    # k = 0, ..., m1 + 1, so that surv[m1 + 2] = 0 marks "never rejects"
    surv_h1 <- c(rev(cumsum(rev(dbinom(0:m1, m1, p1_CC)))), 0)
    surv_h0 <- c(rev(cumsum(rev(dbinom(0:m1, m1, p0_CC)))), 0)

    for (m2 in keep_m2) {

      idx <- .reject_threshold(m1, m2, z_alpha) + 1
      power <- power +
        prob_m1_h1[m1 + 1] * prob_m2_h1[m2 + 1] *
        sum(f2_h1[[m2 + 1]] * surv_h1[idx])
      type1_error <- type1_error +
        prob_m1_h0[m1 + 1] * prob_m2_h0[m2 + 1] *
        sum(f2_h0[[m2 + 1]] * surv_h0[idx])
    }
  }

  result <- data.frame(
    n1 = n1,
    n2 = n2,
    p1 = p1,
    p2 = p2,
    omega1 = omega1,
    omega2 = omega2,
    rho1 = rho1,
    rho2 = rho2,
    p1_CC = p1_CC,
    p2_CC = p2_CC,
    effect_latent = p1 - p2,
    effect_CC = p1_CC - p2_CC,
    alpha = alpha,
    p_null = p0_CC,
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
  cat(sprintf("             n1 = %d\n", x$n1))
  cat(sprintf("             n2 = %d\n", x$n2))
  cat(sprintf("             p1 = %.3f\n", x$p1))
  cat(sprintf("             p2 = %.3f\n", x$p2))
  cat(sprintf("         omega1 = %.3f\n", x$omega1))
  cat(sprintf("         omega2 = %.3f\n", x$omega2))
  cat(sprintf("           rho1 = %.3f\n", x$rho1))
  cat(sprintf("           rho2 = %.3f\n", x$rho2))
  cat(sprintf("          alpha = %.4f\n", x$alpha))

  cat("\nComplete case probabilities:\n")
  cat(sprintf("          p1_CC = %.4f\n", x$p1_CC))
  cat(sprintf("          p2_CC = %.4f\n", x$p2_CC))
  cat(sprintf("  effect_latent = %.4f\n", x$effect_latent))
  cat(sprintf("      effect_CC = %.4f\n", x$effect_CC))

  cat("\nEvaluation results:\n")
  cat(sprintf("         p_null = %.4f\n", x$p_null))
  cat(sprintf("          power = %.4f\n", x$power))
  cat(sprintf("    type1_error = %.4f\n", x$type1_error))

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
#' res <- power_cc(n1 = 40, n2 = 40, p1 = 0.6, p2 = 0.4,
#'                 omega1 = 0.2, omega2 = 0.2)
#' plot(res, vary = "rho", n_points = 11)
#'
#' @importFrom rlang .data
#' @export
plot.power_cc <- function(x, vary = c("rho", "omega"), n_points = 41, ...) {
  vary <- match.arg(vary)
  .plot_power(x, vary = vary, n_points = n_points, fun = power_cc,
              label = "Complete case analysis")
}

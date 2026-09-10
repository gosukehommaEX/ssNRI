#' Identified Range for the Latent Response Probability
#'
#' Given a response rate reported by a historical trial under a stated
#' missing-data handling rule, together with the number randomized and the
#' number of dropouts, returns the range of latent response probabilities
#' compatible with what was reported.
#'
#' @param p_obs Numeric. Response probability as reported by the historical
#'   trial, on the scale implied by \code{method}, 0 <= p_obs <= 1.
#' @param N_randomized Integer. Number of patients randomized in that group.
#' @param N_dropout Integer. Number of dropouts in that group,
#'   0 <= N_dropout <= N_randomized.
#' @param method Character. Missing-data handling rule the historical trial
#'   applied: \code{"nri"} (default), \code{"cc"}, \code{"bc"} (best case) or
#'   \code{"wc"} (worst case).
#' @param group Character. Which group the reported rate refers to,
#'   \code{"treatment"} (default) or \code{"control"}. Relevant only for the
#'   best and worst case rules, which are asymmetric across groups.
#'
#' @return A \code{data.frame} of class \code{latent_bounds} with one row
#' containing the inputs together with \code{omega}, \code{p_nri},
#' \code{p_lower}, \code{p_upper} and \code{p_midpoint}.
#'
#' @details
#' A historical trial reports a rate produced by an analysis rule. The rule is
#' a decision about how to count dropouts; it is not a statement about what the
#' dropouts' outcomes actually were. The latent response probability is
#' therefore not identified by the reported rate, and only a range can be
#' recovered.
#'
#' Every rule is first placed on the NRI scale, that is on the joint probability
#' of responding and completing. Counting dropouts as non-responders reports
#' that joint probability directly, so p_NRI = p_obs; this is what NRI does in
#' both groups, what a best case analysis does in the control group and what a
#' worst case analysis does in the treatment group. Counting dropouts as
#' responders adds the whole dropout probability, so p_NRI = p_obs - omega; this
#' is what a best case analysis does in the treatment group and a worst case
#' analysis in the control group. A complete case analysis reports the rate
#' conditional on completing, so p_NRI = p_obs (1 - omega). Each conversion is an
#' identity and requires no assumption about the missingness mechanism.
#'
#' The latent probability then lies between p_NRI, attained when no dropout
#' would have responded, and p_NRI + omega, attained when every dropout would
#' have responded:
#'
#'   p_lower = p_NRI,   p_upper = p_NRI + omega.
#'
#' The width of the interval is exactly omega, so a historical trial with little
#' dropout pins the latent probability down tightly and one with heavy dropout
#' barely constrains it. Setting a joint cell to zero, as is sometimes done to
#' read a latent probability off a best or worst case rate, selects one endpoint
#' of this interval rather than identifying a point: the value that a best case
#' treatment rate is often taken to give is p_upper, and it carries the
#' assumption that every dropout in that group would have responded.
#'
#' Planning should carry the interval rather than a single recovered value: pass
#' candidate values from this range to \code{\link{infer_rho}} to obtain the
#' corresponding correlation, and to \code{\link{sample_size_nri}} to obtain the
#' corresponding sample size.
#'
#' A bound outside the unit interval indicates that the reported rate and the
#' dropout count are mutually incompatible, usually because the stated rule was
#' not the one actually applied. A warning is issued and the bound is truncated.
#'
#' @examples
#' # A historical trial reporting 40 percent under NRI, with 15 percent dropout
#' latent_bounds(p_obs = 0.40, N_randomized = 200, N_dropout = 30)
#'
#' # The same trial reported on the complete case scale
#' latent_bounds(p_obs = 0.47, N_randomized = 200, N_dropout = 30,
#'               method = "cc")
#'
#' # A best case rate in the treatment group counts dropouts as responders
#' latent_bounds(p_obs = 0.55, N_randomized = 200, N_dropout = 30,
#'               method = "bc", group = "treatment")
#'
#' @seealso \code{\link{infer_rho}}, \code{\link{sample_size_nri}}
#'
#' @export
latent_bounds <- function(p_obs, N_randomized, N_dropout,
                          method = c("nri", "cc", "bc", "wc"),
                          group  = c("treatment", "control")) {

  method <- match.arg(method)
  group  <- match.arg(group)

  # Input validation
  if (p_obs < 0 || p_obs > 1) stop("p_obs must be between 0 and 1")
  if (N_randomized <= 0 || N_randomized != round(N_randomized)) {
    stop("N_randomized must be a positive integer")
  }
  if (N_dropout < 0 || N_dropout > N_randomized ||
      N_dropout != round(N_dropout)) {
    stop("N_dropout must be a non-negative integer not exceeding N_randomized")
  }

  omega <- N_dropout / N_randomized

  # Special case: no dropout, so the latent probability is identified
  if (omega == 0) {
    message("No dropout observed. The latent p equals p_obs with certainty.")
    result <- data.frame(
      method       = method,
      group        = group,
      p_obs        = p_obs,
      N_randomized = N_randomized,
      N_dropout    = 0L,
      omega        = 0,
      p_nri        = p_obs,
      p_lower      = p_obs,
      p_upper      = p_obs,
      p_midpoint   = p_obs,
      stringsAsFactors = FALSE
    )
    class(result) <- c("latent_bounds", "data.frame")
    return(result)
  }

  # How the stated rule counted the dropouts
  counts_dropouts_as_responders <-
    (method == "bc" && group == "treatment") ||
    (method == "wc" && group == "control")

  # Place the reported rate on the NRI scale; each conversion is an identity
  p_nri <- if (method == "cc") {
    p_obs * (1 - omega)
  } else if (counts_dropouts_as_responders) {
    p_obs - omega
  } else {
    p_obs
  }

  # No dropout would have responded, versus every dropout would have
  p_lower <- p_nri
  p_upper <- p_nri + omega

  if (p_lower < 0) {
    warning(paste0(
      "Lower bound p_lower = ", round(p_lower, 3), " is below 0. ",
      "This indicates inconsistency in the reported data: the reported rate ",
      "is smaller than the dropout probability, which the stated method ",
      "makes impossible. Setting p_lower = 0.\n",
      "Please verify: p_obs = ", p_obs, ", omega = ", round(omega, 3),
      ", method = \"", method, "\", group = \"", group, "\""
    ))
    p_lower <- 0
  }

  if (p_upper > 1) {
    warning(paste0(
      "Upper bound p_upper = ", round(p_upper, 3), " exceeds 1. ",
      "This indicates inconsistency in the reported data: the stated ",
      "method may not have been applied, or p_obs and omega are ",
      "incompatible. Setting p_upper = 1.\n",
      "Please verify: p_obs = ", p_obs, ", omega = ", round(omega, 3),
      ", method = \"", method, "\", group = \"", group, "\""
    ))
    p_upper <- 1
  }

  result <- data.frame(
    method       = method,
    group        = group,
    p_obs        = p_obs,
    N_randomized = N_randomized,
    N_dropout    = N_dropout,
    omega        = omega,
    p_nri        = p_nri,
    p_lower      = p_lower,
    p_upper      = p_upper,
    p_midpoint   = (p_lower + p_upper) / 2,
    stringsAsFactors = FALSE
  )

  class(result) <- c("latent_bounds", "data.frame")

  return(result)
}

#' Print Method for latent_bounds Objects
#'
#' @param x An object of class \code{latent_bounds}.
#' @param ... Additional arguments (not used).
#'
#' @return \code{x}, invisibly.
#'
#' @export
print.latent_bounds <- function(x, ...) {
  cat("\nIdentified Range for the Latent Response Probability\n")
  cat("=====================================================\n\n")
  cat("Reported data:\n")
  cat(sprintf("         method = %s\n",     x$method))
  cat(sprintf("          group = %s\n",     x$group))
  cat(sprintf("          p_obs = %.4f\n",   x$p_obs))
  cat(sprintf("   N_randomized = %d\n",     x$N_randomized))
  cat(sprintf("      N_dropout = %d\n",     x$N_dropout))
  cat(sprintf("          omega = %.4f\n",   x$omega))
  cat(sprintf("          p_NRI = %.4f  (reported rate on the NRI scale)\n\n",
              x$p_nri))
  cat("Compatible range for the latent p:\n")
  cat(sprintf("        p_lower = %.4f\n", x$p_lower))
  cat(sprintf("     p_midpoint = %.4f\n", x$p_midpoint))
  cat(sprintf("        p_upper = %.4f\n", x$p_upper))
  cat(sprintf("          width = %.4f\n", x$p_upper - x$p_lower))
  invisible(x)
}

#' Plot Method for latent_bounds Objects
#'
#' Shows what the reported historical rate does and does not determine: across
#' the identified range for the latent response probability, the implied
#' correlation between the latent response and dropout, together with the range
#' that a bivariate Bernoulli distribution admits at each candidate value.
#'
#' @param x An object of class \code{latent_bounds}.
#' @param n_points Integer. Number of grid points across the range
#'   (default = 101).
#' @param ... Additional arguments (not used).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' lb <- latent_bounds(p_obs = 0.40, N_randomized = 200, N_dropout = 30)
#' plot(lb, n_points = 21)
#'
#' @importFrom rlang .data
#' @export
plot.latent_bounds <- function(x, n_points = 101, ...) {

  .require_ggplot2("The plot method")
  if (n_points < 3) stop("n_points must be at least 3")
  if (x$omega == 0) {
    stop("With no dropout the latent probability is identified; ",
         "there is nothing to plot", call. = FALSE)
  }

  omega <- x$omega
  p_nri <- x$p_nri

  grid <- seq(x$p_lower, x$p_upper, length.out = n_points)
  grid <- grid[grid > 0 & grid < 1]
  if (length(grid) < 3) {
    stop("The identified range is too close to 0 or 1 to plot", call. = FALSE)
  }

  implied <- (grid * (1 - omega) - p_nri) /
    sqrt(grid * (1 - grid) * omega * (1 - omega))
  bnd <- vapply(grid, function(g) {
    b <- rho_bounds(p = g, omega = omega)
    c(b$rho_lower, b$rho_upper)
  }, numeric(2))

  dat <- data.frame(
    p = rep(grid, 3),
    rho = c(implied, bnd[1, ], bnd[2, ]),
    series = factor(rep(c("Implied by the reported rate",
                          "Feasible range", "Feasible range"),
                        each = length(grid)),
                    levels = c("Implied by the reported rate",
                               "Feasible range")),
    grp = rep(c("implied", "lower", "upper"), each = length(grid))
  )

  ggplot2::ggplot(
    dat, ggplot2::aes(x = .data$p, y = .data$rho,
                      group = .data$grp, linetype = .data$series)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(
      x = "Candidate latent response probability",
      y = expression(rho), linetype = NULL,
      title = sprintf("Reported %.3f under %s in the %s group, omega = %.3f",
                      x$p_obs, toupper(x$method), x$group, omega),
      subtitle = sprintf("Latent p is identified only within [%.3f, %.3f]",
                         x$p_lower, x$p_upper)
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "top")
}

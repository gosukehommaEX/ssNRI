#' Feasible Range of the Dropout Response Probability
#'
#' Returns the range of values that the probability of responding among those
#' who drop out can take for given marginal probabilities. The range is the
#' Frechet-Prentice restriction of \code{\link{rho_bounds}} expressed on the
#' probability scale, where it reduces to the requirement that the four joint
#' cells be non-negative.
#'
#' @param p Numeric. Latent response probability, 0 < p < 1.
#' @param omega Numeric. Dropout probability, 0 < omega < 1.
#'
#' @return A \code{data.frame} of class \code{gamma_bounds} with one row
#' containing \code{p}, \code{omega}, \code{gamma_lower}, \code{gamma_upper}
#' and the equivalent correlation bounds \code{rho_lower} and \code{rho_upper}.
#'
#' @details
#' With gamma = Pr(R = 1 | D = 1) the four joint cells are
#' gamma omega, (1 - gamma) omega, p - gamma omega and 1 - p - omega +
#' gamma omega. Requiring all four to be non-negative gives
#'
#'   max(0, (p + omega - 1) / omega) <= gamma <= min(1, p / omega),
#'
#' which is the whole of the Frechet-Prentice restriction. The correlation
#' bounds of Prentice (1988) follow by applying \code{\link{gamma_to_rho}} to
#' the two endpoints. The upper bound changes branch at omega = p and the lower
#' bound at omega = 1 - p, exactly where the correlation bounds do.
#'
#' Reading the restriction on this scale makes its content plain. The upper
#' bound binds when dropouts cannot all have responded because there are not
#' enough responders in the group, and the lower bound binds when dropouts
#' cannot all have failed to respond because there are not enough
#' non-responders.
#'
#' @references
#' Prentice, R. L. (1988). Correlated binary regression with covariates specific
#' to each binary observation. \emph{Biometrics}, 44, 1033-1048.
#'
#' @examples
#' gamma_bounds(p = 0.45, omega = 0.10)
#' gamma_bounds(p = 0.20, omega = 0.30)
#'
#' @seealso \code{\link{rho_bounds}}, \code{\link{rho_to_gamma}},
#' \code{\link{gamma_to_rho}}
#'
#' @export
gamma_bounds <- function(p, omega) {

  # Input validation
  if (length(p) != 1 || length(omega) != 1) {
    stop("p and omega must each be of length one")
  }
  if (!is.finite(p) || p <= 0 || p >= 1) {
    stop("p must be strictly between 0 and 1")
  }
  if (!is.finite(omega) || omega <= 0 || omega >= 1) {
    stop("omega must be strictly between 0 and 1")
  }

  # Non-negativity of the four joint cells
  gamma_lower <- max(0, (p + omega - 1) / omega)
  gamma_upper <- min(1, p / omega)

  result <- data.frame(
    p           = p,
    omega       = omega,
    gamma_lower = gamma_lower,
    gamma_upper = gamma_upper,
    rho_lower   = gamma_to_rho(p, omega, gamma_lower),
    rho_upper   = gamma_to_rho(p, omega, gamma_upper)
  )

  class(result) <- c("gamma_bounds", "data.frame")

  return(result)
}

#' Print Method for gamma_bounds Objects
#'
#' @param x An object of class \code{gamma_bounds}.
#' @param ... Additional arguments (not used).
#'
#' @return \code{x}, invisibly.
#'
#' @export
print.gamma_bounds <- function(x, ...) {
  cat("Feasible Range of the Dropout Response Probability\n")
  cat("==================================================\n")
  cat("Parameters:\n")
  cat(sprintf("              p = %.4f\n", x$p))
  cat(sprintf("          omega = %.4f\n", x$omega))
  cat("Feasible range for gamma = Pr(R = 1 | D = 1):\n")
  cat(sprintf("    gamma_lower = %.4f\n", x$gamma_lower))
  cat(sprintf("    gamma_upper = %.4f\n", x$gamma_upper))
  cat("Equivalent correlation range:\n")
  cat(sprintf("      rho_lower = %.4f\n", x$rho_lower))
  cat(sprintf("      rho_upper = %.4f\n", x$rho_upper))
  invisible(x)
}

#' Plot Method for gamma_bounds Objects
#'
#' The exchange rate between the two parameterizations: the probability of
#' responding among dropouts implied by each attainable value of the
#' correlation, at the latent response probability stored in the object and at
#' several dropout probabilities. Each segment spans the feasible range of the
#' correlation at its own dropout probability, so its vertical extent is the
#' attainable range of gamma there.
#'
#' The bounds themselves are not plotted, because they are almost always 0 and
#' 1: gamma is an ordinary conditional probability and the Frechet-Prentice
#' restriction binds only once the dropout probability exceeds p or 1 - p. What
#' does vary, sharply, is the slope. At a small dropout probability a
#' correlation of a few tenths already moves gamma across most of the unit
#' interval, so an assumption that looks mild on the correlation scale is a
#' strong statement about what the dropouts would have done.
#'
#' @param x An object of class \code{gamma_bounds}.
#' @param omega_values Numeric vector of dropout probabilities to draw, each
#'   strictly inside (0, 1). Defaults to the object's own dropout probability
#'   together with 0.05, 0.15 and 0.30.
#' @param n_points Integer. Number of grid points per segment (default = 101).
#' @param ... Additional arguments (not used).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' plot(gamma_bounds(p = 0.45, omega = 0.10), n_points = 21)
#'
#' @importFrom rlang .data
#' @export
plot.gamma_bounds <- function(x, omega_values = NULL, n_points = 101, ...) {

  .require_ggplot2("The plot method")
  if (is.null(omega_values)) {
    omega_values <- c(x$omega, 0.05, 0.15, 0.30)
  }
  omega_values <- sort(unique(omega_values))
  if (any(!is.finite(omega_values)) || any(omega_values <= 0) ||
      any(omega_values >= 1)) {
    stop("omega_values must be strictly inside (0, 1)")
  }
  if (n_points < 3) stop("n_points must be at least 3")

  dat <- do.call(rbind, lapply(omega_values, function(w) {
    b <- rho_bounds(p = x$p, omega = w)
    rho <- seq(b$rho_lower, b$rho_upper, length.out = n_points)
    data.frame(
      omega = w,
      rho = rho,
      gamma = pmin(pmax(rho_to_gamma(p = x$p, omega = w, rho = rho), 0), 1),
      stringsAsFactors = FALSE
    )
  }))
  dat$omega <- factor(formatC(dat$omega, digits = 2, format = "f"),
                      levels = formatC(omega_values, digits = 2, format = "f"))

  ggplot2::ggplot(
    dat, ggplot2::aes(x = .data$rho, y = .data$gamma, colour = .data$omega)
  ) +
    ggplot2::geom_hline(yintercept = x$p, linetype = "dotted") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted") +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      x = expression(rho), y = expression(gamma),
      colour = expression(omega),
      title = sprintf("Correlation and dropout response probability at p = %.2f",
                      x$p)
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "right")
}

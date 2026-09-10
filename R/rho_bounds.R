#' Feasible Range of the Response-Dropout Correlation
#'
#' Returns the feasible range of the correlation between the latent response
#' indicator and the dropout indicator for given marginal probabilities. Because
#' both variables are binary, the correlation is not free to range over the whole
#' interval from -1 to 1; the attainable range follows from the Frechet bounds
#' and was given for binary variables by Prentice (1988).
#'
#' @param p Numeric. Latent response probability, 0 < p < 1.
#' @param omega Numeric. Dropout probability, 0 <= omega < 1.
#'
#' @return A \code{data.frame} of class \code{rho_bounds} with one row containing
#' \code{p}, \code{omega}, \code{rho_lower} and \code{rho_upper}.
#'
#' @details
#' For binary R with P(R = 1) = p and binary D with P(D = 1) = omega, the
#' correlation is confined to
#'
#'   rho_lower = max(-sqrt(p omega / ((1 - p)(1 - omega))),
#'                   -sqrt((1 - p)(1 - omega) / (p omega)))
#'
#'   rho_upper = min( sqrt(p (1 - omega) / (omega (1 - p))),
#'                    sqrt(omega (1 - p) / (p (1 - omega))))
#'
#' When omega = 0 there is no dropout and the correlation is not defined; the
#' degenerate range from 0 to 0 is returned so that callers can treat this case
#' without a special branch.
#'
#' @references
#' Prentice, R. L. (1988). Correlated binary regression with covariates specific
#' to each binary observation. \emph{Biometrics}, 44, 1033-1048.
#'
#' @examples
#' rho_bounds(p = 0.35, omega = 0.10)
#' rho_bounds(p = 0.60, omega = 0.25)
#'
#' @export
rho_bounds <- function(p, omega) {

  # Input validation
  if (p <= 0 || p >= 1) {
    stop("p must be strictly between 0 and 1")
  }
  if (omega < 0 || omega >= 1) {
    stop("omega must be between 0 (inclusive) and 1 (exclusive)")
  }

  # Special case: no dropout
  if (omega == 0) {
    result <- data.frame(
      p         = p,
      omega     = 0,
      rho_lower = 0,
      rho_upper = 0
    )
    class(result) <- c("rho_bounds", "data.frame")
    return(result)
  }

  # Prentice (1988) bounds
  lower1    <- -sqrt(p * omega / ((1 - p) * (1 - omega)))
  lower2    <- -sqrt((1 - p) * (1 - omega) / (p * omega))
  rho_lower <- max(lower1, lower2)

  upper1    <- sqrt(p * (1 - omega) / (omega * (1 - p)))
  upper2    <- sqrt(omega * (1 - p) / (p * (1 - omega)))
  rho_upper <- min(upper1, upper2)

  result <- data.frame(
    p         = p,
    omega     = omega,
    rho_lower = rho_lower,
    rho_upper = rho_upper
  )

  class(result) <- c("rho_bounds", "data.frame")

  return(result)
}

#' Print Method for rho_bounds Objects
#'
#' @param x An object of class \code{rho_bounds}.
#' @param ... Additional arguments (not used).
#'
#' @return \code{x}, invisibly.
#'
#' @export
print.rho_bounds <- function(x, ...) {
  cat("Feasible Correlation Bounds (Prentice 1988)\n")
  cat("============================================\n")
  cat("Parameters:\n")
  cat(sprintf("              p = %.4f\n", x$p))
  cat(sprintf("          omega = %.4f\n", x$omega))
  cat("Feasible range:\n")
  cat(sprintf("      rho_lower = %.4f\n", x$rho_lower))
  cat(sprintf("      rho_upper = %.4f\n", x$rho_upper))
  cat(sprintf("    range_width = %.4f\n", x$rho_upper - x$rho_lower))
  invisible(x)
}

#' Plot Method for rho_bounds Objects
#'
#' Feasible range of the response-dropout correlation as a function of the
#' dropout probability, at the latent response probability stored in the object.
#' A dashed vertical line marks the dropout probability the object itself was
#' evaluated at.
#'
#' @param x An object of class \code{rho_bounds}.
#' @param omega_range Numeric vector of length two. Range of dropout
#'   probabilities to sweep (default = c(0.01, 0.50)).
#' @param n_points Integer. Number of grid points (default = 101).
#' @param ... Additional arguments (not used).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' plot(rho_bounds(p = 0.35, omega = 0.10), n_points = 21)
#'
#' @importFrom rlang .data
#' @export
plot.rho_bounds <- function(x, omega_range = c(0.01, 0.50),
                            n_points = 101, ...) {

  .require_ggplot2("The plot method")
  if (length(omega_range) != 2 || omega_range[1] <= 0 ||
      omega_range[2] >= 1 || omega_range[2] <= omega_range[1]) {
    stop("omega_range must be two increasing values strictly inside (0, 1)")
  }
  if (n_points < 3) stop("n_points must be at least 3")

  grid <- seq(omega_range[1], omega_range[2], length.out = n_points)
  bnd <- vapply(grid, function(g) {
    b <- rho_bounds(p = x$p, omega = g)
    c(b$rho_lower, b$rho_upper)
  }, numeric(2))

  dat <- data.frame(
    omega = rep(grid, 2),
    rho = c(bnd[1, ], bnd[2, ]),
    bound = factor(rep(c("Lower", "Upper"), each = length(grid)),
                   levels = c("Upper", "Lower"))
  )

  out <- ggplot2::ggplot(
    dat, ggplot2::aes(x = .data$omega, y = .data$rho, group = .data$bound)
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
    ggplot2::labs(
      x = expression(omega), y = expression(rho),
      title = sprintf("Feasible correlation range at p = %.2f", x$p)
    ) +
    ggplot2::theme_bw(base_size = 12)

  if (x$omega > 0) {
    out <- out + ggplot2::geom_vline(xintercept = x$omega, linetype = "dashed")
  }

  out
}

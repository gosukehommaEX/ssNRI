#' Convert the Dropout Response Probability to the Response-Dropout Correlation
#'
#' Inverse of \code{\link{rho_to_gamma}}. Converts the probability of responding
#' among those who drop out into the correlation between the latent response
#' indicator and the dropout indicator.
#'
#' @param p Numeric vector. Latent response probability, 0 < p < 1.
#' @param omega Numeric vector. Dropout probability, 0 < omega < 1.
#' @param gamma Numeric vector. Probability of responding among dropouts,
#'   0 <= gamma <= 1.
#'
#' @return A numeric vector of the same length as the recycled arguments,
#' giving the correlation.
#'
#' @details
#' Inverting gamma = p + rho sqrt(p (1 - p) (1 - omega) / omega) gives
#'
#'   rho = (gamma - p) sqrt(omega / (p (1 - p) (1 - omega))).
#'
#' Values of gamma admissible for the given marginal probabilities map into the
#' Frechet-Prentice range returned by \code{\link{rho_bounds}}; see
#' \code{\link{gamma_bounds}}.
#'
#' @examples
#' gamma_to_rho(p = 0.45, omega = 0.10, gamma = 0.45)
#' gamma_to_rho(p = 0.45, omega = 0.10,
#'              gamma = rho_to_gamma(0.45, 0.10, c(-0.3, 0, 0.3)))
#'
#' @seealso \code{\link{rho_to_gamma}}, \code{\link{gamma_bounds}}
#'
#' @export
gamma_to_rho <- function(p, omega, gamma) {

  # Input validation
  if (any(!is.finite(p)) || any(p <= 0) || any(p >= 1)) {
    stop("p must be strictly between 0 and 1")
  }
  if (any(!is.finite(omega)) || any(omega <= 0) || any(omega >= 1)) {
    stop("omega must be strictly between 0 and 1")
  }
  if (any(!is.finite(gamma)) || any(gamma < 0) || any(gamma > 1)) {
    stop("gamma must be between 0 and 1")
  }

  rho <- (gamma - p) * sqrt(omega / (p * (1 - p) * (1 - omega)))

  return(rho)
}

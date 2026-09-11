#' Convert the Response-Dropout Correlation to the Dropout Response Probability
#'
#' Converts the correlation between the latent response indicator and the
#' dropout indicator into the probability of responding among those who drop
#' out. The two parameterizations describe the same bivariate Bernoulli
#' distribution, so the map is one to one for a fixed pair of marginal
#' probabilities.
#'
#' @param pi Numeric vector. Latent response probability, 0 < pi < 1.
#' @param omega Numeric vector. Dropout probability, 0 < omega < 1.
#' @param rho Numeric vector. Correlation between the latent response indicator
#'   and the dropout indicator.
#'
#' @return A numeric vector of the same length as the recycled arguments,
#' giving gamma = Pr(R = 1 | D = 1).
#'
#' @details
#' Writing R for the latent response indicator and D for the dropout indicator,
#' the joint cell in which a patient would have responded but dropped out is
#' pi omega + rho sqrt(pi (1 - pi) omega (1 - omega)). Dividing by omega gives
#'
#'   gamma = pi + rho sqrt(pi (1 - pi) (1 - omega) / omega).
#'
#' Independence corresponds to gamma = pi, a positive correlation to dropouts
#' who respond more often than the group as a whole, and a negative correlation
#' to dropouts who respond less often. The NRI estimand is
#' (pi1 - gamma1 omega1) - (pi0 - gamma0 omega0), so gamma states directly the
#' quantity that non-responder imputation gets wrong: it counts every dropout as
#' a non-responder, and gamma is the fraction of them for whom that is false.
#'
#' No feasibility check is applied. Use \code{\link{gamma_bounds}} for the range
#' that a bivariate Bernoulli distribution admits.
#'
#' @examples
#' rho_to_gamma(pi = 0.45, omega = 0.10, rho = 0)
#' rho_to_gamma(pi = 0.45, omega = 0.10, rho = c(-0.3, 0, 0.3))
#'
#' @seealso \code{\link{gamma_to_rho}}, \code{\link{gamma_bounds}},
#' \code{\link{rho_bounds}}
#'
#' @export
rho_to_gamma <- function(pi, omega, rho) {

  # Input validation
  if (any(!is.finite(pi)) || any(pi <= 0) || any(pi >= 1)) {
    stop("pi must be strictly between 0 and 1")
  }
  if (any(!is.finite(omega)) || any(omega <= 0) || any(omega >= 1)) {
    stop("omega must be strictly between 0 and 1")
  }
  if (any(!is.finite(rho))) {
    stop("rho must be finite")
  }

  gamma <- pi + rho * sqrt(pi * (1 - pi) * (1 - omega) / omega)

  return(gamma)
}

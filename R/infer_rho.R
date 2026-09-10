#' Correlation Implied by a Reported Historical Rate
#'
#' Returns the correlation between the latent response and dropout indicators
#' that a historical trial's reported response rate implies, for an assumed
#' latent response probability.
#'
#' @param p_obs Numeric. Response probability as reported by the historical
#'   trial, on the scale implied by \code{method}, 0 <= p_obs <= 1.
#' @param omega Numeric. Dropout probability in that group,
#'   0 <= omega < 1.
#' @param p Numeric. Assumed latent response probability, 0 < p < 1. Required,
#'   because no missing-data handling rule identifies it; use
#'   \code{\link{latent_bounds}} to obtain the interval it lies in.
#' @param method Character. Missing-data handling rule the historical trial
#'   applied: \code{"nri"} (default), \code{"cc"}, \code{"bc"} (best case) or
#'   \code{"wc"} (worst case).
#' @param group Character. Which group the reported rate refers to,
#'   \code{"treatment"} (default) or \code{"control"}.
#'
#' @return A \code{data.frame} of class \code{infer_rho} with one row containing
#' \code{method}, \code{group}, \code{p_obs}, \code{omega}, \code{p} and
#' \code{rho}.
#'
#' @details
#' Given the marginal probabilities p and omega, one joint cell determines the
#' correlation. The reported rate fixes that cell once it has been placed on the
#' NRI scale, so
#'
#'   rho = (p (1 - omega) - p_NRI) / sqrt(p (1 - p) omega (1 - omega)),
#'
#' with p_NRI equal to p_obs under NRI and to p_obs (1 - omega) under a complete
#' case analysis.
#'
#' The assumed latent probability p is an input, not something the reported rate
#' supplies. A rate reported under any of these rules constrains p only to an
#' interval of width omega; use \code{\link{latent_bounds}} to obtain that
#' interval and treat the resulting correlations as a sensitivity range rather
#' than a point estimate. Carrying a single recovered value forward as though it
#' were estimated would treat the historical trial's counting convention as a
#' statement about its patients' unobserved outcomes.
#'
#' p_NRI is the reported rate placed on the NRI scale, which is an identity in
#' every case: p_obs when the rule counted dropouts as non-responders, which is
#' NRI in both groups, a best case analysis in the control group and a worst
#' case analysis in the treatment group; p_obs - omega when it counted them as
#' responders, which is a best case analysis in the treatment group and a worst
#' case analysis in the control group; and p_obs (1 - omega) under a complete
#' case analysis.
#'
#' No rule identifies the correlation on its own. Setting a joint cell to zero,
#' as is sometimes done to read a correlation off a best or worst case rate,
#' amounts to assuming that every dropout in that group would have responded;
#' that assumption is recovered here by passing the upper endpoint of the
#' interval that \code{\link{latent_bounds}} returns, and it is one point of
#' the interval rather than a consequence of the analysis rule.
#'
#' When the implied correlation falls outside the range a bivariate Bernoulli
#' distribution admits, a warning is issued: that combination of reported rate,
#' dropout probability and assumed p is not attainable.
#'
#' @examples
#' # Reported 40 percent under NRI with 15 percent dropout; the latent p is
#' # identified only within [0.40, 0.55], so sweep that range
#' lb <- latent_bounds(p_obs = 0.40, N_randomized = 200, N_dropout = 30)
#' for (p in c(lb$p_lower + 0.01, lb$p_midpoint, lb$p_upper - 0.01)) {
#'   print(infer_rho(p_obs = 0.40, omega = lb$omega, p = p, method = "nri"))
#' }
#'
#' @seealso \code{\link{latent_bounds}}, \code{\link{rho_bounds}}
#'
#' @export
infer_rho <- function(p_obs, omega, p = NULL,
                      method = c("nri", "cc", "bc", "wc"),
                      group  = c("treatment", "control")) {

  method <- match.arg(method)
  group  <- match.arg(group)

  # Input validation
  if (p_obs < 0 || p_obs > 1) stop("p_obs must be between 0 and 1")
  if (omega < 0 || omega >= 1) {
    stop("omega must be between 0 (inclusive) and 1 (exclusive)")
  }

  # Special case: no dropout
  if (omega == 0) {
    if (!is.null(p) && abs(p_obs - p) > 1e-10) {
      stop("With omega = 0, p_obs must equal p. Check input values.")
    }
    result <- data.frame(
      method = method, group = group,
      p_obs = p_obs, omega = 0,
      p = p_obs, rho = 0,
      stringsAsFactors = FALSE
    )
    class(result) <- c("infer_rho", "data.frame")
    return(result)
  }

  # Place the reported rate on the NRI scale; each conversion is an identity
  counts_dropouts_as_responders <-
    (method == "bc" && group == "treatment") ||
    (method == "wc" && group == "control")

  p_nri <- if (method == "cc") {
    p_obs * (1 - omega)
  } else if (counts_dropouts_as_responders) {
    p_obs - omega
  } else {
    p_obs
  }

  # The reported rate does not identify the latent probability under any rule,
  # so the assumed p is always an input
  if (is.null(p)) {
    stop(paste0(
      "p must be supplied for method = \"", method, "\" and group = \"",
      group, "\". Use latent_bounds() to establish the compatible range."
    ))
  }
  if (p <= 0 || p >= 1) stop("p must be strictly between 0 and 1")

  denominator <- sqrt(p * (1 - p) * omega * (1 - omega))
  if (denominator < 1e-10) {
    stop("Denominator too small. Check input parameters.")
  }
  rho_val      <- (p * (1 - omega) - p_nri) / denominator
  p_for_output <- p

  # Check against the feasible range
  bounds <- rho_bounds(p = p_for_output, omega = omega)
  ctol   <- 1e-6
  if (rho_val < bounds$rho_lower - ctol || rho_val > bounds$rho_upper + ctol) {
    warning(paste0(
      "Inferred rho = ", round(rho_val, 4),
      " is outside the feasible range [",
      round(bounds$rho_lower, 4), ", ", round(bounds$rho_upper, 4), "]. ",
      "Possible causes:\n",
      "  1. The assumed latent p lies outside the interval compatible with ",
      "p_obs, omega and this method; latent_bounds() returns that interval.\n",
      "  2. The stated method was not actually applied in the historical ",
      "trial.\n",
      "  3. There are other sources of inconsistency in the reported data."
    ))
  }

  result <- data.frame(
    method = method,
    group  = group,
    p_obs  = p_obs,
    omega  = omega,
    p      = p_for_output,
    rho    = rho_val,
    stringsAsFactors = FALSE
  )

  class(result) <- c("infer_rho", "data.frame")

  return(result)
}

#' Print Method for infer_rho Objects
#'
#' @param x An object of class \code{infer_rho}.
#' @param ... Additional arguments (not used).
#'
#' @return \code{x}, invisibly.
#'
#' @export
print.infer_rho <- function(x, ...) {
  cat("Correlation Implied by a Reported Historical Rate\n")
  cat("=================================================\n")
  cat("Inputs:\n")
  cat(sprintf("         method = %s\n",   x$method))
  cat(sprintf("          group = %s\n",   x$group))
  cat(sprintf("          p_obs = %.4f\n", x$p_obs))
  cat(sprintf("          omega = %.4f\n", x$omega))
  cat(sprintf("              p = %.4f  (assumed)\n", x$p))
  cat("Implied correlation:\n")
  cat(sprintf("            rho = %.4f\n", x$rho))
  bounds <- rho_bounds(p = x$p, omega = x$omega)
  cat("Feasible range at this p:\n")
  cat(sprintf("      rho_lower = %.4f\n", bounds$rho_lower))
  cat(sprintf("      rho_upper = %.4f\n", bounds$rho_upper))
  invisible(x)
}

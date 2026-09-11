# ssNRI 0.2.0

* The names now transliterate the article's notation, so every argument, column,
  printed label and plot label corresponds one to one with a symbol in the
  article.
* The latent response probability is `pi`, written `pi1` and `pi0` for the two
  groups, matching the article's `pi_j`. `p` becomes `pi`, `p_obs` becomes
  `pi_obs`, `p1_NRI` becomes `pi1_NRI`, `p_lower` becomes `pi_lower`, and so on
  throughout. The package does not use the constant `pi` anywhere, so nothing is
  shadowed in a way that changes a result.
* The experimental group is 1 and the control group is 0, matching the article.
  `omega2` becomes `omega0`, `rho2` becomes `rho0`, `n2` becomes `n0`, and the
  allocation ratio `r` is now `n1 / n0`. The argument order is unchanged, the
  experimental group still coming first, so positional calls keep their meaning.
* The quantity assumed under the null hypothesis is `pi_null`, and the pooled
  values used internally are `pi_bar_CC` in `power_cc()`, `pi_bar` in
  `power_nri()` and `omega_bar` for the pooled dropout probability. A 0 subscript
  now means the control group and nothing else.
* The column `effect_latent` is renamed `effect_full`, matching the article's
  renaming of the dropout-free design from "traditional" to "full data".
* Plot titles, subtitles and axis labels are drawn as plotmath, so they show the
  Greek letters the article uses rather than their ASCII names.
* This is a breaking change for every call that names an argument. No numerical
  result changes; the renaming is confined to names.

# ssNRI 0.1.0

* Initial version.
* `power_cc()` and `power_nri()` evaluate exact power and exact type I error rate
  by complete enumeration of the joint outcome distribution.
* `power_cc()` uses a monotonicity property of the pooled two-sample Z statistic
  to collapse the enumeration over responder pairs into a single pass over the
  control-group counts. The returned values are identical to those of a full
  double enumeration; only the order of summation differs.
* `sample_size_cc()` implements the simple inflation method, and, with both
  dropout probabilities set to zero, the design without dropout.
  `sample_size_nri()` implements the closed-form sample size on the NRI scale.
* `latent_bounds()` returns the range of latent response probabilities
  compatible with a rate reported by a historical trial under a stated
  missing-data handling rule, covering NRI, complete case, best case and worst
  case, the last two of which count dropouts differently in the two groups;
  `infer_rho()` returns the correlation implied by that rate at an assumed
  latent probability. No rule identifies the latent probability, only an
  interval of width omega, so these are intended to be used together as a
  sensitivity range.
* `rho_bounds()` returns the feasible range of the response-dropout correlation
  for given marginal probabilities.
* `rho_to_gamma()`, `gamma_to_rho()` and `gamma_bounds()` express the same
  bivariate Bernoulli distribution through gamma, the probability of responding
  among those who drop out. The NRI response probability is then pi - gamma
  omega, and the Frechet-Prentice restriction becomes the requirement that the
  four joint cells be non-negative.
* `power_cc()` and `power_nri()` take `pi_null`, the common response probability
  assumed under the null hypothesis. The default reproduces the sample-size
  weighted average used previously; supplying a value makes the dependence of
  the size on this nuisance parameter visible.
* `print()` methods for every class, and `plot()` methods for `power_cc`,
  `power_nri`, `sample_size_nri`, `latent_bounds`, `rho_bounds` and
  `gamma_bounds`.
* A vignette works through a whole design on the trial the article uses as its
  application: the two estimands, the two ways of stating the association, the
  required sample size across the feasible range of the association, the exact
  power the simple inflation method actually delivers, and the interval that a
  historical rate leaves the latent response probability in.
* `inst/scripts/` and `inst/sup_info/` regenerate every figure and table of the
  article and its supplementary material, and `inst/scripts/number_check.R`
  re-extracts every number quoted in the article from the regenerated files.

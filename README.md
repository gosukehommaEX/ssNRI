# ssNRI

<!-- badges: start -->
[![R-CMD-check](https://github.com/gosukehommaEX/ssNRI/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gosukehommaEX/ssNRI/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Sample size and operating characteristics for two-group superiority trials with
a binary endpoint in which dropouts are handled by **non-responder imputation
(NRI)**.

The usual practice is to compute the sample size for a trial without dropout and
then inflate it by `1 / (1 - omega)`. That adjustment is justified for a complete
case analysis when dropout is unrelated to the outcome and occurs at the same
rate in both groups. Under NRI it is not, because NRI targets a different
estimand: a dropout is counted as a non-responder, so the effect that the trial
is actually powered to detect is attenuated. `ssNRI` makes that distinction
explicit and provides the corresponding sample size formula.

The latent response indicator and the dropout indicator are modelled jointly as
a bivariate Bernoulli distribution, so the group-specific dropout probabilities,
the association between response and dropout, and the allocation ratio can all
be specified directly.

This package accompanies the article

> Homma, G. A Cautionary Note on Sample Size Inflation for Trials with
> Non-Responder Imputation. (submitted)

## Installation

``` r
# install.packages("pak")
pak::pak("gosukehommaEX/ssNRI")
```

The vignette walks through a whole design on the trial the article uses as its
application. It is not built by the command above, so install with `remotes` if
you want it.

``` r
# install.packages("remotes")
remotes::install_github("gosukehommaEX/ssNRI", build_vignettes = TRUE)
vignette("ssNRI")
```

From a local checkout, `devtools::build_vignettes()` only renders into `doc/`,
which is not installed, so `vignette("ssNRI")` finds nothing until the package
itself is installed with the vignette built.

``` r
devtools::install(build_vignettes = TRUE)
```

## Exact operating characteristics

`power_cc()` and `power_nri()` return exact power and exact type I error rate by
enumerating the joint distribution of the observed outcomes; nothing is
simulated and no normal approximation is used at the evaluation stage.

``` r
library(ssNRI)

power_cc(n1 = 60, n2 = 60, p1 = 0.6, p2 = 0.4,
         omega1 = 0.2, omega2 = 0.2, rho1 = 0, rho2 = 0)
```

For a complete case analysis the enumeration runs over the number of completers
in each group and, within that, over the number of responders among completers.
A direct implementation is O(N^4). `power_cc()` instead uses the fact that the
pooled-variance Z statistic is strictly increasing in the group 1 responder count
wherever it is positive: the rejection region is an upper set, so the inner
double sum collapses to upper binomial tail probabilities evaluated at a single
threshold per control-group count. This is a regrouping of the same finite sum,
so the returned values are the exact ones; only the order of summation differs.

The optional argument `tol` additionally discards completer counts of negligible
probability. It defaults to `0`, which discards nothing. With `tol > 0` the
returned values are lower bounds for the exact ones and differ from them by at
most the discarded probability mass.

## Two ways of writing the same assumption

The association between response and dropout can be given as a correlation or as
`gamma`, the probability that a dropout would have responded. The two are in one
to one correspondence for fixed marginal probabilities, and `gamma` is often the
easier one to elicit, since the NRI response probability is `p - gamma * omega`:
`gamma` is exactly the fraction of dropouts that non-responder imputation
misclassifies.

``` r
rho_to_gamma(p = 0.45, omega = 0.10, rho = 0)      # 0.45, independence
gamma_bounds(p = 0.45, omega = 0.10)               # and the attainable range
```

## What a historical trial's reported rate does not tell you

A rate reported by a historical trial is the output of an analysis rule, not a
statement about the outcomes of its dropouts. `latent_bounds()` returns the
interval of latent response probabilities compatible with the reported value
under NRI, a complete case analysis, or a best or worst case analysis; the
interval always has width `omega`, and only its position depends on the rule.
`infer_rho()` then maps a candidate value from that interval to the correlation
it implies. Neither function identifies a point, and planning should carry the
interval through to the sample size.

``` r
lb <- latent_bounds(p_obs = 0.41, N_randomized = 620, N_dropout = 62)
sample_size_nri(p1 = lb$p_lower + 0.06, p2 = lb$p_lower,
                omega1 = 0.10, omega2 = 0.10)
```

## Reproducing the article

`inst/scripts/` regenerates the figures and tables of the main text and
`inst/sup_info/` those of the supplementary material. Each is split into a
computation stage that writes intermediate results to `data/*.rds` and a
rendering stage that reads those files and writes `results/*.eps`,
`results/*.pdf` and `results/*.tex`, so figures can be adjusted without
re-running any computation. `inst/scripts/number_check.R` re-extracts every
numerical value quoted in the article from the regenerated `.rds` files and
writes them to `inst/scripts/number_check_output.txt`, so that a number
transcribed into the manuscript can be checked against its source rather than
against memory.

Run them in this order, from the directory each script lives in.

``` r
setwd("inst/scripts")
source("data_generate_main.R")          # writes data/*.rds
source("table_and_figure_manuscript.R") # writes results/*.eps, *.pdf, *.tex
setwd("../sup_info")
source("data_generate_supplement.R")
source("table_and_figure_supplement.R")
setwd("../..")
source("inst/scripts/number_check.R")
```

## License

MIT

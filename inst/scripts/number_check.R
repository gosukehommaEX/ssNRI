# =============================================================================
# number_check.R
#
# Re-extraction of every numerical value quoted in the manuscript and the
# supplementary material, taken directly from the regenerated .rds files.
#
# Nothing here recomputes anything: the point is to read back what the
# computation stage actually produced, so that a value transcribed into the
# LaTeX source can be checked against its source rather than against memory.
#
# Usage
#   Run from the ssNRI package root, after running both
#   inst/scripts/data_generate_main.R and inst/sup_info/data_generate_supplement.R.
#
# Output
#   printed to the console and written to inst/scripts/number_check_output.txt
# =============================================================================

data_main <- "inst/scripts/data"
data_si <- "inst/sup_info/data"
out_path <- "inst/scripts/number_check_output.txt"
if (file.exists(out_path)) file.remove(out_path)

emit <- function(...) {
  line <- paste0(...)
  cat(line, "\n", sep = "")
  cat(line, "\n", sep = "", file = out_path, append = TRUE)
}

r2 <- function(x) formatC(round(as.numeric(x), 2), format = "f", digits = 2)
r3 <- function(x) formatC(round(as.numeric(x), 3), format = "f", digits = 3)
r4 <- function(x) formatC(round(as.numeric(x), 4), format = "f", digits = 4)
int <- function(x) formatC(as.numeric(x), format = "d", big.mark = ",")
rng3 <- function(x) paste0("[", r3(min(x)), ", ", r3(max(x)), "]")
rng4 <- function(x) paste0("[", r4(min(x)), ", ", r4(max(x)), "]")
rngi <- function(x) paste0("[", int(min(x)), ", ", int(max(x)), "]")

emit("==================================================================")
emit(" ssNRI: numbers quoted in the manuscript, re-extracted from .rds")
emit(" generated ", format(Sys.time()))
emit("==================================================================")

# =============================================================================
# Section 3.1 / Figure 1
# =============================================================================
f1 <- readRDS(file.path(data_main, "fig1_data.rds"))
d1 <- f1$df

emit("")
emit("################ Section 3.1 / Figure 1 ################")
emit("delta values: ", paste(f1$delta, collapse = ", "))
emit("p1 values: ", paste(f1$p1, collapse = ", "))
emit("omega values (", length(f1$omega_values), " of them): ",
     paste(f1$omega_values, collapse = ", "))
emit("  -> the text must say ", length(f1$omega_values),
     " values, starting at ", min(f1$omega_values))
emit("target power = ", f1$target_power, ", alpha = ", f1$alpha,
     ", r = ", f1$r)
emit("exact power across all scenarios:       ", rng4(d1$power))
emit("exact type I error across all scenarios:", rng4(d1$type1_error))
emit("largest departure from the target power:  ",
     r4(max(abs(d1$power - f1$target_power))))
emit("largest departure from the nominal alpha: ",
     r4(max(abs(d1$type1_error - f1$alpha))))
emit("total sample size range: ", rngi(d1$n_total))
for (dv in sort(unique(d1$delta_val))) {
  s <- d1[d1$delta_val == dv, ]
  emit("  delta = ", dv, ": power ", rng4(s$power),
       ", type I error ", rng4(s$type1_error))
}

# =============================================================================
# Section 3.2 / Figure 2
# =============================================================================
f2 <- readRDS(file.path(data_main, "fig2_data.rds"))
d2 <- f2$df

emit("")
emit("################ Section 3.2 / Figure 2 ################")
emit("delta = ", f2$delta, ", p1 values: ", paste(f2$p1, collapse = ", "),
     ", omega values: ", paste(f2$omega_values, collapse = ", "),
     ", rho step = ", f2$rho_step)
emit("grid points: ", nrow(d2))
emit("feasible rho range by case and omega:")
for (cs in unique(d2$case)) {
  for (om in sort(unique(d2$omega))) {
    s <- d2[d2$case == cs & d2$omega == om, ]
    if (nrow(s) == 0) next
    emit("  ", cs, ", omega = ", om, ": rho in ", rng3(s$rho),
         " (", nrow(s), " points)")
  }
}
emit("N under the NRI formula:      ", rngi(d2$N_proposed))
emit("N under the simple inflation: ", rngi(d2$N_Simple))
emit("exact power, NRI formula:      ", rng4(d2$power_proposed))
emit("exact power, simple inflation: ", rng4(d2$power_Simple))
emit("  -> the simple inflation method falls BELOW the target in ",
     sum(d2$power_Simple < f2$target_power), " of ", nrow(d2), " settings",
     " and ABOVE it in ", sum(d2$power_Simple > f2$target_power))
emit("  -> its worst shortfall is ", r4(min(d2$power_Simple)),
     " and its largest overshoot is ", r4(max(d2$power_Simple)))
worst <- d2[which.min(d2$power_Simple), ]
emit("     worst shortfall at ", worst$case, ", omega = ", worst$omega,
     ", rho = ", r3(worst$rho))
over <- d2[which.max(d2$power_Simple), ]
emit("     largest overshoot at ", over$case, ", omega = ", over$omega,
     ", rho = ", r3(over$rho))
emit("exact type I error, NRI formula:      ", rng4(d2$type1_proposed))
emit("exact type I error, simple inflation: ", rng4(d2$type1_Simple))

# =============================================================================
# Section 3.3 / Figure 3
# =============================================================================
f3 <- readRDS(file.path(data_main, "fig3_data.rds"))
d3 <- f3$df

emit("")
emit("################ Section 3.3 / Figure 3 ################")
emit("delta values: ", paste(f3$delta, collapse = ", "))
emit("p1 in [", min(f3$p1_values), ", ", max(f3$p1_values), "], omega in [",
     min(f3$omega_values), ", ", max(f3$omega_values), "]")
emit("grid points: ", nrow(d3))
emit("relative efficiency range: ", rng3(d3$re))
emit("  -> RE is below 1 in ", sum(d3$re < 1), " of ", nrow(d3),
     " grid points, at or above 1 in ", sum(d3$re >= 1))
for (dv in sort(unique(d3$delta))) {
  s <- d3[d3$delta == dv, ]
  emit("  delta = ", dv, ": RE ", rng3(s$re))
}

# =============================================================================
# Section 4 / Table 2: the cytisine application
# =============================================================================
t2 <- readRDS(file.path(data_main, "table2_data.rds"))
dt2 <- t2$df

emit("")
emit("################ Section 4 / Table 2 (cytisine) ################")
emit("design assumptions: p1 = ", t2$p1, ", p2 = ", t2$p2,
     ", omega = ", t2$omega, ", r = ", t2$r,
     ", alpha = ", t2$alpha, ", target power = ", t2$target_power)
emit("delta_Trad = ", r2(dt2$effect_latent[1]))
emit("N_CC (simple inflation)      = ", int(dt2$N_cc[1]),
     "   [reported by Dogar et al.: ", int(t2$n_reported), "]")
emit("  match: ", identical(as.numeric(dt2$N_cc[1]), as.numeric(t2$n_reported)))
emit("N_CC per group = ", int(dt2$N_cc[1] / 2))
for (k in seq_len(nrow(dt2))) {
  emit("rho = ", dt2$rho_label[k], " (", r3(dt2$rho[k]), "): ",
       "gamma1 = ", r3(dt2$gamma1[k]),
       ", gamma2 = ", r3(dt2$gamma2[k]),
       ", p1_NRI = ", r4(dt2$p1_NRI[k]),
       ", p2_NRI = ", r4(dt2$p2_NRI[k]),
       ", delta_NRI = ", r4(dt2$effect_NRI[k]),
       ", N_NRI = ", int(dt2$N_nri[k]),
       ", RE = ", r2(dt2$RE[k]))
}
mcar <- dt2[dt2$rho_label == "0", ]
emit("shortfall at rho = 0: ", int(mcar$N_nri - mcar$N_cc),
     " patients, i.e. ", r3(100 * (1 - mcar$RE)), " percent")

app <- readRDS(file.path(data_main, "application_data.rds"))
ap <- app$power
ob <- app$observed
emit("")
emit("################ Section 4 / trial results (AE-M5) ################")
emit("planned  ", int(app$n_planned), " (", int(ap$n1[1]), " per group); ",
     "enrolled ", int(app$n_enrolled), " (", int(ap$n1[2]), " and ",
     int(ap$n2[2]), ")")
emit("NRI formula requires ", int(ap$n_total[3]), " (", int(ap$n1[3]),
     " per group)")
emit("delta_NRI under the design assumptions = ", r4(ap$effect_NRI[1]))
for (k in seq_len(nrow(ap))) {
  emit("  exact NRI power at the ", ap$scenario[k], " size (",
       int(ap$n_total[k]), ") = ", r4(ap$power[k]))
}
emit("  -> intended power ", app$target_power, ", attained ",
     r4(ap$power[1]), " at the planned size and ", r4(ap$power[2]),
     " at the size actually enrolled")
emit("observed (Dogar et al.): ", ob$x1, "/", ap$n1[2], " = ", r3(ob$rate1),
     " vs ", ob$x2, "/", ap$n2[2], " = ", r3(ob$rate2))
emit("  risk difference ", r4(ob$risk_difference), " (95% CI ",
     r4(ob$ci_lower), " to ", r4(ob$ci_upper), "), p = ", ob$p_value)
emit("  observed dropout ", r4(ob$omega1_observed), " and ",
     r4(ob$omega2_observed), " against the assumed ", app$omega)
emit("  -> the observed difference ", r3(ob$risk_difference), " is ",
     r2(ob$risk_difference / (app$p1 - app$p2)),
     " times the assumed ", r2(app$p1 - app$p2),
     ", so the trial's null result is mostly a smaller than assumed effect;")
emit("     the design-analysis mismatch is a separate, additional loss of ",
     r3(100 * (app$target_power - ap$power[1])), " percentage points of power")

# =============================================================================
# Supplementary material
# =============================================================================
emit("")
emit("==================================================================")
emit(" Supplementary material")
emit("==================================================================")

s1 <- readRDS(file.path(data_si, "tableS1_data.rds"))$df
emit("")
emit("################ Table S1 (application sensitivity) ################")
emit("scenarios: ", nrow(s1), " rows (",
     length(unique(paste(s1$p1, s1$p2, s1$omega))), " parameter sets x 3 rho)")
emit("N_CC range:  ", rngi(s1$N_cc))
emit("N_NRI range: ", rngi(s1$N_nri))
emit("RE range:    ", rng3(s1$RE))
emit("  -> RE is below 1 in all ", sum(s1$RE < 1), " of ", nrow(s1), " rows: ",
     all(s1$RE < 1))
emit("delta_NRI range: ", rng4(s1$effect_NRI),
     "  (delta_Trad = ", r2(s1$effect_latent[1]), " throughout: ",
     length(unique(s1$effect_latent)) == 1, ")")
key <- s1[s1$p1 == 0.47 & s1$omega == 0.10, ]
emit("cross-check against main-text Table 2 (p1 = 0.47, omega = 0.10):")
for (k in seq_len(nrow(key))) {
  emit("  rho = ", key$rho_label[k], ": N_NRI = ", int(key$N_nri[k]),
       ", RE = ", r2(key$RE[k]))
}

fs1 <- readRDS(file.path(data_si, "figS1_data.rds"))
ds1 <- fs1$df
emit("")
emit("################ Figure S1 (unequal dropout, alternative) ################")
emit("p1 = ", fs1$p1, ", p2 = ", fs1$p2, ", rho = 0")
emit("omega1 values: ", paste(fs1$omega1_values, collapse = ", "),
     "; omega2 in [", min(ds1$omega2), ", ", max(ds1$omega2), "]")
emit("delta_NRI range: ", rng4(ds1$effect_NRI),
     " (delta_Trad = ", r2(ds1$effect_latent[1]), ")")
emit("N under the NRI formula:      ", rngi(ds1$N_proposed))
emit("N under the simple inflation: ", rngi(ds1$N_Simple))
emit("exact power, simple inflation: ", rng4(ds1$power_Simple))
w <- ds1[which.min(ds1$power_Simple), ]
emit("  worst case: omega1 = ", w$omega1, ", omega2 = ", w$omega2,
     " -> power ", r4(w$power_Simple),
     ", N_NRI = ", int(w$N_proposed), " vs N_CC = ", int(w$N_Simple),
     " (ratio ", r2(w$N_proposed / w$N_Simple), ")")
b <- ds1[which.max(ds1$power_Simple), ]
emit("  best case:  omega1 = ", b$omega1, ", omega2 = ", b$omega2,
     " -> power ", r4(b$power_Simple))
emit("exact power, NRI formula: ", rng4(ds1$power_proposed),
     "  (target ", fs1$target_power, ")")

fs2 <- readRDS(file.path(data_si, "figS2_data.rds"))
ds2 <- fs2$df
emit("")
emit("################ Figure S2 (unequal dropout, full-data null) ################")
emit("p1 = p2 in {", paste(fs2$p_values, collapse = ", "), "}, omega1 = ",
     fs2$omega1, ", n = ", fs2$n_per_group, " per group, alpha = ", fs2$alpha)
# seq() with a step produces values that are not exactly representable, so the
# match against omega1 is made with a tolerance rather than with ==
eq <- ds2[abs(ds2$omega2 - fs2$omega1) < 1e-8, ]
if (nrow(eq) == 0) {
  emit("WARNING: the omega2 grid does not contain omega1 = ", fs2$omega1)
} else {
  emit("at omega2 = omega1: delta_NRI = ", rng4(eq$effect_NRI),
       ", rejection probability = ", rng4(eq$reject_full_null))
  emit("  -> at equal dropout the NRI effect is zero and the test is at level: ",
       all(abs(eq$effect_NRI) < 1e-12))
}
for (pv in fs2$p_values) {
  s <- ds2[ds2$p == pv, ]
  emit("  p = ", pv, ": delta_NRI ", rng4(s$effect_NRI),
       ", rejection probability ", rng4(s$reject_full_null))
}
emit("size against the NRI null (should stay near alpha): ",
     rng4(ds2$size_nri_null))

fs3 <- readRDS(file.path(data_si, "figS3_data.rds"))
ds3 <- fs3$df
emit("")
emit("################ Figure S3 (unequal correlations) ################")
emit("p1 = ", fs3$p1, ", p2 = ", fs3$p2, ", omega = ", fs3$omega)
emit("rho1 in ", rng3(ds3$rho1), ", rho2 in ", rng3(ds3$rho2))
emit("N_CC = ", int(ds3$N_cc[1]), " throughout (it does not depend on rho): ",
     length(unique(ds3$N_cc)) == 1)
emit("N_NRI range: ", rngi(ds3$N_nri))
emit("RE range: ", rng3(ds3$re))
emit("  -> RE exceeds 1 in ", sum(ds3$re > 1), " of ", nrow(ds3),
     " grid points (", r3(100 * mean(ds3$re > 1)), " percent)")
emit("  -> at rho1 = rho2 = 0 the RE is ",
     r3(ds3$re[which.min(abs(ds3$rho1) + abs(ds3$rho2))]))
dg3 <- fs3$diagonal
emit("  -> on the diagonal rho1 = rho2 (", nrow(dg3), " points over rho in ",
     rng3(dg3$rho), ") the RE is ", rng3(dg3$re),
     ", below 1 at ", sum(dg3$re < 1), " of ", nrow(dg3))
emit("     maximum on the diagonal ", r3(max(dg3$re)), " at rho = ",
     r3(dg3$rho[which.max(dg3$re)]))
emit("     so every RE above 1 in this figure requires rho1 to differ from ",
     "rho2; the caption must say so")

fs45 <- readRDS(file.path(data_si, "figS45_data.rds"))
ds45 <- fs45$df
emit("")
emit("################ Figures S4 and S5 (effect ordering) ################")
emit("delta_Trad = ", r2(ds45$delta_latent[1]))
for (pn in unique(ds45$panel)) {
  s <- ds45[ds45$panel == pn, ]
  emit("  panel '", pn, "': delta_NRI - delta_Trad in ", rng4(s$diff),
       "; positive at ", r3(100 * mean(s$diff > 0)), " percent of the grid")
}
# The zero contour in the dropout panel is the line omega2 = (p1 / p2) omega1
emit("  dropout panel: the NRI effect exceeds the full-data effect when ",
     "omega2 / omega1 > p1 / p2 = ", r3(fs45$p1 / fs45$p2))

fs6 <- readRDS(file.path(data_si, "figS6_data.rds"))
ds6 <- fs6$df
emit("")
emit("################ Figure S6 (conservative dropout) ################")
emit("p1 = ", fs6$p1, ", p2 = ", fs6$p2, ", rho = 0")
emit("dropout-free sample size N_Trad = ", int(ds6$n_trad[1]))
emit("true omega in [", min(ds6$omega), ", ", max(ds6$omega), "]")
emit("assumption needed, omega*: ", rng3(ds6$omega_star))
emit("gap omega* - omega: ", rng3(ds6$inflation_gap))
for (om in c(0.05, 0.10, 0.20, 0.30)) {
  s <- ds6[abs(ds6$omega - om) < 1e-8, ]
  if (nrow(s) == 0) next
  emit("  true omega = ", om, ": N_CC = ", int(s$n_cc),
       ", N_NRI = ", int(s$n_nri),
       ", omega* = ", r3(s$omega_star),
       " (gap ", r3(s$inflation_gap), ")")
}

ts2 <- readRDS(file.path(data_si, "tableS2_data.rds"))$df
emit("")
emit("################ Table S2 (allocation ratio) ################")
emit("r values: ", paste(sort(unique(ts2$r)), collapse = ", "))
emit("RE range overall: ", rng3(ts2$RE))
for (key in unique(paste(ts2$p1, ts2$p2, ts2$omega))) {
  s <- ts2[paste(ts2$p1, ts2$p2, ts2$omega) == key, ]
  emit("  p = (", s$p1[1], ", ", s$p2[1], "), omega = ", s$omega[1],
       ": RE ", rng3(s$RE), "  (spread ", r3(max(s$RE) - min(s$RE)), ")")
}
emit("  -> largest spread in RE attributable to r, within a parameter set: ",
     r3(max(tapply(ts2$RE, paste(ts2$p1, ts2$p2, ts2$omega),
                   function(z) max(z) - min(z)))))

ts3 <- readRDS(file.path(data_si, "tableS3_data.rds"))
dt3 <- ts3$df
emit("")
emit("################ Table S3 (null nuisance parameter) ################")
emit("design: n1 = n2 = ", ts3$n1, ", p1 = ", ts3$p1, ", p2 = ", ts3$p2,
     ", omega = ", ts3$omega, ", alpha = ", ts3$alpha)
for (an in unique(dt3$analysis)) {
  s <- dt3[dt3$analysis == an, ]
  emit("  ", an, ": type I error ", rng4(s$type1_error),
       " over p_null in ", rng4(s$p_null))
  emit("    exceeds alpha at ", sum(s$type1_error > ts3$alpha), " of ",
       nrow(s), " values; maximum ", r4(max(s$type1_error)),
       " (", r3(100 * (max(s$type1_error) / ts3$alpha - 1)),
       " percent above nominal)")
  emit("    value at the weighted average used in the main text: ",
       r4(s$type1_error[s$is_default]))
}
emit("  -> the Z-test is not exact, so its size varies with the nuisance")
emit("     parameter and can sit slightly above the nominal level. This is a")
emit("     property of the test statistic, not of the sample size formula.")

ts4 <- readRDS(file.path(data_si, "tableS4_data.rds"))
dt4 <- ts4$df
emit("")
emit("################ Table S4 (unequal dropout counterexample) ################")
emit("design: p1 = ", ts4$p1, ", p2 = ", ts4$p2,
     ", rho1 = rho2 = 0, alpha = ", ts4$alpha,
     ", target power = ", ts4$target_power)
for (k in seq_len(nrow(dt4))) {
  emit("  omega = (", r2(dt4$omega1[k]), ", ", r2(dt4$omega2[k]), "): ",
       "n_complete = ", int(dt4$n_complete[k]),
       ", n = (", int(dt4$n1[k]), ", ", int(dt4$n2[k]), "), ",
       "power = ", r4(dt4$power[k]),
       ", type I error = ", r4(dt4$type1_error[k]))
}
emit("exact power across the counterexamples: ", rng4(dt4$power))
emit("  -> the target power ", ts4$target_power,
     " is attained however unequal the two dropout probabilities are,")
emit("     so equal dropout probabilities are sufficient but not necessary")

f7 <- readRDS(file.path(data_si, "figS7_data.rds"))
d7 <- f7$df
emit("")
emit("################ Figure S7 (partial identification) ################")
emit("historical rate ", f7$p_obs, " with omega = ", f7$omega,
     " from ", int(f7$n_hist), " randomized; targeted difference ",
     f7$delta)
for (rl in unique(d7$rule)) {
  sub <- d7[d7$rule == rl, ]
  emit("  ", rl, ": latent p2 in [", r3(sub$p_lower[1]), ", ",
       r3(sub$p_upper[1]), "], width ",
       r3(sub$p_upper[1] - sub$p_lower[1]))
  for (cv in unique(sub$convention)) {
    z <- sub[sub$convention == cv & sub$feasible, ]
    emit("      ", cv, ": N_CC ", rngi(z$n_cc), ", N_NRI ", rngi(z$n_nri),
         ", RE ", rng3(z$RE), " (", nrow(z), " of ",
         sum(sub$convention == cv), " grid points attainable)")
  }
}
fea <- d7[d7$feasible, ]
emit("RE below 1 at ", sum(fea$RE < 1), " of ", nrow(fea),
     " attainable grid points; overall RE ", rng3(fea$RE))
for (cv in unique(fea$convention)) {
  z <- fea[fea$convention == cv, ]
  emit("  ", cv, ": RE ", rng3(z$RE), ", below 1 at ", sum(z$RE < 1),
       " of ", nrow(z))
}
emit("  -> where the latent probability sits inside its interval barely")
emit("     moves the answer; which counting rule produced the historical")
emit("     rate, and how the two groups are tied together, move it a lot")

emit("")
emit("==================================================================")
emit(" end")
emit("==================================================================")

message("Numbers written to ", out_path)

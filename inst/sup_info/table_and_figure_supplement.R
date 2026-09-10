# =============================================================================
# table_and_figure_supplement.R
#
# Rendering stage for the supplementary material of:
#   "A Cautionary Note on Sample Size Inflation for Trials with
#    Non-Responder Imputation"
#
# Reads the intermediate data written by data_generate_supplement.R
# (data/*.rds) and produces the figures and LaTeX tables. No numerical
# computation happens here.
#
# Prerequisite
#   Run data_generate_supplement.R first, from this same working directory.
#
# Required packages: ggplot2, dplyr, tidyr, metR
#
# Outputs (relative to the working directory)
#   results/FigureS1.eps ... results/FigureS7.eps  (and matching .pdf)
#   results/tableS1_sensitivity.tex
#   results/tableS2_allocation.tex
#   results/tableS3_null_prevalence.tex
#   results/tableS4_unequal_dropout.tex
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(metR)

data_dir <- "data"
results_dir <- "results"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

tableS1_data <- readRDS(file.path(data_dir, "tableS1_data.rds"))
tableS2_data <- readRDS(file.path(data_dir, "tableS2_data.rds"))
tableS3_data <- readRDS(file.path(data_dir, "tableS3_data.rds"))
tableS4_data <- readRDS(file.path(data_dir, "tableS4_data.rds"))
figS1_data <- readRDS(file.path(data_dir, "figS1_data.rds"))
figS2_data <- readRDS(file.path(data_dir, "figS2_data.rds"))
figS3_data <- readRDS(file.path(data_dir, "figS3_data.rds"))
figS45_data <- readRDS(file.path(data_dir, "figS45_data.rds"))
figS6_data <- readRDS(file.path(data_dir, "figS6_data.rds"))
figS7_data <- readRDS(file.path(data_dir, "figS7_data.rds"))

OKABE <- c("#E69F00", "#56B4E9", "#009E73",
           "#F0E442", "#0072B2", "#D55E00")

# Diverging palette for the signed effect-difference maps: purple below zero,
# white at zero, green above
DIVERGING <- c("#40004B", "#762A83", "#9970AB", "#C2A5CF", "#E7D4E8",
               "#FFFFFF",
               "#D9F0D3", "#A6DBA0", "#5AAE61", "#1B7837", "#00441B")

# Diverging palette for the relative efficiency map, centred at RE = 1: orange
# where the simple inflation method falls short, blue where it overshoots
RE_PALETTE <- c("#7F3B08", "#B35806", "#E08214", "#FDB863", "#FEE0B6",
                "#FFFFFF",
                "#D8DAEB", "#B2ABD2", "#8073AC", "#542788", "#2D004B")

theme_panel <- function(base = 22) {
  theme_bw(base_size = base) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.key.width = unit(2.5, "cm"),
      strip.text = element_text(size = base + 2, face = "bold"),
      axis.title = element_text(size = base + 2, face = "bold"),
      legend.title = element_text(size = base + 2, face = "bold")
    )
}

save_figure <- function(plot, name, width = 16, height = 9, dpi = 600) {
  for (ext in c("eps", "pdf")) {
    ggsave(
      filename = file.path(results_dir, paste(name, ext, sep = ".")),
      plot = plot, device = ext, dpi = dpi, width = width, height = height
    )
  }
  message("wrote ", name, ".eps and ", name, ".pdf")
}

# A tabular whose column specification disagrees with the number of ampersands
# in its rows fails only at compile time, in a file nobody reads. The check runs
# on every table this script writes, and stops rather than warns: a table that
# cannot compile is not an output.
check_tabular <- function(lines, path) {
  i <- grep("\\begin{tabular}", lines, fixed = TRUE)
  if (length(i) != 1) return(invisible(NULL))
  spec <- sub("\\}.*$", "", sub("^.*\\{tabular\\}\\{", "", lines[i]))
  n_spec <- nchar(gsub("[^lcr]", "", spec))
  rows <- lines[endsWith(lines, "\\\\")]
  if (length(rows) == 0) return(invisible(NULL))
  n_row <- vapply(rows, function(r) {
    m <- gregexpr("&", r, fixed = TRUE)[[1]]
    sum(m > 0) + 1L
  }, integer(1))
  bad <- which(n_row != n_spec)
  if (length(bad) > 0) {
    stop(sprintf(
      "%s: the tabular specification '%s' has %d columns, but %d row(s) have %s. LaTeX would not compile this.",
      path, spec, n_spec, length(bad),
      paste(sort(unique(n_row[bad])), collapse = " or ")
    ), call. = FALSE)
  }
  invisible(NULL)
}

write_tex <- function(lines, path) {
  check_tabular(lines, path)
  writeLines(lines, path)
  message("wrote ", path)
}

fmt4 <- function(x) formatC(x, digits = 4, format = "f")
fmt2 <- function(x) formatC(x, digits = 2, format = "f")
fmt3 <- function(x) formatC(x, digits = 3, format = "f")
fmt_i <- function(x) formatC(x, format = "d", big.mark = ",")

# Shared builder for the filled contour maps.
#
# The breaks are centred on the value that separates the two regimes, and the
# palette diverges from that centre, so the threshold is visible without
# reading the colour bar. A ratio is passed in on the log scale so that the
# centre sits at zero and the two directions are treated alike; label_fun then
# converts the tick positions back to the original scale.
#
# Two things keep the colour bar legible. The step is chosen so that the
# number of bands stays in single figures, and once there are more than nine
# ticks only every second one is labelled, which halves the text without
# coarsening the colours. The check at the end is a guard, not a suggestion:
# if it fires the figure needs a larger step, not a warning in the log.
contour_map <- function(df, zvar, palette, step, centre, legend_title,
                        xlab, ylab, digits = 2, label_fun = NULL) {
  z <- df[[zvar]]
  half <- ceiling(max(abs(z - centre), na.rm = TRUE) / step) * step
  brks <- seq(centre - half, centre + half, by = step)
  if (length(brks) > 15) {
    stop("contour_map: ", length(brks), " breaks for '", zvar,
         "'. The colour bar would be unreadable; increase step, or pass the ",
         "quantity on a log scale if it is a ratio.", call. = FALSE)
  }

  if (is.null(label_fun)) {
    label_fun <- function(b) formatC(b, digits = digits, format = "f")
  }
  thin <- function(b) {
    txt <- label_fun(b)
    if (length(b) > 9) txt[seq_along(b) %% 2 == 0] <- ""
    txt
  }

  df$.z <- z
  ggplot(df, aes(x = x, y = y, z = .z)) +
    metR::geom_contour_fill(breaks = brks, na.fill = FALSE) +
    geom_contour(breaks = centre, colour = "black", linewidth = 1.2) +
    scale_fill_stepsn(
      colors = palette, breaks = brks,
      limits = c(centre - half, centre + half), na.value = "gray80",
      name = legend_title,
      labels = thin,
      guide = guide_colorsteps(title.position = "top", title.hjust = 0.5,
                               barwidth = 36, barheight = 1.5,
                               even.steps = TRUE)
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = xlab, y = ylab) +
    theme_panel() +
    theme(panel.grid = element_blank())
}

# =============================================================================
# Figure S1: unequal dropout probabilities under the alternative (R2-M2b)
# =============================================================================

s1 <- figS1_data$df
target_power <- figS1_data$target_power

s1_long <- bind_rows(
  s1 %>%
    select(omega1, omega2, N_proposed, N_Simple) %>%
    pivot_longer(c(N_proposed, N_Simple), names_to = "method_raw",
                 values_to = "value") %>%
    mutate(metric = "Total sample size"),
  s1 %>%
    select(omega1, omega2, power_proposed, power_Simple) %>%
    pivot_longer(c(power_proposed, power_Simple), names_to = "method_raw",
                 values_to = "value") %>%
    mutate(metric = "Exact power under NRI")
) %>%
  mutate(
    method = ifelse(grepl("Simple", method_raw), "Simple inflation",
                    "NRI formula"),
    metric = factor(metric, levels = c("Total sample size",
                                       "Exact power under NRI")),
    omega1_factor = factor(omega1)
  )

omega1_labels <- lapply(figS1_data$omega1_values,
                        function(x) bquote(omega[1] == .(x)))

figS1 <- ggplot(s1_long, aes(x = omega2, y = value,
                             color = omega1_factor, linetype = method)) +
  geom_line(linewidth = 1.2) +
  geom_hline(
    data = filter(s1_long, metric == "Exact power under NRI"),
    aes(yintercept = target_power),
    color = "gray", linetype = "longdash", linewidth = 1
  ) +
  facet_wrap(vars(metric), scales = "free_y", nrow = 1) +
  scale_color_manual(values = OKABE[seq_along(figS1_data$omega1_values)],
                     breaks = figS1_data$omega1_values,
                     labels = omega1_labels) +
  scale_linetype_manual(values = c("NRI formula" = "solid",
                                   "Simple inflation" = "dashed")) +
  labs(
    x = expression(paste("Control group dropout probability (",
                         omega[2], ")")),
    y = NULL,
    color = "Treatment group dropout",
    linetype = "Method"
  ) +
  theme_panel()

save_figure(figS1, "FigureS1")

# =============================================================================
# Figure S2: unequal dropout probabilities under the full-data null (R2-M2b)
# =============================================================================

s2 <- figS2_data$df
alpha <- figS2_data$alpha

s2_long <- bind_rows(
  s2 %>% transmute(p, omega2, value = effect_NRI, metric = "NRI effect"),
  s2 %>% transmute(p, omega2, value = reject_full_null,
                   metric = "Rejection probability")
) %>%
  mutate(
    metric = factor(metric, levels = c("NRI effect", "Rejection probability")),
    p_factor = factor(p)
  )

p_labels <- lapply(figS2_data$p_values, function(x) bquote(p == .(x)))

figS2 <- ggplot(s2_long, aes(x = omega2, y = value, color = p_factor)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = figS2_data$omega1, color = "gray40",
             linetype = "longdash", linewidth = 1) +
  geom_hline(
    data = filter(s2_long, metric == "Rejection probability"),
    aes(yintercept = alpha),
    color = "gray", linetype = "longdash", linewidth = 1
  ) +
  geom_hline(
    data = filter(s2_long, metric == "NRI effect"),
    aes(yintercept = 0),
    color = "gray", linetype = "longdash", linewidth = 1
  ) +
  facet_wrap(vars(metric), scales = "free_y", nrow = 1) +
  scale_color_manual(values = OKABE[seq_along(figS2_data$p_values)],
                     breaks = figS2_data$p_values,
                     labels = p_labels) +
  labs(
    x = expression(paste("Control group dropout probability (",
                         omega[2], ")")),
    y = NULL,
    color = "Common latent response probability"
  ) +
  theme_panel()

save_figure(figS2, "FigureS2")

# =============================================================================
# Figure S3: relative efficiency over the two correlations (R1-M3a)
# =============================================================================

s3 <- figS3_data$df

# The relative efficiency spans more than an order of magnitude over the
# feasible rectangle, from about 0.13 to about 3.4, so it is mapped on the log
# scale. That puts RE = 1 at the centre, gives a halving and a doubling the
# same visual weight, and keeps the number of colour bands small. The colour
# bar is labelled with the relative efficiency itself, not its logarithm.
s3_map <- data.frame(x = s3$rho1, y = s3$rho2, log2re = log2(s3$re))

figS3 <- contour_map(
  s3_map, "log2re", palette = RE_PALETTE, step = 0.5, centre = 0,
  legend_title = expression(RE == N[CC] / N[NRI]),
  xlab = expression(paste("Treatment group correlation (", rho[1], ")")),
  ylab = expression(paste("Control group correlation (", rho[2], ")")),
  label_fun = function(b) formatC(2 ^ b, digits = 2, format = "f")
) +
  # Every point with a relative efficiency above one lies off this line, so the
  # dashed diagonal is what keeps the map from reading as a contradiction of
  # the main text, which assumes a common correlation
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "black", linewidth = 1) +
  annotate("point", x = 0, y = 0, size = 4, shape = 21,
           fill = "white", colour = "black")

save_figure(figS3, "FigureS3", width = 12, height = 12)

# =============================================================================
# Figure S4: effect ordering over the two dropout probabilities (R2-M1)
# =============================================================================

s4 <- filter(figS45_data$df, panel == "dropout")

figS4 <- contour_map(
  s4, "diff", palette = DIVERGING, step = 0.04, centre = 0,
  legend_title = expression(delta[NRI] - delta[Trad]),
  xlab = expression(paste("Treatment group dropout probability (",
                          omega[1], ")")),
  ylab = expression(paste("Control group dropout probability (",
                          omega[2], ")"))
)

save_figure(figS4, "FigureS4", width = 12, height = 12)

# =============================================================================
# Figure S5: effect ordering over the two correlations (R2-M1, R1-M3a)
# =============================================================================

s5 <- filter(figS45_data$df, panel == "correlation")

figS5 <- contour_map(
  s5, "diff", palette = DIVERGING, step = 0.025, centre = 0,
  legend_title = expression(delta[NRI] - delta[Trad]),
  xlab = expression(paste("Treatment group correlation (", rho[1], ")")),
  ylab = expression(paste("Control group correlation (", rho[2], ")"))
) +
  annotate("point", x = 0, y = 0, size = 4, shape = 21,
           fill = "white", colour = "black")

save_figure(figS5, "FigureS5", width = 12, height = 12)

# =============================================================================
# Figure S6: conservative dropout assumptions (R1-M4, AE-m15)
# =============================================================================

s6 <- figS6_data$df

s6_long <- bind_rows(
  s6 %>% transmute(omega, value = omega_star,
                   series = "Assumption needed to match the NRI sample size",
                   metric = "Dropout probability"),
  s6 %>% transmute(omega, value = omega,
                   series = "True dropout probability",
                   metric = "Dropout probability"),
  s6 %>% transmute(omega, value = n_nri, series = "NRI formula",
                   metric = "Total sample size"),
  s6 %>% transmute(omega, value = n_cc, series = "Simple inflation",
                   metric = "Total sample size")
) %>%
  mutate(metric = factor(metric, levels = c("Dropout probability",
                                            "Total sample size")),
         series = factor(series, levels = c(
           "Assumption needed to match the NRI sample size",
           "True dropout probability",
           "NRI formula",
           "Simple inflation"
         )))

# The two panels carry different pairs of series, so line type alone leaves
# two indistinguishable keys in the legend. Colour separates all four, as in
# the other supplementary figures; line type is kept as a redundant cue.
S6_COLOURS <- c(
  "Assumption needed to match the NRI sample size" = "#E69F00",
  "True dropout probability" = "#999999",
  "NRI formula" = "#0072B2",
  "Simple inflation" = "#D55E00"
)
S6_LINETYPES <- c(
  "Assumption needed to match the NRI sample size" = "solid",
  "True dropout probability" = "dashed",
  "NRI formula" = "solid",
  "Simple inflation" = "dashed"
)

figS6 <- ggplot(s6_long, aes(x = omega, y = value,
                             colour = series, linetype = series)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(vars(metric), scales = "free_y", nrow = 1) +
  scale_colour_manual(values = S6_COLOURS) +
  scale_linetype_manual(values = S6_LINETYPES) +
  labs(
    x = expression(paste("True dropout probability (", omega, ")")),
    y = NULL,
    colour = NULL,
    linetype = NULL
  ) +
  guides(
    colour = guide_legend(nrow = 2, byrow = TRUE),
    linetype = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme_panel() +
  theme(legend.text = element_text(size = 18))

save_figure(figS6, "FigureS6")

# =============================================================================
# Figure S7: what a historical rate does and does not pin down (R1-M2, R2-M3)
#
# Each curve spans its own rule's identified interval, so the horizontal extent
# of a curve is that interval. Under a common correlation the curves stop short
# of the right-hand endpoint: the correlation the control group's assumption
# implies is then larger than the treatment group can attain, which the gamma
# scale shows directly, since it would require more than every dropout in the
# treatment group to have responded.
# =============================================================================

s7 <- figS7_data$df %>%
  mutate(
    rule = factor(rule, levels = c("Dropouts counted as non-responders",
                                   "Complete case",
                                   "Dropouts counted as responders")),
    convention = factor(convention,
                        levels = c("Common correlation",
                                   "Common dropout response probability"))
  )

S7_COLOURS <- c(
  "Dropouts counted as non-responders" = "#0072B2",
  "Complete case"                      = "#009E73",
  "Dropouts counted as responders"     = "#D55E00"
)

figS7 <- ggplot(s7, aes(x = p2, y = RE, colour = rule)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 1.2, na.rm = TRUE) +
  facet_wrap(vars(convention), nrow = 1) +
  scale_colour_manual(values = S7_COLOURS) +
  labs(
    x = expression(paste("Candidate latent response probability ", p[2])),
    y = expression(paste("RE = ", N[CC], " / ", N[NRI])),
    colour = NULL
  ) +
  guides(colour = guide_legend(nrow = 3, byrow = TRUE)) +
  theme_panel() +
  theme(legend.text = element_text(size = 18))

save_figure(figS7, "FigureS7")

# =============================================================================
# Table S1: sensitivity grid for the cytisine application (AE-M4)
# =============================================================================

t1 <- tableS1_data$df

rho_cell <- function(label, value) {
  if (label == "L") {
    sprintf("$L$ (%.3f)", value)
  } else if (label == "U") {
    sprintf("$U$ (%.3f)", value)
  } else {
    "$0$"
  }
}

tabS1 <- c(
  "\\begin{table}[htbp]",
  "  \\centering",
  "  \\footnotesize",
  paste0(
    "  \\caption{Required sample sizes and relative efficiency for the ",
    "cytisine trial across a range of assumed response and dropout ",
    "probabilities. Fixed parameters: $r = ", tableS1_data$r,
    "$, $\\alpha = ", tableS1_data$alpha, "$ (one-sided), target power $= ",
    tableS1_data$target_power, "$.}"
  ),
  "  \\label{tab:application-sensitivity}",
  "  \\begin{threeparttable}",
  "  \\begin{tabular}{cccccccccrrr}",
  "    \\toprule",
  paste0(
    "    $p_{1}$ & $p_{2}$ & $\\omega$ & $\\rho$ & ",
    "$\\gamma_{1}$ & $\\gamma_{2}$ & ",
    "$p_{1,\\mathrm{NRI}}$ & $p_{2,\\mathrm{NRI}}$ & ",
    "$\\delta_{\\mathrm{NRI}}$ & $N_{\\mathrm{CC}}$ & ",
    "$N_{\\mathrm{NRI}}$ & RE \\\\"
  ),
  "    \\midrule"
)

prev_key <- ""
for (k in seq_len(nrow(t1))) {
  key <- paste(t1$p1[k], t1$p2[k], t1$omega[k])
  first <- key != prev_key
  if (first && k > 1) tabS1 <- c(tabS1, "    \\addlinespace")
  tabS1 <- c(tabS1, paste0(
    "    ",
    if (first) fmt2(t1$p1[k]) else "", " & ",
    if (first) fmt2(t1$p2[k]) else "", " & ",
    if (first) fmt2(t1$omega[k]) else "", " & ",
    rho_cell(t1$rho_label[k], t1$rho[k]), " & ",
    fmt3(t1$gamma1[k]), " & ", fmt3(t1$gamma2[k]), " & ",
    fmt4(t1$p1_NRI[k]), " & ", fmt4(t1$p2_NRI[k]), " & ",
    fmt4(t1$effect_NRI[k]), " & ",
    if (first) fmt_i(t1$N_cc[k]) else "", " & ",
    fmt_i(t1$N_nri[k]), " & ", fmt2(t1$RE[k]), " \\\\"
  ))
  prev_key <- key
}

tabS1 <- c(
  tabS1,
  "    \\bottomrule",
  "  \\end{tabular}",
  "  \\begin{tablenotes}[flushleft]\\footnotesize",
  paste0(
    "    \\item $L$ and $U$ denote the lower and upper attainable values of ",
    "$\\rho$ given the marginal probabilities, and ",
    "$\\gamma_{j} = \\Pr(R_{j} = 1 \\mid D_{j} = 1)$ is the same ",
    "assumption expressed as the probability that a dropout in group $j$ ",
    "would have responded. ",
    "RE $= N_{\\mathrm{CC}} / N_{\\mathrm{NRI}}$."
  ),
  "  \\end{tablenotes}",
  "  \\end{threeparttable}",
  "\\end{table}"
)

write_tex(tabS1, file.path(results_dir, "tableS1_sensitivity.tex"))

# =============================================================================
# Table S2: allocation ratios other than one (R1-M3b)
# =============================================================================

t2 <- tableS2_data$df

tabS2 <- c(
  "\\begin{table}[htbp]",
  "  \\centering",
  paste0(
    "  \\caption{Effect of the allocation ratio on the required sample size, ",
    "under within-group independence between response and dropout. ",
    "Fixed parameters: $\\alpha = ", tableS2_data$alpha,
    "$ (one-sided), target power $= ", tableS2_data$target_power, "$.}"
  ),
  "  \\label{tab:allocation}",
  "  \\begin{threeparttable}",
  "  \\begin{tabular}{ccccrrr}",
  "    \\toprule",
  paste0(
    "    $p_{1}$ & $p_{2}$ & $\\omega$ & $r$ & ",
    "$N_{\\mathrm{CC}}$ & $N_{\\mathrm{NRI}}$ & RE \\\\"
  ),
  "    \\midrule"
)

prev_key <- ""
for (k in seq_len(nrow(t2))) {
  key <- paste(t2$p1[k], t2$p2[k], t2$omega[k])
  first <- key != prev_key
  if (first && k > 1) tabS2 <- c(tabS2, "    \\addlinespace")
  tabS2 <- c(tabS2, paste0(
    "    ",
    if (first) fmt2(t2$p1[k]) else "", " & ",
    if (first) fmt2(t2$p2[k]) else "", " & ",
    if (first) fmt2(t2$omega[k]) else "", " & ",
    formatC(t2$r[k], format = "g"), " & ",
    fmt_i(t2$N_cc[k]), " & ", fmt_i(t2$N_nri[k]), " & ",
    fmt2(t2$RE[k]), " \\\\"
  ))
  prev_key <- key
}

tabS2 <- c(
  tabS2,
  "    \\bottomrule",
  "  \\end{tabular}",
  "  \\begin{tablenotes}[flushleft]\\footnotesize",
  "    \\item $r = n_{1} / n_{2}$. RE $= N_{\\mathrm{CC}} / N_{\\mathrm{NRI}}$.",
  "  \\end{tablenotes}",
  "  \\end{threeparttable}",
  "\\end{table}"
)

write_tex(tabS2, file.path(results_dir, "tableS2_allocation.tex"))

# =============================================================================
# Table S3: type I error rate at several null nuisance parameters (R2-m1)
# =============================================================================

t3 <- tableS3_data$df

tabS3 <- c(
  "\\begin{table}[htbp]",
  "  \\centering",
  paste0(
    "  \\caption{Exact type I error rate as a function of the common response ",
    "probability assumed under the null hypothesis, at $n_{1} = n_{2} = ",
    tableS3_data$n1, "$ with $p_{1} = ", fmt2(tableS3_data$p1),
    "$, $p_{2} = ", fmt2(tableS3_data$p2),
    "$, $\\omega_{1} = \\omega_{2} = ", fmt2(tableS3_data$omega),
    "$, $\\rho_{1} = \\rho_{2} = 0$ and $\\alpha = ", tableS3_data$alpha,
    "$ (one-sided). The size of a two-sample binomial test is not constant in ",
    "this nuisance parameter, so no single value establishes size control.}"
  ),
  "  \\label{tab:null-prevalence}",
  "  \\begin{threeparttable}",
  "  \\begin{tabular}{llr}",
  "    \\toprule",
  "    Analysis & Assumed null probability & Type I error rate \\\\",
  "    \\midrule"
)

prev_analysis <- ""
for (k in seq_len(nrow(t3))) {
  first <- t3$analysis[k] != prev_analysis
  if (first && k > 1) tabS3 <- c(tabS3, "    \\addlinespace")
  tabS3 <- c(tabS3, paste0(
    "    ",
    if (first) t3$analysis[k] else "", " & ",
    fmt4(t3$p_null[k]),
    if (t3$is_default[k]) " (weighted average)" else "", " & ",
    fmt4(t3$type1_error[k]), " \\\\"
  ))
  prev_analysis <- t3$analysis[k]
}

tabS3 <- c(
  tabS3,
  "    \\bottomrule",
  "  \\end{tabular}",
  "  \\begin{tablenotes}[flushleft]\\footnotesize",
  paste0(
    "    \\item The sample size is the one the NRI formula gives under these ",
    "assumptions. For the NRI analysis the assumed probability is on the NRI ",
    "scale; for the complete case analysis it is the conditional response ",
    "probability among completers. The row marked as the weighted average is ",
    "the convention used in the main text."
  ),
  "  \\end{tablenotes}",
  "  \\end{threeparttable}",
  "\\end{table}"
)

write_tex(tabS3, file.path(results_dir, "tableS3_null_prevalence.tex"))

# =============================================================================
# Table S4: the three conditions are sufficient but not necessary (R2-m3, AE-M2)
# =============================================================================

t4 <- tableS4_data$df

tabS4 <- c(
  "\\begin{table}[htbp]",
  "  \\centering",
  paste0(
    "  \\caption{Exact power of a complete case analysis when the two dropout ",
    "probabilities differ and each group is inflated by its own dropout ",
    "probability. Design parameters: $p_{1} = ",
    fmt2(tableS4_data$p1), "$, $p_{2} = ", fmt2(tableS4_data$p2),
    "$, $\\rho_{1} = \\rho_{2} = 0$, $r = 1$, $\\alpha = ",
    tableS4_data$alpha, "$ (one-sided) and target power $= ",
    tableS4_data$target_power,
    "$. The target power is attained throughout, so equal dropout ",
    "probabilities are not necessary for the simple inflation method to be ",
    "valid.}"
  ),
  "  \\label{tab:unequal-dropout}",
  "  \\begin{threeparttable}",
  "  \\begin{tabular}{ccrrrrr}",
  "    \\toprule",
  paste0(
    "    $\\omega_{1}$ & $\\omega_{2}$ & $n_{\\mathrm{complete}}$ & ",
    "$n_{1}$ & $n_{2}$ & Power & Type I error \\\\"
  ),
  "    \\midrule"
)

for (k in seq_len(nrow(t4))) {
  tabS4 <- c(tabS4, paste0(
    "    ",
    fmt2(t4$omega1[k]), " & ", fmt2(t4$omega2[k]), " & ",
    fmt_i(t4$n_complete[k]), " & ",
    fmt_i(t4$n1[k]), " & ", fmt_i(t4$n2[k]), " & ",
    fmt4(t4$power[k]), " & ", fmt4(t4$type1_error[k]), " \\\\"
  ))
}

tabS4 <- c(
  tabS4,
  "    \\bottomrule",
  "  \\end{tabular}",
  "  \\begin{tablenotes}[flushleft]\\footnotesize",
  paste0(
    "    \\item $n_{\\mathrm{complete}}$ is the number of completers per ",
    "group the unadjusted formula requires, and $n_{j} = \\lceil ",
    "n_{\\mathrm{complete}} / (1 - \\omega_{j}) \\rceil$ is the number ",
    "randomized to group $j$. Power and type I error rate are exact, ",
    "computed by complete enumeration over the numbers of completers and ",
    "responders rather than by simulation."
  ),
  "  \\end{tablenotes}",
  "  \\end{threeparttable}",
  "\\end{table}"
)

write_tex(tabS4, file.path(results_dir, "tableS4_unequal_dropout.tex"))

message("All supplementary figures and tables written to '", results_dir, "/'.")

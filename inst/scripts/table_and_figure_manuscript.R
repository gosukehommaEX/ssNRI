# =============================================================================
# table_and_figure_manuscript.R
#
# Rendering stage for the main-text figures and table of:
#   "A Cautionary Note on Sample Size Inflation for Trials with
#    Non-Responder Imputation"
#
# Reads the intermediate data written by data_generate_main.R (data/*.rds) and
# produces the figures and the LaTeX table. No numerical computation happens
# here, so labels, scales, colours and captions can be adjusted without
# re-running anything.
#
# Prerequisite
#   Run data_generate_main.R first, from this same working directory.
#
# Required packages: ggplot2, dplyr, tidyr, ggh4x, metR
#
# Outputs (relative to the working directory)
#   results/Figure1.eps, results/Figure1.pdf
#   results/Figure2.eps, results/Figure2.pdf
#   results/Figure3.eps, results/Figure3.pdf
#   results/table2_application.tex
#
# The .eps files go to Overleaf. The .pdf files carry the same content and
# exist so that the output can be inspected directly.
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggh4x)
library(metR)

data_dir <- "data"
results_dir <- "results"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

fig1_data <- readRDS(file.path(data_dir, "fig1_data.rds"))
fig2_data <- readRDS(file.path(data_dir, "fig2_data.rds"))
fig3_data <- readRDS(file.path(data_dir, "fig3_data.rds"))
table2_data <- readRDS(file.path(data_dir, "table2_data.rds"))

# --------------------------------------------------------------------------- #
# Shared presentation constants
# --------------------------------------------------------------------------- #

# Okabe-Ito palette (colour-blind friendly)
OKABE <- c("#E69F00", "#56B4E9", "#009E73",
           "#F0E442", "#0072B2", "#D55E00")

METRIC_LEVELS <- c("Total Sample Size (N)", "Power", "Type I Error Rate")

# Shared theme for the three-row panel figures
theme_panel <- function() {
  theme_bw() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.key.width = unit(2.5, "cm"),
      strip.text = element_text(size = 26, face = "bold"),
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 28, face = "bold"),
      axis.title.x = element_text(size = 28, face = "bold"),
      legend.text = element_text(size = 30),
      legend.title = element_text(size = 30, face = "bold")
    )
}

# Save a figure to the results directory in both EPS and PDF
save_figure <- function(plot, name, width = 16, height = 16, dpi = 600) {
  for (ext in c("eps", "pdf")) {
    ggsave(
      filename = file.path(results_dir, paste(name, ext, sep = ".")),
      plot = plot,
      device = ext,
      dpi = dpi,
      width = width,
      height = height
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

# =============================================================================
# Figure 1: validity of the simple inflation method under a CC analysis
# =============================================================================

f1 <- fig1_data$df
target_power <- fig1_data$target_power
alpha <- fig1_data$alpha

fig1_long <- f1 %>%
  tidyr::pivot_longer(
    cols = c("n_total", "power", "type1_error"),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    metric = factor(
      metric,
      levels = c("n_total", "power", "type1_error"),
      labels = METRIC_LEVELS
    )
  )

case_levels <- unique(fig1_long$case)
case_labels_expr <- lapply(unique(fig1_long$pi1_value), function(pi_val) {
  case_idx <- which(fig1_data$pi1 == pi_val)
  bquote(paste("Case ", .(LETTERS[case_idx]), " (", pi[1] == .(pi_val), ")"))
})
fig1_long$case <- factor(fig1_long$case,
                         levels = case_levels,
                         labels = case_labels_expr)

delta_breaks <- sort(unique(fig1_long$delta_val))
delta_labels <- lapply(delta_breaks, function(x) bquote(delta[Full] == .(x)))
fig1_long$delta_factor <- factor(fig1_long$delta_val, levels = delta_breaks)

fig1 <- ggplot(fig1_long, aes(x = omega, y = value,
                              color = delta_factor)) +
  geom_line(linewidth = 1.2) +
  geom_hline(
    data = dplyr::filter(fig1_long, metric == "Power"),
    aes(yintercept = target_power),
    color = "gray", linetype = "longdash", linewidth = 1
  ) +
  geom_hline(
    data = dplyr::filter(fig1_long, metric == "Type I Error Rate"),
    aes(yintercept = alpha),
    color = "gray", linetype = "longdash", linewidth = 1
  ) +
  facet_grid(
    rows = vars(metric),
    cols = vars(case),
    scales = "free_y",
    labeller = labeller(case = label_parsed)
  ) +
  ggh4x::scale_y_facet(metric == "Total Sample Size (N)", limits = NULL) +
  ggh4x::scale_y_facet(metric == "Power", limits = c(0.7, 0.9)) +
  ggh4x::scale_y_facet(metric == "Type I Error Rate", limits = c(0.02, 0.03)) +
  scale_color_manual(
    values = OKABE[seq_along(delta_breaks)],
    breaks = delta_breaks,
    labels = delta_labels
  ) +
  scale_x_continuous(
    breaks = seq(0, max(fig1_data$omega_values), by = 0.1)
  ) +
  labs(
    x = expression(paste("Dropout Probability (", omega, ")")),
    y = NULL,
    color = "Treatment Effect"
  ) +
  theme_panel()

save_figure(fig1, "Figure1")

# =============================================================================
# Figure 2: proposed NRI formula against the simple inflation method
# =============================================================================

f2 <- fig2_data$df
target_power <- fig2_data$target_power
alpha <- fig2_data$alpha

to_long <- function(df, cols, metric_name) {
  df %>%
    dplyr::select(dplyr::all_of(c("rho", "case", "pi1_value", "omega", cols))) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(cols),
                        names_to = "method_raw", values_to = "value") %>%
    dplyr::mutate(
      method = ifelse(grepl("proposed", method_raw),
                      "Proposed", "Simple inflation"),
      metric = metric_name
    )
}

fig2_long <- dplyr::bind_rows(
  to_long(f2, c("N_proposed", "N_Simple"), "Total Sample Size (N)"),
  to_long(f2, c("power_proposed", "power_Simple"), "Power"),
  to_long(f2, c("type1_proposed", "type1_Simple"), "Type I Error Rate")
) %>%
  dplyr::mutate(
    omega_factor = factor(omega),
    metric = factor(metric, levels = METRIC_LEVELS)
  )

case_levels2 <- unique(fig2_long$case)
case_labels2 <- vapply(case_levels2, function(label) {
  parts <- strsplit(label, "_")[[1]]
  deparse(bquote(Case ~ .(parts[2]) ~ "(" * pi[1] == .(as.numeric(parts[4])) * ")"))
}, character(1))
fig2_long$case <- factor(fig2_long$case,
                         levels = case_levels2, labels = case_labels2)

omega_breaks <- sort(unique(f2$omega))
omega_labels <- lapply(omega_breaks, function(x) bquote(omega == .(x)))

rho_abs_max <- max(abs(range(fig2_long$rho)))
rho_limit <- min(ceiling(rho_abs_max * 5) / 5, 1.0)

# The feasible correlation range differs by (case, omega); annotate it inside
# the type I error panels so the reader can see how far each curve extends
rho_ranges <- fig2_long %>%
  dplyr::filter(metric == "Type I Error Rate") %>%
  dplyr::group_by(case, omega) %>%
  dplyr::summarise(rho_min = min(rho), rho_max = max(rho),
                   .groups = "drop") %>%
  dplyr::mutate(
    label = paste0("[", sprintf("%.3f", rho_min), ", ",
                   sprintf("%.3f", rho_max), "]"),
    metric = factor("Type I Error Rate", levels = METRIC_LEVELS)
  ) %>%
  dplyr::group_by(case) %>%
  dplyr::mutate(y_pos = 0.0202 + (dplyr::row_number() - 1) * 0.0008) %>%
  dplyr::ungroup()

rho_range_header <- rho_ranges %>%
  dplyr::group_by(case) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    label = "rho~range",
    y_pos = max(rho_ranges$y_pos) + 0.0008,
    omega = NA_real_
  )

fig2 <- ggplot(fig2_long, aes(x = rho, y = value,
                              color = omega_factor,
                              linetype = method)) +
  geom_line(linewidth = 1.5) +
  geom_hline(
    data = dplyr::filter(fig2_long, metric == "Power"),
    aes(yintercept = target_power),
    color = "gray", linetype = "longdash", linewidth = 1
  ) +
  geom_hline(
    data = dplyr::filter(fig2_long, metric == "Type I Error Rate"),
    aes(yintercept = alpha),
    color = "gray", linetype = "longdash", linewidth = 1
  ) +
  geom_text(
    data = rho_ranges,
    aes(x = -rho_limit * 0.95, y = y_pos, label = label,
        color = factor(omega)),
    hjust = 0, vjust = 0, size = 7, fontface = "bold",
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  geom_text(
    data = rho_range_header,
    aes(x = -rho_limit * 0.95, y = y_pos, label = label),
    hjust = 0, vjust = 0, size = 7, fontface = "bold",
    color = "black", parse = TRUE,
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  facet_grid(
    rows = vars(metric),
    cols = vars(case),
    scales = "free_y",
    labeller = labeller(case = label_parsed)
  ) +
  ggh4x::scale_y_facet(metric == "Total Sample Size (N)", limits = NULL) +
  ggh4x::scale_y_facet(metric == "Power", limits = c(0.7, 0.9)) +
  ggh4x::scale_y_facet(metric == "Type I Error Rate", limits = c(0.02, 0.03)) +
  scale_linetype_manual(
    values = c("Proposed" = "solid", "Simple inflation" = "dashed")
  ) +
  scale_color_manual(
    values = OKABE[seq_along(omega_breaks)],
    breaks = omega_breaks,
    labels = omega_labels
  ) +
  scale_x_continuous(
    breaks = round(seq(-rho_limit, rho_limit, by = 0.2), 1),
    limits = c(-rho_limit, rho_limit)
  ) +
  labs(
    x = expression(paste("Correlation (", rho, ")")),
    y = NULL,
    color = "Dropout Probability",
    linetype = "Method"
  ) +
  theme_panel()

save_figure(fig2, "Figure2")

# =============================================================================
# Figure 3: relative efficiency over the (pi1, omega) plane
# =============================================================================

f3 <- fig3_data$df

contour_step <- fig3_data$contour_step
re_min <- floor(min(f3$re, na.rm = TRUE) / contour_step) * contour_step
re_max <- ceiling(max(f3$re, na.rm = TRUE) / contour_step) * contour_step
re_breaks <- seq(re_min, re_max, by = contour_step)

# White at low RE, through yellow and orange, to vermilion at RE = 1 and dark
# red above it
RE_PALETTE <- c("#FFFFFF", "#F0E442", "#E69F00", "#D55E00", "#7B2D00")

f3$delta_label <- factor(f3$delta_label,
                         levels = paste0("delta[Full] == ", fig3_data$delta))

fig3 <- ggplot(f3, aes(x = omega, y = pi1, z = re)) +
  metR::geom_contour_fill(breaks = re_breaks, na.fill = FALSE) +
  facet_wrap(vars(delta_label), ncol = 2, labeller = label_parsed) +
  scale_fill_stepsn(
    colors = RE_PALETTE,
    breaks = re_breaks,
    limits = c(re_min, re_max),
    na.value = "gray80",
    name = expression(RE == N[CC] / N[NRI]),
    guide = guide_colorsteps(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 50,
      barheight = 2,
      even.steps = TRUE
    )
  ) +
  scale_x_continuous(breaks = seq(0, max(f3$omega), by = 0.10),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), expand = c(0, 0)) +
  labs(
    x = expression(paste("Dropout Probability (", omega, ")")),
    y = expression(paste("Response Probability (", pi[1], ")"))
  ) +
  theme_panel() +
  theme(
    plot.margin = margin(10, 30, 10, 10),
    panel.grid = element_blank()
  )

save_figure(fig3, "Figure3")

# =============================================================================
# Table 2: application to the cytisine smoking cessation trial
# =============================================================================

t2 <- table2_data$df

# formatC keeps the sign of a value that rounds to zero, which prints as "-0.000"
zap_neg_zero <- function(x, digits) ifelse(abs(x) < 0.5 * 10 ^ (-digits), 0, x)
fmt_p <- function(x) formatC(zap_neg_zero(x, 4), digits = 4, format = "f")
fmt_i <- function(x) formatC(x, format = "d", big.mark = ",")
fmt_re <- function(x) formatC(zap_neg_zero(x, 2), digits = 2, format = "f")
fmt_g <- function(x) formatC(zap_neg_zero(x, 3), digits = 3, format = "f")

rho_cell <- function(label, value) {
  if (label == "L") {
    sprintf("$L$ (%.3f)", value)
  } else if (label == "U") {
    sprintf("$U$ (%.3f)", value)
  } else {
    "$0$"
  }
}

# Measured against the manuscript class (article, 12pt, a4paper: textwidth
# 390.0pt): at the default column separation this table is 389.5pt wide, which
# fits by half a point and would overflow on any later edit. At 4pt separation
# it is 353.5pt, so it stays at the body font size with room to spare.
tab2_lines <- c(
  "\\begin{table}[htbp]",
  "  \\centering",
  "  \\setlength{\\tabcolsep}{4pt}",
  paste0(
    "  \\caption{Required sample sizes for the cytisine trial under the design ",
    "assumptions of the original trial: $\\pi_{1} = ",
    formatC(table2_data$pi1, digits = 2, format = "f"),
    "$, $\\pi_{0} = ",
    formatC(table2_data$pi0, digits = 2, format = "f"),
    "$, $\\omega_{1} = \\omega_{0} = ",
    formatC(table2_data$omega, digits = 2, format = "f"),
    "$, $r = ", table2_data$r,
    "$, $\\alpha = ", table2_data$alpha,
    "$ (one-sided) and target power $= ", table2_data$target_power,
    "$. The simple inflation method gives $N_{\\mathrm{CC}} = ",
    fmt_i(t2$N_cc[1]),
    "$, which is the sample size the trial enrolled.}"
  ),
  "  \\label{tab:application}",
  "  \\begin{threeparttable}",
  "  \\begin{tabular}{lccccccrr}",
  "    \\toprule",
  paste0(
    "    $\\rho$ & $\\gamma_{1}$ & $\\gamma_{0}$ & ",
    "$\\pi_{1,\\mathrm{NRI}}$ & $\\pi_{0,\\mathrm{NRI}}$ & ",
    "$\\delta_{\\mathrm{NRI}}$ & $N_{\\mathrm{CC}}$ & ",
    "$N_{\\mathrm{NRI}}$ & RE \\\\"
  ),
  "    \\midrule"
)

for (k in seq_len(nrow(t2))) {
  tab2_lines <- c(tab2_lines, paste0(
    "    ",
    rho_cell(t2$rho_label[k], t2$rho[k]), " & ",
    fmt_g(t2$gamma1[k]), " & ",
    fmt_g(t2$gamma0[k]), " & ",
    fmt_p(t2$pi1_NRI[k]), " & ",
    fmt_p(t2$pi0_NRI[k]), " & ",
    fmt_p(t2$effect_NRI[k]), " & ",
    fmt_i(t2$N_cc[k]), " & ",
    fmt_i(t2$N_nri[k]), " & ",
    fmt_re(t2$RE[k]),
    " \\\\"
  ))
}

tab2_lines <- c(
  tab2_lines,
  "    \\bottomrule",
  "  \\end{tabular}",
  "  \\begin{tablenotes}[flushleft]\\footnotesize",
  paste0(
    "    \\item $L$ and $U$ denote the lower and upper attainable values of ",
    "$\\rho$ given the marginal probabilities. ",
    "$\\gamma_{j} = \\Pr(R_{ij} = 1 \\mid D_{ij} = 1)$ is the probability ",
    "that a dropout in group $j$ would have responded, which is the same ",
    "assumption as $\\rho$ expressed on the probability scale; a common ",
    "$\\rho$ implies different $\\gamma_{j}$ in the two groups because the ",
    "response probabilities differ. ",
    "$\\pi_{j,\\mathrm{NRI}} = \\Pr(R_{ij} = 1, D_{ij} = 0)$ is the response ",
    "probability an NRI analysis observes, and ",
    "$\\delta_{\\mathrm{NRI}} = \\pi_{1,\\mathrm{NRI}} - \\pi_{0,\\mathrm{NRI}}$ ",
    "is the effect it is powered against; the corresponding full-data effect ",
    "is $\\delta_{\\mathrm{Full}} = ",
    formatC(t2$effect_full[1], digits = 2, format = "f"),
    "$. RE $= N_{\\mathrm{CC}} / N_{\\mathrm{NRI}}$."
  ),
  "  \\end{tablenotes}",
  "  \\end{threeparttable}",
  "\\end{table}"
)

write_tex(tab2_lines, file.path(results_dir, "table2_application.tex"))

message("All figures and tables written to '", results_dir, "/'.")

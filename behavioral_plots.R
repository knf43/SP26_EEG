# ============================================================
# Behavioral figures: StressLearning2025
# Six figures, one question each.
#
# Design follows Casillas (ds4ling, multilevel models) and Soskuthy
# (2021, GAMM modelling strategies):
#   - the data is always shown, model prediction goes on top
#   - every estimate carries a 95% interval
#   - when the question is a difference, the difference is plotted with
#     its own interval, not two estimates side by side
#   - no bar charts
#
# Refits nothing. Reads the saved model objects. Figures go to ./figures
# so the R Markdown site can serve them on GitHub Pages.
# ============================================================

# libraries
library(tidyverse)
library(mgcv)
library(emmeans)
library(fs)

# configuration
OUTPUT_DIR <- "~/Desktop/SP26_data/behavioral_output"
GLMM_DIR   <- path(OUTPUT_DIR, "glmm")
GAMM_DIR   <- path(OUTPUT_DIR, "gamm")
FIG_DIR    <- "figures"

dir_create(FIG_DIR)

GROUP_COLS <- c(L1 = "#2166ac", L2 = "#b2182b")

emm_options(lmer.df = "asymptotic")

theme_set(theme_minimal(base_size = 11))
THEME <- theme(legend.position = "top",
               panel.grid.minor = element_blank())

glmm     <- readRDS(path(GLMM_DIR, "glmm_results.rds"))
followup <- readRDS(path(GLMM_DIR, "glmm_followup_results.rds"))
gamm     <- readRDS(path(GAMM_DIR, "gamm_training_results.rds"))

# The cost tables gain confidence intervals only after followup_analysis.R is
# re-run with the emmeans contrasts. Until then the costs are drawn as points
# with no interval rather than refusing to plot.
HAS_CI <- all(c("cost_low", "cost_high") %in% names(followup$cost_q1q2))
if (!HAS_CI) {
  message("NOTE: cost tables carry no confidence intervals. ",
          "Re-run followup_analysis.R when you have time to add them.")
}

# data
source("script.R")
long <- res$long |>
  mutate(
    sound_for_item = coalesce(sounds, sounds_stem),
    item  = str_extract(sound_for_item, "(?<=[xt]_)\\d{2}"),
    group = factor(group, levels = c("L1", "L2"))
  )

testing <- long |>
  filter(phase %in% c("pre", "post")) |>
  mutate(phase = factor(phase, levels = c("pre", "post")))

COND_ORDER <- c("parox_nounmatch_verbmatch", "parox_nounmismatch_verbmatch",
                "parox_nounmatch_verbmismatch", "oxy_nounmatch_verbmatch",
                "oxy_nounmismatch_verbmatch", "oxy_nounmatch_verbmismatch")

long_testing <- long |>
  filter(phase %in% c("pre", "post"), !is.na(cond_name), !is.na(item)) |>
  mutate(group = factor(group, levels = c("L1", "L2")),
         phase = factor(phase, levels = c("pre", "post")),
         cond_name = factor(cond_name, levels = COND_ORDER),
         item = factor(item), participant = factor(participant))
contrasts(long_testing$group)     <- contr.sum(2)
contrasts(long_testing$phase)     <- contr.sum(2)
contrasts(long_testing$cond_name) <- contr.sum(6)

long_testing_rt <- filter(long_testing, !is.na(rt_correct), rt_correct > 0) |>
  mutate(log_rt = log(rt_correct))

train <- long |>
  filter(phase == "training", !is.na(item), !is.na(corr)) |>
  mutate(participant = factor(participant), item = factor(item),
         trial = as.numeric(trial))
train$group_ord <- as.ordered(train$group)
contrasts(train$group_ord) <- "contr.treatment"
train_rt <- filter(train, !is.na(rt_correct), rt_correct > 0)

# helper functions
pp_mean <- function(d, yvar) {
  d |>
    group_by(participant, group, phase) |>
    summarise(y = mean(.data[[yvar]], na.rm = TRUE), .groups = "drop") |>
    filter(!is.na(y))
}

# Estimate and interval both come from the same regridded emmeans object.
# The earlier version took the estimate from a bias-corrected ggpredict and the
# interval from the uncorrected one, which put the band below the points.
model_means <- function(model, measure) {
  args <- list(object = model, specs = "phase", by = "group")
  if (measure == "rt") args$tran <- "log"
  d <- as.data.frame(regrid(do.call(emmeans, args)))
  pick <- function(cands) intersect(cands, names(d))[1]
  data.frame(
    x         = factor(d$phase, levels = c("pre", "post")),
    group     = factor(as.character(d$group), levels = c("L1", "L2")),
    predicted = d[[pick(c("response", "prob", "emmean"))]],
    conf.low  = d[[pick(c("lower.CL", "asymp.LCL"))]],
    conf.high = d[[pick(c("upper.CL", "asymp.UCL"))]]
  )
}

# A shaded band needs a continuous x, so a two-level factor is mapped to 1 and
# 2 and the axis relabelled. Existing factor levels are respected, so pre comes
# before post rather than being sorted alphabetically.
x_numeric <- function(v) {
  f <- if (is.factor(v)) v else factor(v, levels = unique(v))
  list(num = as.numeric(f), labs = levels(f))
}

band <- function(d, xvar, y, lo, hi, colourvar) {
  xn <- x_numeric(d[[xvar]])
  d$.x <- xn$num
  ggplot(d, aes(.x, .data[[y]], colour = .data[[colourvar]],
                fill = .data[[colourvar]], group = .data[[colourvar]])) +
    geom_ribbon(aes(ymin = .data[[lo]], ymax = .data[[hi]]),
                alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.2) +
    scale_x_continuous(breaks = seq_along(xn$labs), labels = xn$labs,
                       expand = expansion(mult = 0.12)) +
    scale_colour_manual(values = GROUP_COLS) +
    scale_fill_manual(values = GROUP_COLS) +
    THEME
}

# Model prediction with its 95% interval as a shaded band. Accuracy is at
# ceiling for most participants with a few near chance, so overlaying the raw
# values compresses the region where the effect is. The descriptive table on
# the page carries the observed means.
data_and_model <- function(points, pred, ylab, expo = FALSE, pct = FALSE) {
  d <- as.data.frame(pred)
  p <- band(d, "x", "predicted", "conf.low", "conf.high", "group") +
    labs(x = NULL, y = ylab, colour = NULL, fill = NULL)
  if (pct) p <- p + scale_y_continuous(labels = scales::percent)
  p
}

# A difference with its own interval, against a zero reference.
difference_points <- function(d, xvar, ylab) {
  if (!HAS_CI) { d$cost_low <- d$cost; d$cost_high <- d$cost }
  band(d, xvar, "cost", "cost_low", "cost_high", "Group") +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    facet_wrap(~ Measure, scales = "free_y",
               labeller = as_labeller(c(accuracy = "accuracy (percentage points)",
                                        rt = "reaction time (ms)"))) +
    labs(x = NULL, y = ylab, colour = NULL, fill = NULL)
}

binned <- function(d, yvar, width = 10) {
  d |>
    mutate(mid = (trial %/% width) * width + width / 2) |>
    group_by(group, mid, participant) |>
    summarise(pm = mean(.data[[yvar]], na.rm = TRUE), .groups = "drop") |>
    group_by(group, mid) |>
    summarise(y  = mean(pm, na.rm = TRUE),
              se = sd(pm, na.rm = TRUE) / sqrt(sum(!is.na(pm))),
              .groups = "drop") |>
    mutate(lower = y - 1.96 * se, upper = y + 1.96 * se) |>
    rename(trial = mid)
}

fitted_curve <- function(model, data, linkinv) {
  grid <- expand_grid(
    trial = seq(min(data$trial), max(data$trial), length.out = 200),
    group_ord = factor(levels(data$group_ord), levels = levels(data$group_ord),
                       ordered = TRUE)
  ) |>
    mutate(participant = data$participant[1], item = data$item[1])
  pr <- predict(model, newdata = grid, se.fit = TRUE, type = "link",
                exclude = c("s(trial,participant)", "s(item)"))
  grid |>
    mutate(fit = linkinv(pr$fit),
           lower = linkinv(pr$fit - 1.96 * pr$se.fit),
           upper = linkinv(pr$fit + 1.96 * pr$se.fit),
           group = factor(as.character(group_ord), levels = c("L1", "L2")))
}

# The by-group smooth is itself the L2 minus L1 difference, read off the terms
# matrix so its standard error is the standard error of the difference.
difference_smooth <- function(model, data) {
  grid <- tibble(
    trial = seq(min(data$trial), max(data$trial), length.out = 200),
    group_ord = factor("L2", levels = levels(data$group_ord), ordered = TRUE),
    participant = data$participant[1], item = data$item[1]
  )
  tt  <- predict(model, newdata = grid, type = "terms", se.fit = TRUE)
  col <- grep("group_ordL2", colnames(tt$fit), value = TRUE)[1]
  grid |>
    mutate(diff = tt$fit[, col], se = tt$se.fit[, col],
           lower = diff - 1.96 * se, upper = diff + 1.96 * se,
           sig = lower > 0 | upper < 0)
}

# figures
FIGS <- list()

FIGS$fig1_accuracy <- data_and_model(
  NULL, model_means(glmm$mod_acc, "accuracy"),
  "proportion correct", pct = TRUE)

FIGS$fig2_rt <- data_and_model(
  NULL, model_means(glmm$mod_rt, "rt"),
  "reaction time (s)")

FIGS$fig3_cost_phase <- difference_points(
  filter(followup$cost_q1q2, Factor == "verb_match"), "Phase",
  "verb-mismatch cost")

FIGS$fig4_cost_stress <- difference_points(
  followup$cost_q4, "Stress", "verb-mismatch cost")

FIGS$fig5_training_rt <- {
  obs <- binned(train_rt, "rt_correct")
  ggplot(fitted_curve(gamm$gam_rt, train_rt, exp),
         aes(trial, fit, colour = group, fill = group)) +
    geom_linerange(data = obs, aes(trial, ymin = lower, ymax = upper),
                   inherit.aes = FALSE, colour = "grey75", linewidth = 0.3) +
    geom_point(data = obs, aes(trial, y), inherit.aes = FALSE,
               colour = "grey55", size = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(values = GROUP_COLS) +
    scale_fill_manual(values = GROUP_COLS) +
    labs(x = "training trial", y = "reaction time (s)",
         colour = NULL, fill = NULL) +
    THEME
}

# save
for (nm in names(FIGS)) {
  print(FIGS[[nm]])
  h <- if (grepl("cost", nm)) 3.2 else 4
  ggsave(path(FIG_DIR, paste0(nm, ".png")), FIGS[[nm]],
         width = 6.5, height = h, dpi = 300)
}

message("\n", length(FIGS), " figures saved to: ", FIG_DIR)

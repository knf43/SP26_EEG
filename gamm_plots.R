# ============================================================
# Training GAMM figures: StressLearning2025
# Reads gamm_training_results.rds, refits nothing.
#
# Two figures per measure:
#   1. fitted curves by group, population level, with binned raw data
#   2. the L2 minus L1 difference smooth, with the regions where its
#      95% CI excludes zero marked
#
# The difference smooth is the group contrast the model actually tests,
# so it is plotted rather than left to be read off two overlaid curves.
# ============================================================

# libraries
library(tidyverse)
library(mgcv)
library(fs)

# configuration
OUTPUT_DIR <- "~/Desktop/SP26_data/behavioral_output"
GAMM_DIR   <- path(OUTPUT_DIR, "gamm")

GROUP_COLS <- c(L1 = "#2166ac", L2 = "#b2182b")

theme_set(theme_minimal(base_size = 11))
THEME <- theme(legend.position = "top")

fit <- readRDS(path(GAMM_DIR, "gamm_training_results.rds"))

# data
source("script.R")

train <- res$long |>
  mutate(
    sound_for_item = coalesce(sounds, sounds_stem),
    item = str_extract(sound_for_item, "(?<=[xt]_)\\d{2}")
  ) |>
  filter(phase == "training", !is.na(item), !is.na(corr)) |>
  mutate(
    group       = factor(group, levels = c("L1", "L2")),
    participant = factor(participant),
    item        = factor(item),
    trial       = as.numeric(trial)
  )
train$group_ord <- as.ordered(train$group)
contrasts(train$group_ord) <- "contr.treatment"

train_rt <- filter(train, !is.na(rt_correct), rt_correct > 0)

# helper functions
# Population-level fit, with the by-participant and by-item terms excluded.
fitted_curve <- function(model, data, linkinv) {
  grid <- expand_grid(
    trial     = seq(min(data$trial), max(data$trial), length.out = 200),
    group_ord = factor(levels(data$group_ord), levels = levels(data$group_ord),
                       ordered = TRUE)
  ) |>
    mutate(participant = data$participant[1], item = data$item[1])

  pr <- predict(model, newdata = grid, se.fit = TRUE, type = "link",
                exclude = c("s(trial,participant)", "s(item)"))
  grid |>
    mutate(fit   = linkinv(pr$fit),
           lower = linkinv(pr$fit - 1.96 * pr$se.fit),
           upper = linkinv(pr$fit + 1.96 * pr$se.fit),
           group = factor(as.character(group_ord), levels = c("L1", "L2")))
}

# The by-group smooth term is itself the L2 minus L1 difference, so it is read
# straight off the terms matrix rather than reconstructed by subtraction.
difference_smooth <- function(model, data) {
  grid <- tibble(
    trial       = seq(min(data$trial), max(data$trial), length.out = 200),
    group_ord   = factor("L2", levels = levels(data$group_ord), ordered = TRUE),
    participant = data$participant[1],
    item        = data$item[1]
  )
  tt  <- predict(model, newdata = grid, type = "terms", se.fit = TRUE)
  col <- grep("group_ordL2", colnames(tt$fit), value = TRUE)[1]
  grid |>
    mutate(diff  = tt$fit[, col],
           se    = tt$se.fit[, col],
           lower = diff - 1.96 * se,
           upper = diff + 1.96 * se,
           sig   = lower > 0 | upper < 0)
}

binned <- function(data, yvar, width = 10) {
  data |>
    mutate(bin = (trial %/% width) * width + width / 2) |>
    group_by(group, bin) |>
    summarise(y = mean(.data[[yvar]], na.rm = TRUE), .groups = "drop")
}

plot_fitted <- function(curve, points, ylab, title, hline = NA) {
  p <- ggplot(curve, aes(trial, fit, colour = group, fill = group))
  if (!is.na(hline)) {
    p <- p + geom_hline(yintercept = hline, colour = "grey60", linewidth = 0.3)
  }
  p +
    geom_point(data = points, aes(trial, y), inherit.aes = FALSE,
               colour = "grey55", size = 0.9, alpha = 0.7) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.7) +
    scale_colour_manual(values = GROUP_COLS) +
    scale_fill_manual(values = GROUP_COLS) +
    labs(title = title, x = "training trial", y = ylab,
         colour = NULL, fill = NULL) +
    THEME
}

plot_difference <- function(d, ylab, title) {
  runs <- d |>
    mutate(run = cumsum(sig != lag(sig, default = first(sig)))) |>
    filter(sig) |>
    group_by(run) |>
    summarise(xmin = min(trial), xmax = max(trial), .groups = "drop")

  p <- ggplot(d, aes(trial, diff))
  if (nrow(runs) > 0) {
    p <- p + geom_rect(data = runs, inherit.aes = FALSE,
                       aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                       fill = "grey85", alpha = 0.6)
  }
  p +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#4d4d4d", alpha = 0.18) +
    geom_line(linewidth = 0.7, colour = "#4d4d4d") +
    labs(title = title, x = "training trial", y = ylab,
         caption = if (nrow(runs) > 0)
           "shaded: 95% CI excludes zero" else
           "95% CI includes zero throughout") +
    THEME
}

# figures
curve_acc <- fitted_curve(fit$gam_acc, train,    plogis)
curve_rt  <- fitted_curve(fit$gam_rt,  train_rt, exp)
diff_acc  <- difference_smooth(fit$gam_acc, train)
diff_rt   <- difference_smooth(fit$gam_rt,  train_rt)

p1 <- plot_fitted(curve_acc, binned(train, "corr"),
                  "proportion correct", "training accuracy", hline = 0.5)
p2 <- plot_fitted(curve_rt, binned(train_rt, "rt_correct"),
                  "rt (s)", "training rt, correct trials")
p3 <- plot_difference(diff_acc, "L2 minus L1 (log odds)",
                      "training accuracy, group difference smooth")
p4 <- plot_difference(diff_rt, "L2 minus L1 (log rt)",
                      "training rt, group difference smooth")

figs <- list(gamm_training_accuracy = p1,
             gamm_training_rt       = p2,
             gamm_diff_accuracy     = p3,
             gamm_diff_rt           = p4)

for (nm in names(figs)) {
  print(figs[[nm]])
  ggsave(path(GAMM_DIR, paste0(nm, ".png")), figs[[nm]],
         width = 6.5, height = 4, dpi = 300)
}

message("\nFigures saved to: ", GAMM_DIR)

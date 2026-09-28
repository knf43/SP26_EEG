# ============================================================
# Training-phase GAMMs: StressLearning2025
# Non-linear change in accuracy and RT across the 120 training trials
#
# Pre-specified in the analysis plan as the alternative to the linear
# growth curve in statistical_analysis.R, to be reported when the linear
# trial slope is a poor description of the learning curve.
#
# group_ord is an ordered factor, so s(trial) is the L1 curve and
# s(trial, by = group_ord) is the L2 departure from it. A significant
# difference smooth is the group by learning effect, tested directly
# rather than by comparing two curves by eye.
# ============================================================

# libraries
library(tidyverse)
library(mgcv)
library(fs)

# configuration
OUTPUT_DIR <- "~/Desktop/SP26_data/behavioral_output"
GAMM_DIR   <- path(OUTPUT_DIR, "gamm")
K_TRIAL    <- 8        # 120 trials per participant; larger k chases noise

dir_create(GAMM_DIR)
theme_set(theme_minimal(base_size = 12))

# data preparation
source("script.R")
long <- res$long

train <- long |>
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

message("Training data: ", nrow(train), " accuracy trials, ",
        nrow(train_rt), " RT trials, ",
        n_distinct(train$participant), " participants, ",
        n_distinct(train$item), " items")

# models
message("\n=== training accuracy GAMM ===")
gam_acc <- bam(
  corr ~ group_ord +
    s(trial, k = K_TRIAL) +
    s(trial, by = group_ord, k = K_TRIAL) +
    s(trial, participant, bs = "fs", m = 1, k = 5) +
    s(item, bs = "re"),
  data = train, family = binomial, method = "fREML", discrete = TRUE
)

message("\n=== training RT GAMM ===")
gam_rt <- bam(
  rt_correct ~ group_ord +
    s(trial, k = K_TRIAL) +
    s(trial, by = group_ord, k = K_TRIAL) +
    s(trial, participant, bs = "fs", m = 1, k = 5) +
    s(item, bs = "re"),
  data = train_rt, family = Gamma(link = "log"),
  method = "fREML", discrete = TRUE
)

# model checks
# edf near 1 means the curve is effectively a straight line, and the linear
# growth curve in statistical_analysis.R is then the right model to report.
for (nm in c("gam_acc", "gam_rt")) {
  cat("\n\n========== ", nm, " ==========\n")
  print(summary(get(nm)))
  cat("\n--- gam.check ---\n")
  print(gam.check(get(nm)))
}

edf_table <- bind_rows(
  as.data.frame(summary(gam_acc)$s.table) |>
    rownames_to_column("smooth") |> mutate(model = "accuracy"),
  as.data.frame(summary(gam_rt)$s.table) |>
    rownames_to_column("smooth") |> mutate(model = "rt")
) |>
  filter(str_detect(smooth, fixed("s(trial)")) |
         str_detect(smooth, fixed("group_ord"))) |>
  relocate(model, smooth)

cat("\n=== trial smooths: edf near 1 means linear ===\n")
print(edf_table)

# results
# Population-level curves, with the by-participant and by-item terms excluded.
curve_for <- function(model, data, linkinv) {
  grid <- expand_grid(
    trial     = seq(min(data$trial), max(data$trial), length.out = 120),
    group_ord = factor(levels(data$group_ord), levels = levels(data$group_ord),
                       ordered = TRUE)
  ) |>
    mutate(participant = data$participant[1], item = data$item[1])

  pr <- predict(model, newdata = grid, se.fit = TRUE, type = "link",
                exclude = c("s(trial,participant)", "s(item)"))
  grid |>
    mutate(
      fit   = linkinv(pr$fit),
      lower = linkinv(pr$fit - 1.96 * pr$se.fit),
      upper = linkinv(pr$fit + 1.96 * pr$se.fit),
      group = factor(as.character(group_ord), levels = c("L1", "L2"))
    )
}

curve_acc <- curve_for(gam_acc, train,    plogis)
curve_rt  <- curve_for(gam_rt,  train_rt, exp)

# figures
plot_curve <- function(d, ylab, title, hline = NA) {
  p <- ggplot(d, aes(trial, fit, colour = group, fill = group)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.18, colour = NA) +
    geom_line(linewidth = 1) +
    labs(x = "Training trial", y = ylab, colour = "Group", fill = "Group",
         title = title)
  if (!is.na(hline)) {
    p <- p + geom_hline(yintercept = hline, linetype = "dashed",
                        colour = "grey50")
  }
  p
}

p_acc <- plot_curve(curve_acc, "Predicted accuracy",
                    "Training accuracy across trials (GAMM)", hline = 0.5)
p_rt  <- plot_curve(curve_rt, "Predicted RT (s)",
                    "Training RT across trials (GAMM, correct only)")

print(p_acc)
print(p_rt)

ggsave(path(GAMM_DIR, "gamm_training_accuracy.png"), p_acc,
       width = 7, height = 4.5, dpi = 150)
ggsave(path(GAMM_DIR, "gamm_training_rt.png"), p_rt,
       width = 7, height = 4.5, dpi = 150)

# save
saveRDS(list(
  gam_acc   = gam_acc,
  gam_rt    = gam_rt,
  edf_table = edf_table,
  curve_acc = curve_acc,
  curve_rt  = curve_rt
), path(GAMM_DIR, "gamm_training_results.rds"))

message("\nGAMM outputs saved to: ", GAMM_DIR)

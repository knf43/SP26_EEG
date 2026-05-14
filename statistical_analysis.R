# ============================================================
# Statistical analysis: StressLearning2025
# GLMMs for accuracy and RT, growth curve for training
# Following Casillas (ds4ling) framework: lme4, sum contrasts,
# nested model comparisons, DHARMa diagnostics
# ============================================================

library(tidyverse)
library(lme4)
library(lmerTest)        # gives p-values for lmer (Satterthwaite df)
library(DHARMa)          # GLMM residual diagnostics
library(MuMIn)           # marginal & conditional R^2
library(broom.mixed)     # tidy() for mixed models
library(ggeffects)       # model-predicted means
library(openxlsx)
library(fs)

# ---- CONFIG -------------------------------------------------
DATA_DIR   <- "~/Desktop/SP26_data/2_psychopy_data"
OUTPUT_DIR <- "~/Desktop/SP26_data/behavioral_output"
GLMM_DIR   <- path(OUTPUT_DIR, "glmm")
RT_LOWER   <- 0.200
RT_UPPER   <- 5.000

dir_create(GLMM_DIR)

# ---- Re-source the data prep from the descriptive script ----
# We need the long-format trial data. Rather than duplicate the
# extraction code, source script.R (which runs main() and leaves
# `res$long` in the environment).
source("script.R")
long <- res$long

# =============================================================
# SECTION 1: DATA PREP
# =============================================================
# - Set factor levels
# - Apply sum-to-zero contrasts (Casillas-style, for factorial designs)
# - Extract item ID from sound filename
# - Build separate testing and training datasets

# Extract item ID (the "01" in "x_01a_buscar.wav" or "t_01a_bajar_cut.wav")
long <- long |>
  mutate(
    sound_for_item = coalesce(sounds, sounds_stem),
    item = str_extract(sound_for_item, "(?<=[xt]_)\\d{2}")
  )

# Testing data: pre + post phases
long_testing <- long |>
  filter(phase %in% c("pre", "post"), !is.na(cond_name), !is.na(item)) |>
  mutate(
    group     = factor(group, levels = c("L1", "L2")),
    phase     = factor(phase, levels = c("pre", "post")),
    cond_name = factor(cond_name, levels = COND_NAME_ORDER),
    item      = factor(item),
    participant = factor(participant)
  )

# Sum-to-zero contrasts (effect coding)
contrasts(long_testing$group)     <- contr.sum(2)
contrasts(long_testing$phase)     <- contr.sum(2)
contrasts(long_testing$cond_name) <- contr.sum(6)

# Training data
long_training <- long |>
  filter(phase == "training", !is.na(cond_name), !is.na(item)) |>
  mutate(
    group       = factor(group, levels = c("L1", "L2")),
    cond_name   = factor(cond_name, levels = c("parox", "oxy")),
    item        = factor(item),
    participant = factor(participant),
    trial_z     = as.numeric(scale(trial))  # standardize for the growth curve
  )

contrasts(long_training$group)     <- contr.sum(2)
contrasts(long_training$cond_name) <- contr.sum(2)

message("Testing data:  ", nrow(long_testing), " trials, ",
        n_distinct(long_testing$participant), " participants, ",
        n_distinct(long_testing$item), " items")
message("Training data: ", nrow(long_training), " trials, ",
        n_distinct(long_training$participant), " participants, ",
        n_distinct(long_training$item), " items")

# =============================================================
# SECTION 2: MODELS
# =============================================================
# Path A: simplified fixed structure (group*phase + cond_name)
# Random: (1 + phase | participant) + (1 | item)
# Optimizer fallback for convergence issues

ctrl_bin <- glmerControl(optimizer = "bobyqa",
                         optCtrl   = list(maxfun = 2e5))

# --- Model 1: Accuracy (testing) ---
message("\n=== Fitting Model 1: Accuracy (testing) ===")
mod_acc <- glmer(
  corr ~ group * phase + cond_name +
    (1 + phase | participant) + (1 | item),
  data    = long_testing,
  family  = binomial(link = "logit"),
  control = ctrl_bin
)

# --- Model 2: RT (testing, correct only) ---
message("\n=== Fitting Model 2: RT (testing, correct only) ===")
long_testing_rt <- filter(long_testing, !is.na(rt_correct), rt_correct > 0)
mod_rt <- glmer(
  rt_correct ~ group * phase + cond_name +
    (1 + phase | participant) + (1 | item),
  data    = long_testing_rt,
  family  = Gamma(link = "log"),
  control = ctrl_bin
)

# --- Model 3: Training (growth curve) ---
message("\n=== Fitting Model 3: Training growth curve ===")
mod_training <- glmer(
  corr ~ group * trial_z +
    (1 + trial_z | participant) + (1 | item),
  data    = long_training,
  family  = binomial(link = "logit"),
  control = ctrl_bin
)

# =============================================================
# SECTION 3: DIAGNOSTICS
# =============================================================
# DHARMa scaled residuals (the GLMM analog of Casillas's lm() residual checks)

run_diagnostics <- function(mod, label) {
  message("\n--- Diagnostics: ", label, " ---")
  sim <- simulateResiduals(mod, n = 1000)

  # Save DHARMa summary plot (the headline "is the model OK?" figure)
  png(path(GLMM_DIR, paste0("dharma_", label, ".png")),
      width = 1200, height = 600, res = 120)
  plot(sim)
  dev.off()

  # Tests
  ks   <- testUniformity(sim, plot = FALSE)
  disp <- testDispersion(sim, plot = FALSE)

  data.frame(
    model = label,
    KS_pvalue = signif(ks$p.value, 3),
    dispersion = signif(disp$statistic, 3),
    dispersion_p = signif(disp$p.value, 3)
  )
}

diag_acc      <- run_diagnostics(mod_acc, "accuracy")
diag_rt       <- run_diagnostics(mod_rt, "rt")
diag_training <- run_diagnostics(mod_training, "training")

diagnostics_summary <- bind_rows(diag_acc, diag_rt, diag_training)
print(diagnostics_summary)

# =============================================================
# SECTION 4: NESTED MODEL COMPARISONS
# =============================================================
# Drop each fixed-effect term one at a time, compare to full model

nmc_compare <- function(full_model, drop_term, label, data) {
  full_formula <- formula(full_model)
  reduced_formula <- update(full_formula, paste(". ~ . -", drop_term))

  family_info <- family(full_model)
  ctrl <- if (inherits(full_model, "glmerMod")) ctrl_bin else lmerControl()

  if (inherits(full_model, "glmerMod")) {
    reduced_mod <- glmer(reduced_formula, data = data,
                         family = family_info, control = ctrl)
  } else {
    reduced_mod <- lmer(reduced_formula, data = data, control = ctrl)
  }

  cmp <- anova(reduced_mod, full_model)
  data.frame(
    model = label,
    effect = drop_term,
    chisq = round(cmp$Chisq[2], 2),
    df = cmp$`Chi Df`[2],
    p = signif(cmp$`Pr(>Chisq)`[2], 3)
  )
}

nmc_acc <- bind_rows(
  nmc_compare(mod_acc, "group:phase", "accuracy", long_testing),
  nmc_compare(mod_acc, "cond_name",   "accuracy", long_testing)
)
# For accuracy main effects (group, phase) we drop them along with the interaction
# in a two-step process to preserve marginality
mod_acc_nointer <- update(mod_acc, . ~ . - group:phase)
nmc_acc <- bind_rows(
  nmc_acc,
  data.frame(
    model = "accuracy", effect = "group",
    chisq = round(anova(update(mod_acc_nointer, . ~ . - group), mod_acc_nointer)$Chisq[2], 2),
    df    = anova(update(mod_acc_nointer, . ~ . - group), mod_acc_nointer)$`Chi Df`[2],
    p     = signif(anova(update(mod_acc_nointer, . ~ . - group), mod_acc_nointer)$`Pr(>Chisq)`[2], 3)
  ),
  data.frame(
    model = "accuracy", effect = "phase",
    chisq = round(anova(update(mod_acc_nointer, . ~ . - phase), mod_acc_nointer)$Chisq[2], 2),
    df    = anova(update(mod_acc_nointer, . ~ . - phase), mod_acc_nointer)$`Chi Df`[2],
    p     = signif(anova(update(mod_acc_nointer, . ~ . - phase), mod_acc_nointer)$`Pr(>Chisq)`[2], 3)
  )
)

nmc_rt <- bind_rows(
  nmc_compare(mod_rt, "group:phase", "rt", long_testing_rt),
  nmc_compare(mod_rt, "cond_name",   "rt", long_testing_rt)
)
mod_rt_nointer <- update(mod_rt, . ~ . - group:phase)
nmc_rt <- bind_rows(
  nmc_rt,
  data.frame(
    model = "rt", effect = "group",
    chisq = round(anova(update(mod_rt_nointer, . ~ . - group), mod_rt_nointer)$Chisq[2], 2),
    df    = anova(update(mod_rt_nointer, . ~ . - group), mod_rt_nointer)$`Chi Df`[2],
    p     = signif(anova(update(mod_rt_nointer, . ~ . - group), mod_rt_nointer)$`Pr(>Chisq)`[2], 3)
  ),
  data.frame(
    model = "rt", effect = "phase",
    chisq = round(anova(update(mod_rt_nointer, . ~ . - phase), mod_rt_nointer)$Chisq[2], 2),
    df    = anova(update(mod_rt_nointer, . ~ . - phase), mod_rt_nointer)$`Chi Df`[2],
    p     = signif(anova(update(mod_rt_nointer, . ~ . - phase), mod_rt_nointer)$`Pr(>Chisq)`[2], 3)
  )
)

# Training NMC
nmc_train <- bind_rows(
  nmc_compare(mod_training, "group:trial_z", "training", long_training)
)
mod_train_nointer <- update(mod_training, . ~ . - group:trial_z)
nmc_train <- bind_rows(
  nmc_train,
  data.frame(
    model = "training", effect = "group",
    chisq = round(anova(update(mod_train_nointer, . ~ . - group), mod_train_nointer)$Chisq[2], 2),
    df    = anova(update(mod_train_nointer, . ~ . - group), mod_train_nointer)$`Chi Df`[2],
    p     = signif(anova(update(mod_train_nointer, . ~ . - group), mod_train_nointer)$`Pr(>Chisq)`[2], 3)
  ),
  data.frame(
    model = "training", effect = "trial_z",
    chisq = round(anova(update(mod_train_nointer, . ~ . - trial_z), mod_train_nointer)$Chisq[2], 2),
    df    = anova(update(mod_train_nointer, . ~ . - trial_z), mod_train_nointer)$`Chi Df`[2],
    p     = signif(anova(update(mod_train_nointer, . ~ . - trial_z), mod_train_nointer)$`Pr(>Chisq)`[2], 3)
  )
)

nmc_all <- bind_rows(nmc_acc, nmc_rt, nmc_train)
print(nmc_all)

# =============================================================
# SECTION 5: FIXED EFFECTS TABLES (beta, SE, CI, p)
# =============================================================
get_fixed <- function(mod, label) {
  tidy(mod, effects = "fixed", conf.int = TRUE) |>
    mutate(model = label,
           estimate = round(estimate, 3),
           std.error = round(std.error, 3),
           conf.low = round(conf.low, 3),
           conf.high = round(conf.high, 3),
           p.value = signif(p.value, 3)) |>
    select(model, term, estimate, std.error, conf.low, conf.high, p.value)
}

fixed_acc      <- get_fixed(mod_acc, "accuracy")
fixed_rt       <- get_fixed(mod_rt, "rt")
fixed_training <- get_fixed(mod_training, "training")
fixed_all      <- bind_rows(fixed_acc, fixed_rt, fixed_training)
print(fixed_all)

# =============================================================
# SECTION 6: MODEL FIT (marginal & conditional R^2)
# =============================================================
r2_acc      <- r.squaredGLMM(mod_acc)
r2_rt       <- r.squaredGLMM(mod_rt)
r2_training <- r.squaredGLMM(mod_training)

# For binomial, r.squaredGLMM returns a matrix with theoretical & delta versions
r2_summary <- data.frame(
  model         = c("accuracy", "rt", "training"),
  R2_marginal   = c(r2_acc[1, "R2m"],    r2_rt[1, "R2m"],    r2_training[1, "R2m"]),
  R2_conditional= c(r2_acc[1, "R2c"],    r2_rt[1, "R2c"],    r2_training[1, "R2c"])
) |>
  mutate(across(starts_with("R2"), \(x) round(x, 3)))
print(r2_summary)

# =============================================================
# SECTION 7: INTERPRETABLE EFFECTS (probabilities & ms)
# =============================================================
# For accuracy: report population-level predicted probabilities by group x phase
# For RT: report population-level predicted RT in ms

pred_acc <- ggpredict(mod_acc, terms = c("phase", "group"))
pred_rt  <- ggpredict(mod_rt,  terms = c("phase", "group"))
pred_train <- ggpredict(mod_training, terms = c("trial_z [all]", "group"))

interp_acc <- as.data.frame(pred_acc) |>
  rename(phase = x, group = group, predicted_prob = predicted,
         CI_low = conf.low, CI_high = conf.high) |>
  mutate(predicted_pct = sprintf("%.1f%%", 100 * predicted_prob),
         CI_pct = sprintf("[%.1f%%, %.1f%%]", 100 * CI_low, 100 * CI_high))

interp_rt <- as.data.frame(pred_rt) |>
  rename(phase = x, group = group, predicted_rt = predicted,
         CI_low = conf.low, CI_high = conf.high) |>
  mutate(predicted_ms = sprintf("%d ms", round(predicted_rt * 1000)),
         CI_ms = sprintf("[%d, %d]",
                         round(CI_low * 1000), round(CI_high * 1000)))

# =============================================================
# SECTION 8: PLOTS
# =============================================================
theme_set(theme_minimal(base_size = 12))

# --- Plot 1: model-predicted means by group x phase ---
p_pred_acc <- ggplot(as.data.frame(pred_acc),
                    aes(x = x, y = predicted, colour = group, group = group)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                  position = position_dodge(0.2), size = 0.8) +
  geom_line(position = position_dodge(0.2)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  coord_cartesian(ylim = c(0, 1.05)) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Model-predicted accuracy by group and phase",
       x = "Phase", y = "Predicted P(correct)",
       caption = "Error bars = 95% CI from GLMM")
ggsave(path(GLMM_DIR, "pred_accuracy.png"), p_pred_acc,
       width = 7, height = 5, dpi = 150)

p_pred_rt <- ggplot(as.data.frame(pred_rt),
                   aes(x = x, y = predicted * 1000,
                       colour = group, group = group)) +
  geom_pointrange(aes(ymin = conf.low * 1000, ymax = conf.high * 1000),
                  position = position_dodge(0.2), size = 0.8) +
  geom_line(position = position_dodge(0.2)) +
  labs(title = "Model-predicted RT by group and phase",
       x = "Phase", y = "Predicted RT (ms)",
       caption = "Error bars = 95% CI from GLMM (Gamma-log)")
ggsave(path(GLMM_DIR, "pred_rt.png"), p_pred_rt,
       width = 7, height = 5, dpi = 150)

# --- Plot 2: sleepstudy-style trajectories (training) ---
# Per-participant predicted curves from the training model
train_pred <- long_training |>
  mutate(fitted = fitted(mod_training))

p_train_traj <- ggplot(train_pred, aes(x = trial, y = fitted,
                                      group = participant, colour = group)) +
  geom_line(alpha = 0.5, linewidth = 0.5) +
  facet_wrap(~ group) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  coord_cartesian(ylim = c(0, 1.05)) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Training trajectories: per-participant model fits",
       x = "Trial number", y = "Predicted P(correct)",
       caption = "Each line = one participant's fitted curve from the growth-curve GLMM")
ggsave(path(GLMM_DIR, "training_trajectories.png"), p_train_traj,
       width = 10, height = 5, dpi = 150)

# --- Plot 3: forest plot of fixed effects ---
forest_data <- fixed_all |>
  filter(term != "(Intercept)") |>
  mutate(sig = case_when(p.value < .001 ~ "p < .001",
                         p.value < .01  ~ "p < .01",
                         p.value < .05  ~ "p < .05",
                         TRUE           ~ "n.s."),
         model = factor(model, levels = c("accuracy", "rt", "training")))

p_forest <- ggplot(forest_data,
                  aes(x = estimate, y = term, colour = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high), size = 0.7) +
  facet_wrap(~ model, scales = "free", ncol = 1) +
  scale_colour_manual(values = c("p < .001" = "#cc0000",
                                  "p < .01"  = "#ee7700",
                                  "p < .05"  = "#ffaa00",
                                  "n.s."     = "grey50")) +
  labs(title = "Fixed-effect estimates (β with 95% CI)",
       x = "Estimate (log-odds for accuracy/training; log-RT for rt)",
       y = NULL, colour = "Significance")
ggsave(path(GLMM_DIR, "forest_plot.png"), p_forest,
       width = 9, height = 11, dpi = 150)

# =============================================================
# SECTION 9: SAVE EVERYTHING
# =============================================================
write.xlsx(
  list(
    fixed_effects        = fixed_all,
    nested_comparisons   = nmc_all,
    model_fit_R2         = r2_summary,
    diagnostics_summary  = diagnostics_summary,
    predicted_accuracy   = interp_acc,
    predicted_rt         = interp_rt
  ),
  path(GLMM_DIR, "glmm_results.xlsx"),
  overwrite = TRUE
)

saveRDS(list(
  mod_acc = mod_acc, mod_rt = mod_rt, mod_training = mod_training,
  fixed_all = fixed_all, nmc_all = nmc_all, r2_summary = r2_summary,
  diagnostics_summary = diagnostics_summary,
  interp_acc = interp_acc, interp_rt = interp_rt,
  pred_acc = pred_acc, pred_rt = pred_rt, pred_train = pred_train
), path(GLMM_DIR, "glmm_results.rds"))

message("\nAll GLMM outputs saved to: ", GLMM_DIR)
message("- glmm_results.xlsx (tables)")
message("- glmm_results.rds (full objects for Rmd)")
message("- pred_accuracy.png, pred_rt.png, training_trajectories.png, forest_plot.png")
message("- dharma_*.png (one per model)")

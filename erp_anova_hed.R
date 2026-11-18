# erp anovas

# libraries
library(tidyverse)
library(fs)
library(afex)
library(lmerTest)
library(emmeans)

# configuration
OUTPUT_DIR <- "~/Desktop/SP26_data/erp_output"
TRIALS_FILE <- "~/Desktop/SP26_data/behavioral_output/initial/behavioral_trials.csv"
ANOVA_DIR  <- path(OUTPUT_DIR, "anova")
BLOCK_SIZE <- 30

dir_create(ANOVA_DIR)

# data
source("erp_script.R")

dat <- erp$erp |>
  mutate(group       = factor(group, levels = c("L1", "L2")),
         verb_stress = factor(verb_stress, levels = c("paroxytone", "oxytone")),
         verb_match  = factor(verb_match, levels = c("Match", "Mismatch")),
         noun_match  = factor(noun_match, levels = c("Match", "Mismatch")),
         splice      = factor(splice, levels = c("a", "c")))

prepost <- dat |>
  filter(erp_phase == "experimental") |>
  mutate(phase = factor(phase, levels = c("exp_pre", "exp_post"),
                        labels = c("pre", "post")))

training_position <- load_audio(AUDIO_DIR) |>
  filter(phase == "training") |>
  group_by(participant) |>
  mutate(position = rank(trial)) |>
  ungroup() |>
  select(participant, trial, position)

training <- dat |>
  filter(erp_phase == "training", component == "PrAN_verb") |>
  left_join(training_position, by = c("participant", "trial")) |>
  mutate(block = factor(paste0("block", ceiling(position / BLOCK_SIZE))))

# models
VERB_CONDITIONS <- c("a", "c", "e", "f")
NOUN_CONDITIONS <- c("a", "b", "d", "f")
ALL_CONDITIONS  <- c("a", "b", "c", "d", "e", "f")

MODELS <- list(
  PrAN_noun = list(conditions = ALL_CONDITIONS,
                   within = c("phase", "splice"),
                   effect = "splice"),
  PrAN_noun_early = list(conditions = ALL_CONDITIONS,
                         within = c("phase", "splice"),
                         effect = "splice"),
  PrAN_noun_early_left = list(conditions = ALL_CONDITIONS,
                              within = c("phase", "splice"),
                              effect = "splice"),
  PrAN_verb = list(conditions = NOUN_CONDITIONS,
                   within = c("phase", "verb_stress", "noun_match"),
                   effect = "verb_stress"),
  PrAN_verb_left = list(conditions = NOUN_CONDITIONS,
                        within = c("phase", "verb_stress", "noun_match"),
                        effect = "verb_stress"),
  LAN_noun  = list(conditions = NOUN_CONDITIONS,
                   within = c("phase", "verb_stress", "noun_match"),
                   effect = "noun_match"),
  N400_noun = list(conditions = NOUN_CONDITIONS,
                   within = c("phase", "verb_stress", "noun_match"),
                   effect = "noun_match"),
  P600_noun = list(conditions = NOUN_CONDITIONS,
                   within = c("phase", "verb_stress", "noun_match"),
                   effect = "noun_match"),
  LAN_verb  = list(conditions = VERB_CONDITIONS,
                   within = c("phase", "verb_stress", "verb_match"),
                   effect = "verb_match"),
  N400_verb = list(conditions = VERB_CONDITIONS,
                   within = c("phase", "verb_stress", "verb_match"),
                   effect = "verb_match"),
  P600_verb = list(conditions = VERB_CONDITIONS,
                   within = c("phase", "verb_stress", "verb_match"),
                   effect = "verb_match")
)

# helper functions
cell_means <- function(data, within) {
  data |>
    group_by(participant, group, across(all_of(within))) |>
    summarise(amplitude = mean(amplitude), trials = n(), .groups = "drop")
}

describe <- function(means, label, effect, time) {
  means |>
    group_by(participant, group, time = .data[[time]], level = .data[[effect]]) |>
    summarise(amplitude = mean(amplitude), .groups = "drop") |>
    group_by(group, time, level) |>
    summarise(mean = round(mean(amplitude), 2), sd = round(sd(amplitude), 2),
              participants = n(), .groups = "drop") |>
    mutate(model = label, time = as.character(time),
           level = as.character(level), .before = 1)
}

fit_anova <- function(means, within) {
  suppressMessages(
    aov_ez(id = "participant", dv = "amplitude", data = means,
           between = "group", within = within,
           anova_table = list(es = "pes"))
  )
}

anova_table <- function(fit, label) {
  as.data.frame(fit$anova_table) |>
    rownames_to_column("effect") |>
    transmute(model = label, effect,
              df1 = round(`num Df`, 2), df2 = round(`den Df`, 2),
              F = round(F, 2), p = round(`Pr(>F)`, 4), pes = round(pes, 3))
}

followup <- function(fit, label, effect, by) {
  contrasts <- emmeans(fit, as.formula(paste("~", effect, "|",
                                             paste(by, collapse = " * "))),
                       model = "multivariate") |>
    pairs(reverse = TRUE)
  interval <- confint(contrasts, adjust = "none")

  summary(contrasts, by = NULL, adjust = "fdr") |>
    as_tibble() |>
    mutate(lower = interval$lower.CL, upper = interval$upper.CL) |>
    mutate(across(c(estimate, SE, t.ratio, lower, upper), \(x) round(x, 2)),
           p.value = round(p.value, 4)) |>
    rename(time = all_of(by[2])) |>
    mutate(model = label, time = as.character(time), .before = 1)
}

# cell means
prepost_means <- imap(MODELS, \(model, name) {
  prepost |>
    filter(component == name, cond %in% model$conditions) |>
    cell_means(model$within)
})

training_means <- cell_means(training, c("block", "verb_stress"))

# descriptive statistics
prepost_descriptives <- imap(MODELS, \(model, name) {
  describe(prepost_means[[name]], name, model$effect, "phase")
})

training_descriptives <- describe(training_means, "PrAN_training",
                                  "verb_stress", "block")

descriptives <- bind_rows(prepost_descriptives) |>
  bind_rows(training_descriptives)

# pre and post anovas
prepost_fits <- imap(MODELS, \(model, name) {
  fit_anova(prepost_means[[name]], model$within)
})

prepost_anova <- imap(prepost_fits, \(fit, name) anova_table(fit, name))

prepost_followup <- imap(prepost_fits, \(fit, name) {
  followup(fit, name, MODELS[[name]]$effect, c("group", "phase"))
})

# training anova
training_fit <- fit_anova(training_means, c("block", "verb_stress"))

training_anova    <- anova_table(training_fit, "PrAN_training")
training_followup <- followup(training_fit, "PrAN_training", "verb_stress",
                              c("group", "block"))

# verb stress by noun match
stress_by_noun <- map_dfr(c("PrAN_verb", "PrAN_verb_left"), \(name) {
  emmeans(prepost_fits[[name]], ~ verb_stress | group * phase * noun_match,
          model = "multivariate") |>
    pairs(reverse = TRUE) |>
    summary(by = NULL, adjust = "fdr") |>
    as_tibble() |>
    transmute(model = name, group, time = as.character(phase), noun_match,
              estimate = round(estimate, 2), se = round(SE, 2),
              t = round(t.ratio, 2), p = round(p.value, 4))
})

# verb stress by noun cue, all six conditions
CUE_COMPONENTS <- c("PrAN_verb", "PrAN_verb_left")

cue_means <- map(set_names(CUE_COMPONENTS), \(name) {
  prepost |>
    filter(component == name) |>
    cell_means(c("phase", "splice", "verb_stress"))
})

cue_fits <- map(cue_means, \(means) {
  fit_anova(means, c("phase", "splice", "verb_stress"))
})

cue_anova <- imap(cue_fits, \(fit, name) anova_table(fit, name)) |>
  bind_rows()

stress_by_cue <- imap(cue_fits, \(fit, name) {
  emmeans(fit, ~ verb_stress | group * phase * splice,
          model = "multivariate") |>
    pairs(reverse = TRUE) |>
    summary(by = NULL, adjust = "fdr") |>
    as_tibble() |>
    transmute(model = name, group, time = as.character(phase),
              noun_cue = if_else(splice == "a", "paroxytone", "oxytone"),
              estimate = round(estimate, 2), se = round(SE, 2),
              t = round(t.ratio, 2), p = round(p.value, 4))
}) |>
  bind_rows()

cue_means_all <- bind_rows(cue_means, .id = "model")

# pre against post for each stress pattern
phase_by_stress <- map_dfr(c("PrAN_verb", "PrAN_verb_left"), \(name) {
  emmeans(prepost_fits[[name]], ~ phase | group * verb_stress,
          model = "multivariate") |>
    pairs(reverse = TRUE) |>
    summary(by = NULL, adjust = "fdr") |>
    as_tibble() |>
    transmute(model = name, group, verb_stress,
              estimate = round(estimate, 2), se = round(SE, 2),
              t = round(t.ratio, 2), p = round(p.value, 4))
})

# late amplitude and reaction time in training
late <- dat |>
  filter(erp_phase == "training", component == "late_verb") |>
  inner_join(read_csv(TRIALS_FILE, show_col_types = FALSE) |>
               filter(phase == "training", !is.na(rt)) |>
               select(participant, sound, rt),
             by = c("participant", "sound")) |>
  mutate(log_rt        = log10(rt) - mean(log10(rt)),
         group_c       = if_else(group == "L2", 0.5, -0.5),
         verb_stress_c = if_else(verb_stress == "oxytone", 0.5, -0.5))

late_means <- cell_means(late, "verb_stress")
late_fit   <- fit_anova(late_means, "verb_stress")
late_anova <- anova_table(late_fit, "late_training")

late_followup <- emmeans(late_fit, ~ verb_stress | group, model = "multivariate") |>
  pairs(reverse = TRUE) |>
  summary(by = NULL, adjust = "fdr") |>
  as_tibble() |>
  transmute(model = "late_training", group,
            estimate = round(estimate, 2), se = round(SE, 2),
            t = round(t.ratio, 2), p = round(p.value, 4))

late_rt_model <- lmer(amplitude ~ group_c * verb_stress_c * log_rt +
                        (1 | participant) + (1 | item),
                      data = late)

late_rt_effects <- as.data.frame(coef(summary(late_rt_model))) |>
  rownames_to_column("term") |>
  transmute(term, estimate = round(Estimate, 3), se = round(`Std. Error`, 3),
            df = round(df, 1), t = round(`t value`, 2), p = round(`Pr(>|t|)`, 4))

late_rt_group <- emtrends(late_rt_model, ~ group_c, var = "log_rt",
                          at = list(group_c = c(-0.5, 0.5)),
                          lmer.df = "asymptotic") |>
  summary(infer = TRUE) |>
  as_tibble() |>
  transmute(group = if_else(group_c < 0, "L1", "L2"),
            estimate = round(log_rt.trend, 2), se = round(SE, 2),
            z = round(z.ratio, 2), p = round(p.value, 4))

late_rt_means <- late |>
  group_by(participant, group, verb_stress) |>
  summarise(rt = mean(rt), amplitude = mean(amplitude), trials = n(),
            .groups = "drop")

# results
anova_all    <- bind_rows(prepost_anova) |> bind_rows(training_anova)
followup_all <- bind_rows(prepost_followup) |> bind_rows(training_followup)
means_all    <- bind_rows(prepost_means, .id = "model") |>
  bind_rows(mutate(training_means, model = "PrAN_training"))

cat("\ndescriptive statistics\n")
print(descriptives, n = Inf)

cat("\nanova table\n")
print(anova_all, row.names = FALSE)

cat("\npairwise comparisons, fdr corrected\n")
print(followup_all, n = Inf)

cat("\nverb stress by noun match, fdr corrected\n")
print(stress_by_noun, n = Inf)

cat("\nverb stress by noun cue, all six conditions, anova\n")
print(cue_anova, row.names = FALSE)

cat("\nverb stress by noun cue, all six conditions, fdr corrected\n")
print(stress_by_cue, n = Inf)

cat("\npost minus pre for each stress pattern, fdr corrected\n")
print(phase_by_stress, n = Inf)

cat("\nlate amplitude in training, 600 to 800 ms\n")
print(late_anova, row.names = FALSE)
print(late_followup)
print(late_rt_effects, row.names = FALSE)
print(late_rt_group, n = Inf)

cat("\neffects at p < .05\n")
print(filter(anova_all, p < .05), row.names = FALSE)

# save
write_csv(descriptives, path(ANOVA_DIR, "anova_descriptives.csv"))
write_csv(anova_all,    path(ANOVA_DIR, "anova_table.csv"))
write_csv(followup_all, path(ANOVA_DIR, "anova_followup.csv"))
write_csv(means_all,    path(ANOVA_DIR, "anova_cell_means.csv"))
write_csv(stress_by_noun, path(ANOVA_DIR, "anova_stress_by_noun_match.csv"))
write_csv(cue_anova,      path(ANOVA_DIR, "anova_noun_cue_table.csv"))
write_csv(stress_by_cue,  path(ANOVA_DIR, "anova_stress_by_noun_cue.csv"))
write_csv(cue_means_all,  path(ANOVA_DIR, "anova_noun_cue_cell_means.csv"))
write_csv(phase_by_stress, path(ANOVA_DIR, "anova_phase_by_stress.csv"))
write_csv(late_anova,      path(ANOVA_DIR, "late_training_anova.csv"))
write_csv(late_followup,   path(ANOVA_DIR, "late_training_followup.csv"))
write_csv(late_rt_effects, path(ANOVA_DIR, "late_training_rt_model.csv"))
write_csv(late_rt_group,   path(ANOVA_DIR, "late_training_rt_by_group.csv"))
write_csv(late_rt_means,   path(ANOVA_DIR, "late_training_means.csv"))

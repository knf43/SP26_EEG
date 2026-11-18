# behavioral analysis

# libraries
library(tidyverse)
library(fs)
library(readxl)
library(tuneR)
library(lme4)
library(lmerTest)
library(emmeans)

# configuration
DATA_DIR      <- "~/Desktop/SP26_data/2_psychopy_data"
AUDIO_DIR     <- "~/Kaylee_SP26/trigger_audio"
SOUND_DIR     <- "~/Kaylee_SP26/sounds"
TIMING_FILE   <- "~/Kaylee_SP26/conditions_stress_training_master.xlsx"
TIMING_SHEET  <- "5_experimental_psychopy_excel"
OUTPUT_DIR    <- "~/Desktop/SP26_data/behavioral_output/initial"
FIGURE_DIR    <- "figures"
RT_MINIMUM    <- 200
RT_MAXIMUM    <- 5000
EXCLUDED_STEM <- "t_52b_tocar_cut.wav"
BLOCK_SIZE    <- 30
GROUP_COLORS  <- c(L1 = "#2166ac", L2 = "#b2182b")

PHASES <- list(
  pre      = c(loop = "exp_pre_loop", response = "resp_e",
               sound = "sound_e", file = "sounds"),
  training = c(loop = "training_loop", response = "training_resp",
               sound = "sound_t_stem", file = "sounds_stem"),
  post     = c(loop = "exp_post_loop", response = "resp_e",
               sound = "sound_e", file = "sounds"))

CUES <- list(verb = list(conditions = c("a", "c", "e", "f"), match = "verb_match"),
             noun = list(conditions = c("a", "b", "d", "f"), match = "noun_match"))

dir_create(OUTPUT_DIR)
dir_create(FIGURE_DIR)
theme_set(theme_minimal(base_size = 11))
emm_options(lmer.df = "asymptotic")

# helper functions
read_phase <- function(raw, phase) {
  columns  <- PHASES[[phase]]
  loop     <- columns[["loop"]]
  response <- paste0(loop, ".", columns[["response"]])
  clock    <- paste0(columns[["response"]], ".started")
  sound    <- paste0(columns[["sound"]], ".started")

  raw |>
    filter(!is.na(.data[[paste0(loop, ".thisN")]])) |>
    transmute(participant  = as.integer(first(na.omit(raw$participant))),
              phase        = phase,
              trial        = as.numeric(.data[[paste0(loop, ".thisN")]]) + 1,
              correct      = as.numeric(.data[[paste0(response, ".corr")]]),
              rt           = suppressWarnings(as.numeric(.data[[paste0(response, ".rt")]])) * 1000,
              clock_offset = (as.numeric(.data[[clock]]) - as.numeric(.data[[sound]])) * 1000,
              sound        = path_file(.data[[columns[["file"]]]]),
              verb_match, noun_match)
}

read_participant <- function(file) {
  raw <- suppressMessages(suppressWarnings(
    read_csv(file, show_col_types = FALSE, guess_max = 5000)))
  map_dfr(names(PHASES), \(phase) read_phase(raw, phase))
}

read_audio_delay <- function(file) {
  read_csv(file, show_col_types = FALSE) |>
    transmute(participant = as.integer(str_extract(path_file(file), "\\d{4}")),
              phase       = str_remove(phase, "exp_"),
              sound,
              audio_delay = used_lag_ms)
}

stem_duration <- function(file) {
  header <- readWave(file, header = TRUE)
  header$samples / header$sample.rate * 1000
}

deviation <- function(x, levels) {
  x <- factor(x, levels = levels)
  contrasts(x) <- matrix(c(-0.5, 0.5), dimnames = list(levels, levels[2]))
  x
}

fixed_effects <- function(model, label) {
  table <- as.data.frame(coef(summary(model)))
  tibble(model     = label,
         term      = rownames(table),
         estimate  = round(table$Estimate, 4),
         se        = round(table$`Std. Error`, 4),
         statistic = round(table[[ncol(table) - 1]], 2),
         p         = round(table[[ncol(table)]], 4))
}

estimated_means <- function(model, label, effect, by) {
  emmeans(model, as.formula(paste("~", effect, "|", paste(by, collapse = " * "))),
          weights = "proportional", type = "response") |>
    as_tibble() |>
    rename(estimate = any_of(c("emmean", "prob")),
           lower = asymp.LCL, upper = asymp.UCL) |>
    select(all_of(c(by, effect)), estimate, lower, upper) |>
    mutate(model = label, .before = 1)
}

followup <- function(model, label, effect, by) {
  emmeans(model, as.formula(paste("~", effect, "|", paste(by, collapse = " * "))),
          weights = "proportional") |>
    pairs(reverse = TRUE) |>
    summary(by = NULL, adjust = "fdr") |>
    as_tibble() |>
    transmute(model = label, across(all_of(by)), contrast,
              estimate = round(estimate, 4), se = round(SE, 4),
              z = round(z.ratio, 2), p = round(p.value, 4))
}

describe <- function(data, by) {
  data |>
    group_by(participant, across(all_of(by))) |>
    summarise(accuracy = mean(correct) * 100,
              rt = mean(rt, na.rm = TRUE), .groups = "drop") |>
    group_by(across(all_of(by))) |>
    summarise(accuracy_mean = round(mean(accuracy), 1),
              accuracy_sd   = round(sd(accuracy), 1),
              rt_mean       = round(mean(rt)),
              rt_sd         = round(sd(rt)),
              participants  = n(), .groups = "drop")
}

block_means <- function(data, variable) {
  data |>
    group_by(participant, group, block, verb_stress) |>
    summarise(value = mean(.data[[variable]], na.rm = TRUE), .groups = "drop") |>
    group_by(group, block, verb_stress) |>
    summarise(estimate = mean(value),
              error    = 1.96 * sd(value) / sqrt(n()), .groups = "drop") |>
    mutate(lower = estimate - error, upper = estimate + error)
}

plot_lines <- function(data, x, line, y_label) {
  dodge <- position_dodge(0.15)
  ggplot(data, aes(.data[[x]], estimate, colour = group,
                   linetype = .data[[line]], group = .data[[line]])) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.12,
                  linetype = "solid", position = dodge) +
    geom_line(linewidth = 0.8, position = dodge) +
    geom_point(size = 2.2, position = dodge) +
    scale_colour_manual(values = GROUP_COLORS, guide = "none") +
    labs(x = NULL, y = y_label, linetype = NULL) +
    theme(legend.position = "top", panel.grid.minor = element_blank())
}

milliseconds <- sec_axis(~ 10^., name = "reaction time (ms)")

# timing
audio_delay <- dir_ls(AUDIO_DIR, glob = "*audio_0*.csv") |>
  map_dfr(read_audio_delay)

suffix_onset <- bind_rows(
  read_excel(path_expand(TIMING_FILE), sheet = TIMING_SHEET) |>
    transmute(sound = path_file(sounds), suffix_onset = trigger_suffix * 1000),
  tibble(file = dir_ls(SOUND_DIR, regexp = "/t_.*_cut\\.wav$")) |>
    transmute(sound = path_file(file), suffix_onset = map_dbl(file, stem_duration)))

# data
long <- dir_ls(DATA_DIR, glob = "*StressLearning2025*.csv") |>
  map_dfr(read_participant) |>
  filter(sound != EXCLUDED_STEM) |>
  left_join(audio_delay, by = c("participant", "phase", "sound")) |>
  left_join(suffix_onset, by = "sound") |>
  group_by(participant) |>
  mutate(audio_delay = coalesce(audio_delay, median(audio_delay, na.rm = TRUE))) |>
  ungroup() |>
  mutate(group  = if_else(participant < 200, "L1", "L2"),
         item   = factor(str_extract(sound, "(?<=[xt]_)\\d{2}")),
         cond   = str_extract(sound, "(?<=_\\d{2})[a-f]"),
         rt     = rt + clock_offset - audio_delay - suffix_onset,
         rt     = if_else(rt >= RT_MINIMUM & rt <= RT_MAXIMUM, rt, NA_real_),
         log_rt = log10(rt),
         participant = factor(participant),
         group  = deviation(group, c("L1", "L2")))

prepost <- long |>
  filter(phase != "training") |>
  mutate(phase    = deviation(phase, c("pre", "post")),
         accuracy = deviation(if_else(correct == 1, "correct", "incorrect"),
                              c("correct", "incorrect")))

cue_data <- imap(CUES, \(cue, name) {
  prepost |>
    filter(cond %in% cue$conditions) |>
    mutate(match = deviation(.data[[cue$match]], c("Match", "Mismatch")))
})

training <- long |>
  filter(phase == "training") |>
  mutate(verb_stress = deviation(if_else(cond == "a", "paroxytone", "oxytone"),
                                 c("paroxytone", "oxytone")),
         block        = ceiling(trial / BLOCK_SIZE),
         trial_scaled = (trial - 1) / (max(trial) - 1))

training_correct <- filter(training, correct == 1)

# descriptive statistics
prepost_descriptives <- imap(cue_data, \(data, name) {
  describe(data, c("group", "phase", "match")) |> mutate(cue = name, .before = 1)
}) |>
  bind_rows()

training_descriptives <- left_join(
  describe(training, c("group", "block", "verb_stress")) |>
    select(group, block, verb_stress, accuracy_mean, accuracy_sd, participants),
  describe(training_correct, c("group", "block", "verb_stress")) |>
    select(group, block, verb_stress, rt_mean, rt_sd),
  by = c("group", "block", "verb_stress"))

# pre and post models
rt_models <- map(cue_data, \(data) {
  lmer(log_rt ~ group * phase * match + group * match * accuracy +
         (1 + phase | participant) + (1 | item),
       data = data, control = lmerControl(optimizer = "bobyqa"))
})

accuracy_models <- map(cue_data, \(data) {
  glmer(correct ~ group * phase * match +
          (1 + phase | participant) + (1 | item),
        data = data, family = binomial,
        control = glmerControl(optimizer = "bobyqa"))
})

# training models
training_rt_model <- lmer(
  log_rt ~ group * trial_scaled * verb_stress +
    (1 + trial_scaled | participant) + (1 | item),
  data = training_correct, control = lmerControl(optimizer = "bobyqa"))

training_accuracy_model <- glmer(
  correct ~ group * trial_scaled * verb_stress +
    (1 + trial_scaled | participant) + (1 | item),
  data = training, family = binomial,
  control = glmerControl(optimizer = "bobyqa"))

# results
MODELS <- c(rt_verb           = rt_models$verb,
            rt_noun           = rt_models$noun,
            accuracy_verb     = accuracy_models$verb,
            accuracy_noun     = accuracy_models$noun,
            rt_training       = training_rt_model,
            accuracy_training = training_accuracy_model)

fixed_all <- imap(MODELS, fixed_effects) |> bind_rows()

followup_all <- imap(MODELS[1:4], \(model, name) {
  followup(model, name, "match", c("group", "phase"))
}) |>
  bind_rows()

means_all <- imap(MODELS[1:4], \(model, name) {
  estimated_means(model, name, "match", c("group", "phase"))
}) |>
  bind_rows() |>
  separate(model, c("measure", "cue")) |>
  mutate(cue = factor(paste(cue, "cue"), levels = c("verb cue", "noun cue")))

cat("\ndescriptive statistics, pre and post\n")
print(prepost_descriptives, n = Inf)

cat("\ndescriptive statistics, training\n")
print(training_descriptives, n = Inf)

cat("\nfixed effects\n")
print(fixed_all, n = Inf)

cat("\nmismatch minus match, fdr corrected\n")
print(followup_all, n = Inf)

cat("\neffects at p < .05\n")
print(filter(fixed_all, p < .05, term != "(Intercept)"), n = Inf)

# save
write_csv(prepost_descriptives,  path(OUTPUT_DIR, "behavioral_descriptives_prepost.csv"))
write_csv(training_descriptives, path(OUTPUT_DIR, "behavioral_descriptives_training.csv"))
write_csv(fixed_all,             path(OUTPUT_DIR, "behavioral_fixed_effects.csv"))
write_csv(followup_all,          path(OUTPUT_DIR, "behavioral_followup.csv"))
write_csv(means_all,             path(OUTPUT_DIR, "behavioral_estimated_means.csv"))
write_csv(select(long, participant, group, phase, trial, sound, correct, rt),
          path(OUTPUT_DIR, "behavioral_trials.csv"))

# plots
FIGURES <- list(
  behavioral_rt = means_all |>
    filter(measure == "rt") |>
    plot_lines("phase", "match", "log-transformed reaction time") +
    facet_grid(cue ~ group) +
    scale_y_continuous(sec.axis = milliseconds),
  behavioral_accuracy = means_all |>
    filter(measure == "accuracy") |>
    plot_lines("phase", "match", "proportion correct") +
    facet_grid(cue ~ group),
  behavioral_training_rt = block_means(training_correct, "log_rt") |>
    plot_lines("block", "verb_stress", "log-transformed reaction time") +
    facet_wrap(~ group) +
    scale_y_continuous(sec.axis = milliseconds) +
    labs(x = "training block"),
  behavioral_training_accuracy = block_means(training, "correct") |>
    plot_lines("block", "verb_stress", "proportion correct") +
    facet_wrap(~ group) +
    labs(x = "training block")
)

for (figure_name in names(FIGURES)) {
  print(FIGURES[[figure_name]])
  ggsave(path(FIGURE_DIR, paste0(figure_name, ".png")), FIGURES[[figure_name]],
         width = 6.5, height = if (grepl("training", figure_name)) 3.2 else 5,
         dpi = 300)
}

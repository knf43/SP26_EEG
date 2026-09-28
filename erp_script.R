# erp data preparation

# libraries
library(tidyverse)
library(fs)

# configuration
BASELINE   <- "event"   # "event" or "sentence"
ERP_FILE   <- if (BASELINE == "event") {
  "~/Kaylee_SP26/sp26_erp_long_event.csv"
} else {
  "~/Kaylee_SP26/sp26_erp_long.csv"
}
AUDIO_DIR  <- "~/Kaylee_SP26/trigger_audio"
OUTPUT_DIR <- "~/Desktop/SP26_data/erp_output"
AMP_LIMIT  <- 100       # uV; NA disables the trial-level amplitude screen

dir_create(OUTPUT_DIR)

# condition legend
COND_LEGEND_PREPOST <- tibble::tribble(
  ~cond, ~cond_name,                       ~verb_stress, ~verb_match, ~suffix_tense, ~noun_match, ~splice,
  "a",   "parox_nounmatch_verbmatch",      "paroxytone", "Match",     "present",     "Match",     "a",
  "b",   "oxy_nounmismatch_verbmatch",     "oxytone",    "Match",     "past",        "Mismatch",  "a",
  "c",   "parox_nounmatch_verbmismatch",   "paroxytone", "Mismatch",  "past",        "Match",     "a",
  "d",   "parox_nounmismatch_verbmatch",   "paroxytone", "Match",     "present",     "Mismatch",  "c",
  "e",   "oxy_nounmatch_verbmismatch",     "oxytone",    "Mismatch",  "present",     "Match",     "c",
  "f",   "oxy_nounmatch_verbmatch",        "oxytone",    "Match",     "past",        "Match",     "c"
)

COND_LEGEND_TRAINING <- tibble::tribble(
  ~cond, ~cond_name, ~verb_stress, ~tense,
  "a",   "parox",    "paroxytone", "present",
  "b",   "oxy",      "oxytone",    "past"
)

COND_NAME_ORDER <- c("parox_nounmatch_verbmatch",
                     "parox_nounmismatch_verbmatch",
                     "parox_nounmatch_verbmismatch",
                     "oxy_nounmatch_verbmatch",
                     "oxy_nounmismatch_verbmatch",
                     "oxy_nounmatch_verbmismatch")

# helper functions
get_group <- function(pid) {
  s <- as.character(as.integer(pid))
  case_when(str_starts(s, "1") ~ "L1",
            str_starts(s, "2") ~ "L2",
            TRUE ~ "Unknown")
}

# the audio csvs carry the pre/post distinction; matlab only knows
# experimental vs training
load_audio <- function(audio_dir) {
  files <- dir_ls(audio_dir, glob = "*audio_0*.csv")
  if (length(files) == 0) stop("no audio csvs in ", audio_dir)
  map_dfr(files, function(f) {
    pid_str <- str_extract(path_file(f), "(?<=audio_)\\d+")
    read_csv(f, show_col_types = FALSE, guess_max = 5000) |>
      transmute(pid_str, participant = as.integer(pid_str),
                trial = as.integer(trial), phase, sound,
                audio_ok = as.integer(ok), audio_corr = as.numeric(corr))
  })
}

load_erp <- function(erp_file) {
  read_csv(erp_file, show_col_types = FALSE,
           col_types = cols(participant = col_character(),
                            group = col_character(), .default = col_guess())) |>
    transmute(pid_str = participant, participant = as.integer(participant),
              group, erp_phase = phase, trial = as.integer(trial),
              event, component, roi, window_start, window_end, amplitude)
}

# main
main <- function(erp_file = ERP_FILE, audio_dir = AUDIO_DIR) {

  erp_raw <- load_erp(erp_file)
  audio   <- load_audio(audio_dir)

  message("baseline: ", BASELINE, "  (", erp_file, ")")
  message("ERP rows: ", nrow(erp_raw),
          " | participants: ", n_distinct(erp_raw$participant),
          " | components: ", n_distinct(erp_raw$component))

  joined <- erp_raw |> left_join(audio, by = c("pid_str", "participant", "trial"))

  n_unmatched <- sum(is.na(joined$sound))
  if (n_unmatched > 0) {
    message("  WARNING: ", n_unmatched, " ERP rows found no audio row")
  } else {
    message("  join: every ERP row matched an audio row")
  }

  mism <- joined |> filter(!is.na(phase)) |>
    filter((erp_phase == "training") != (phase == "training"))
  if (nrow(mism) > 0)
    message("  WARNING: ", nrow(mism), " rows disagree on phase")

  # cond is the letter after the item number; item is the number itself
  joined <- joined |>
    mutate(cond = str_extract(sound, "(?<=_\\d{2})[a-f]"),
           item = str_extract(sound, "(?<=_)\\d{2}"),
           verb = str_extract(sound, "(?<=[a-f]_)[a-z]+(?=\\.wav)"))

  legend_pp <- COND_LEGEND_PREPOST |>
    select(cond, cond_name, verb_stress, verb_match, suffix_tense, noun_match, splice)
  legend_tr <- COND_LEGEND_TRAINING |>
    select(cond, cond_name_tr = cond_name, verb_stress_tr = verb_stress, tense)

  erp <- joined |>
    left_join(legend_pp, by = "cond") |>
    left_join(legend_tr, by = "cond") |>
    mutate(cond_name   = if_else(erp_phase == "training", cond_name_tr, cond_name),
           verb_stress = if_else(erp_phase == "training", verb_stress_tr, verb_stress),
           phase       = factor(phase, levels = c("exp_pre", "training", "exp_post"))) |>
    select(-cond_name_tr, -verb_stress_tr)

  # amplitude screen
  # a pre-sentence baseline can leave a window mean displaced with no local
  # artifact; a per-event baseline mostly removes that, so this should be tiny
  if (!is.na(AMP_LIMIT)) {
    n_before <- nrow(erp)
    over <- erp |> filter(abs(amplitude) > AMP_LIMIT)
    erp  <- erp |> filter(abs(amplitude) <= AMP_LIMIT)
    message("  amplitude screen at +-", AMP_LIMIT, " uV: removed ",
            n_before - nrow(erp), " of ", n_before, " rows (",
            sprintf("%.2f", 100 * (n_before - nrow(erp)) / n_before), "%)")
    if (nrow(over) > 0)
      message("    by participant: ",
              over |> count(participant) |>
                mutate(s = paste0(participant, ":", n)) |> pull(s) |> paste(collapse = " "))
  }

  # summaries
  by_pp <- erp |>
    group_by(participant, group, phase, component) |>
    summarise(n_trials = n(), mean_amp = mean(amplitude, na.rm = TRUE),
              sd_amp = sd(amplitude, na.rm = TRUE), .groups = "drop")

  by_grp <- erp |>
    group_by(participant, group, phase, component, cond_name) |>
    summarise(amp = mean(amplitude, na.rm = TRUE), .groups = "drop") |>
    group_by(group, phase, component, cond_name) |>
    summarise(mean_amp = mean(amp, na.rm = TRUE),
              se_amp = sd(amp, na.rm = TRUE) / sqrt(n()),
              n_participants = n_distinct(participant), .groups = "drop")

  trials_per <- erp |>
    count(participant, group, phase, component, name = "n_trials") |>
    arrange(participant, phase, component)

  message("\n  rows after screening: ", nrow(erp))

  cat("\n--- trials per component, by phase (mean across participants) ---\n")
  trials_per |>
    group_by(phase, component) |>
    summarise(mean_trials = round(mean(n_trials)), min_trials = min(n_trials),
              .groups = "drop") |>
    print(n = 40)

  invisible(list(baseline = BASELINE, erp = erp, by_pp = by_pp, by_grp = by_grp,
                 trials_per = trials_per,
                 cond_legend_prepost = COND_LEGEND_PREPOST,
                 cond_legend_training = COND_LEGEND_TRAINING,
                 amp_limit = AMP_LIMIT))
}

erp <- main()

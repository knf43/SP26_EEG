# erp waveforms at single electrodes

# libraries
library(tidyverse)
library(fs)

# configuration
WAVE_FILE  <- "~/Kaylee_SP26/sp26_waveforms_channels_event.csv"
FIGURE_DIR <- "figures"
COLORS     <- c("#2166ac", "#b2182b")

CENTROFRONTAL  <- c("F3", "FCz", "FC1", "Cz")
LEFT_ANTERIOR  <- c("F7", "F3", "FC5")
CENTROPARIETAL <- c("Cz", "CP1", "CP2", "Pz")

CONTRASTS <- list(
  noun_cue = list(
    event = "noun_F0", type = "cue",
    levels = list(paroxytone = 1:3, oxytone = 4:6),
    electrodes = CENTROFRONTAL, peak = "F3"),
  verb_stress = list(
    event = "verb_F0", type = "cue",
    levels = list(paroxytone = c(1, 4), oxytone = c(2, 6)),
    electrodes = CENTROFRONTAL, peak = "F3"),
  noun_violation = list(
    event = "verb_F0", type = "violation",
    levels = list(match = c(1, 6), mismatch = c(2, 4)),
    electrodes = c(LEFT_ANTERIOR, CENTROPARIETAL), peak = c("F7", "Pz")),
  verb_violation = list(
    event = "suffix", type = "violation",
    levels = list(match = c(1, 6), mismatch = c(3, 5)),
    electrodes = c(LEFT_ANTERIOR, CENTROPARIETAL), peak = c("F3", "Pz"))
)

WINDOWS <- bind_rows(
  tibble(type = "cue", channel = CENTROFRONTAL, start = 136, end = 280),
  tibble(type = "violation", channel = LEFT_ANTERIOR, start = 225, end = 300),
  tibble(type = "violation", channel = CENTROPARIETAL, start = 300, end = 500),
  tibble(type = "violation", channel = CENTROPARIETAL, start = 400, end = 700)
)

dir_create(FIGURE_DIR)
theme_set(theme_minimal(base_size = 11))

# data
waves <- read_csv(WAVE_FILE, show_col_types = FALSE) |>
  mutate(group   = factor(group, levels = c("L1", "L2")),
         session = factor(session, levels = c("exp_pre", "exp_post"),
                          labels = c("pre", "post")))

# helper functions
grand_average <- function(contrast) {
  levels <- enframe(contrast$levels, "level", "cond") |> unnest(cond)

  waves |>
    filter(event == contrast$event, channel %in% contrast$electrodes) |>
    inner_join(levels, by = "cond") |>
    group_by(participant, group, session, channel, time_ms, level) |>
    summarise(amplitude = mean(amplitude), .groups = "drop") |>
    group_by(group, session, channel, time_ms, level) |>
    summarise(mean = mean(amplitude), error = sd(amplitude) / sqrt(n()),
              .groups = "drop") |>
    mutate(level   = factor(level, levels = names(contrast$levels)),
           channel = factor(channel, levels = contrast$electrodes))
}

plot_waves <- function(average, type) {
  windows <- WINDOWS |>
    filter(type == !!type, channel %in% average$channel) |>
    mutate(channel = factor(channel, levels = levels(average$channel)))

  ggplot(average, aes(time_ms, mean, colour = level, fill = level)) +
    geom_rect(data = windows, inherit.aes = FALSE,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              fill = "grey85", alpha = 0.5) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_ribbon(aes(ymin = mean - error, ymax = mean + error),
                alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.6) +
    scale_colour_manual(values = COLORS) +
    scale_fill_manual(values = COLORS) +
    labs(x = "time (ms)", y = "amplitude (µV)", colour = NULL, fill = NULL) +
    theme(legend.position = "top", panel.grid.minor = element_blank())
}

# grand averages
averages <- map(CONTRASTS, grand_average)

# plots of every roi electrode
ALL_FIGURES <- imap(CONTRASTS, \(contrast, name) {
  plot_waves(averages[[name]], contrast$type) +
    facet_grid(group + session ~ channel)
})

for (figure_name in names(ALL_FIGURES)) {
  print(ALL_FIGURES[[figure_name]])
  ggsave(path(FIGURE_DIR, paste0("erp_electrodes_", figure_name, ".png")),
         ALL_FIGURES[[figure_name]], dpi = 300, height = 7,
         width = 2.2 * length(CONTRASTS[[figure_name]]$electrodes) + 1)
}

# plots of the electrode with the largest effect
PEAK_FIGURES <- imap(CONTRASTS, \(contrast, name) {
  map(set_names(contrast$peak, paste(name, contrast$peak, sep = "_")), \(electrode) {
    averages[[name]] |>
      filter(channel == electrode) |>
      plot_waves(contrast$type) +
      facet_grid(group ~ session) +
      labs(title = electrode)
  })
}) |>
  list_flatten(name_spec = "{inner}")

for (figure_name in names(PEAK_FIGURES)) {
  print(PEAK_FIGURES[[figure_name]])
  ggsave(path(FIGURE_DIR, paste0("erp_wave_", figure_name, ".png")),
         PEAK_FIGURES[[figure_name]], width = 6.5, height = 5, dpi = 300)
}

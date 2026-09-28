# grand average waveform figures

# libraries
library(tidyverse)
library(fs)

# configuration
WAVE_FILE  <- "~/Kaylee_SP26/sp26_waveforms_event.csv"
OUTPUT_DIR <- "~/Desktop/SP26_data/erp_output"
FIG_DIR    <- path(OUTPUT_DIR, "waveforms")
dir_create(FIG_DIR)

# condition 1-6 maps to letters a-f
COND_LETTER <- c("a","b","c","d","e","f")

# measurement windows, for shading
WINDOWS <- tibble::tribble(
  ~event,    ~roi,             ~component,  ~w0, ~w1,
  "noun_F0", "centrofrontal",  "PrAN",      136, 280,
  "verb_F0", "centrofrontal",  "PrAN",      136, 280,
  "verb_F0", "left_anterior",  "LAN",       225, 300,
  "verb_F0", "centroparietal", "N400",      300, 500,
  "verb_F0", "left_anterior",  "P600",      400, 700,
  "suffix",  "left_anterior",  "LAN",       225, 300,
  "suffix",  "centroparietal", "N400",      300, 500,
  "suffix",  "left_anterior",  "P600",      400, 700
)

# data
w <- read_csv(WAVE_FILE, show_col_types = FALSE) |>
  mutate(letter   = COND_LETTER[cond],
         stress   = if_else(letter %in% c("a","c","d"), "paroxytone", "oxytone"),
         session  = factor(session, levels = c("exp_pre","exp_post"),
                           labels = c("pre","post")),
         group    = factor(group, levels = c("L1","L2")),
         stress   = factor(stress, levels = c("paroxytone","oxytone")))

message(nrow(w), " rows, ", n_distinct(w$participant), " participants")

# helper functions
# verb validity: a,f valid / c,e invalid.  noun validity: a,f valid / b,d invalid.
grand_avg <- function(d, ...) {
  d |>
    group_by(participant, group, session, time_ms, ...) |>
    summarise(amp = mean(amplitude), .groups = "drop") |>
    group_by(group, session, time_ms, ...) |>
    summarise(mean_amp = mean(amp),
              se = sd(amp) / sqrt(n()),
              n = n(), .groups = "drop")
}

shade <- function(ev, roi_name) {
  WINDOWS |> filter(event == ev, roi == roi_name)
}

wave_plot <- function(ga, ev, roi_name, title) {
  bands <- shade(ev, roi_name)
  ggplot(ga, aes(time_ms, mean_amp, colour = validity, fill = validity)) +
    geom_rect(data = bands, inherit.aes = FALSE,
              aes(xmin = w0, xmax = w1, ymin = -Inf, ymax = Inf),
              fill = "grey85", alpha = 0.5) +
    geom_text(data = bands, inherit.aes = FALSE,
              aes(x = (w0 + w1) / 2, y = Inf, label = component),
              vjust = 1.4, size = 3, colour = "grey35") +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_ribbon(aes(ymin = mean_amp - se, ymax = mean_amp + se),
                alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.7) +
    facet_grid(group + session ~ stress) +
    scale_colour_manual(values = c(valid = "#2166ac", invalid = "#b2182b")) +
    scale_fill_manual(values = c(valid = "#2166ac", invalid = "#b2182b")) +
    labs(title = title, x = "time (ms)", y = "amplitude (uV)",
         colour = NULL, fill = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top")
}

figs <- list()

# verb violation at the suffix
verb_suffix <- w |>
  filter(event == "suffix", letter %in% c("a","c","e","f")) |>
  mutate(validity = factor(if_else(letter %in% c("a","f"), "valid", "invalid"),
                           levels = c("valid","invalid")))

for (r in c("centroparietal", "left_anterior")) {
  ga <- grand_avg(filter(verb_suffix, roi == r), stress, validity)
  p <- wave_plot(ga, "suffix", r,
                 paste0("verb violation at the suffix, ", r))
  figs[[paste0("verb_suffix_", r)]] <- p
  ggsave(path(FIG_DIR, paste0("verb_suffix_", r, ".png")), p,
         width = 10, height = 9, dpi = 150)
}

# noun violation at verb F0
noun_verbF0 <- w |>
  filter(event == "verb_F0", letter %in% c("a","b","d","f")) |>
  mutate(validity = factor(if_else(letter %in% c("a","f"), "valid", "invalid"),
                           levels = c("valid","invalid")))

for (r in c("centroparietal", "left_anterior")) {
  ga <- grand_avg(filter(noun_verbF0, roi == r), stress, validity)
  p <- wave_plot(ga, "verb_F0", r,
                 paste0("noun violation at verb F0, ", r))
  figs[[paste0("noun_verbF0_", r)]] <- p
  ggsave(path(FIG_DIR, paste0("noun_verbF0_", r, ".png")), p,
         width = 10, height = 9, dpi = 150)
}

# verb stress PrAN at verb F0
pran <- w |>
  filter(event == "verb_F0", roi == "centrofrontal",
         letter %in% c("a","b","d","f"))

ga_pran <- pran |>
  group_by(participant, group, session, time_ms, stress) |>
  summarise(amp = mean(amplitude), .groups = "drop") |>
  group_by(group, session, time_ms, stress) |>
  summarise(mean_amp = mean(amp), se = sd(amp) / sqrt(n()), .groups = "drop")

bands <- shade("verb_F0", "centrofrontal")
p_pran <- ggplot(ga_pran, aes(time_ms, mean_amp, colour = stress, fill = stress)) +
  geom_rect(data = bands, inherit.aes = FALSE,
            aes(xmin = w0, xmax = w1, ymin = -Inf, ymax = Inf),
            fill = "grey85", alpha = 0.5) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_ribbon(aes(ymin = mean_amp - se, ymax = mean_amp + se),
              alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.7) +
  facet_grid(group ~ session) +
  scale_colour_manual(values = c(paroxytone = "#1b7837", oxytone = "#762a83")) +
  scale_fill_manual(values = c(paroxytone = "#1b7837", oxytone = "#762a83")) +
  labs(title = "verb stress PrAN at verb F0, centrofrontal",
       x = "time (ms)", y = "amplitude (uV)", colour = NULL, fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")
figs[["pran_verb_stress"]] <- p_pran
ggsave(path(FIG_DIR, "pran_verb_stress.png"), p_pran,
       width = 9, height = 7, dpi = 150)

# noun cue PrAN at noun F0
cue <- w |>
  filter(event == "noun_F0", roi == "centrofrontal") |>
  mutate(cue = factor(if_else(letter %in% c("a","b","c"), "a", "c")))

ga_cue <- cue |>
  group_by(participant, group, session, time_ms, cue) |>
  summarise(amp = mean(amplitude), .groups = "drop") |>
  group_by(group, session, time_ms, cue) |>
  summarise(mean_amp = mean(amp), se = sd(amp) / sqrt(n()), .groups = "drop")

bands <- shade("noun_F0", "centrofrontal")
p_cue <- ggplot(ga_cue, aes(time_ms, mean_amp, colour = cue, fill = cue)) +
  geom_rect(data = bands, inherit.aes = FALSE,
            aes(xmin = w0, xmax = w1, ymin = -Inf, ymax = Inf),
            fill = "grey85", alpha = 0.5) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_ribbon(aes(ymin = mean_amp - se, ymax = mean_amp + se),
              alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.7) +
  facet_grid(group ~ session) +
  labs(title = "noun cue PrAN at noun F0, centrofrontal",
       x = "time (ms)", y = "amplitude (uV)", colour = NULL, fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")
figs[["pran_noun_cue"]] <- p_cue
ggsave(path(FIG_DIR, "pran_noun_cue.png"), p_cue,
       width = 9, height = 7, dpi = 150)

# print
for (nm in names(figs)) print(figs[[nm]])

# results
message("saved to ", FIG_DIR)
message(length(figs), " figures. call figs$<name> to redraw one:")
message("  ", paste(names(figs), collapse = ", "))

# erp topographies

# libraries
library(tidyverse)
library(fs)

# configuration
MEANS_FILE <- "~/Kaylee_SP26/sp26_topography_means.csv"
FIGURE_DIR <- "figures"
GRID_STEP  <- 0.01

TITLES <- c(
  PrAN_noun = "noun cue: oxytone minus paroxytone, 136 to 280 ms",
  PrAN_verb = "verb stress: oxytone minus paroxytone, 136 to 280 ms",
  LAN_noun  = "noun mismatch minus match, 225 to 300 ms",
  N400_noun = "noun mismatch minus match, 300 to 500 ms",
  P600_noun = "noun mismatch minus match, 400 to 700 ms",
  LAN_verb  = "verb mismatch minus match, 225 to 300 ms",
  N400_verb = "verb mismatch minus match, 300 to 500 ms",
  P600_verb = "verb mismatch minus match, 400 to 700 ms")

CENTROFRONTAL  <- c("F3", "FCz", "FC1", "Cz")
LEFT_ANTERIOR  <- c("F7", "F3", "FC5")
CENTROPARIETAL <- c("Cz", "CP1", "CP2", "Pz")

ROIS <- list(PrAN_noun = CENTROFRONTAL,  PrAN_verb = CENTROFRONTAL,
             LAN_noun  = LEFT_ANTERIOR,  N400_noun = CENTROPARIETAL,
             P600_noun = CENTROPARIETAL, LAN_verb  = LEFT_ANTERIOR,
             N400_verb = CENTROPARIETAL, P600_verb = CENTROPARIETAL)

dir_create(FIGURE_DIR)

# data
means <- read_csv(MEANS_FILE, show_col_types = FALSE) |>
  mutate(group   = factor(group, levels = c("L1", "L2")),
         session = factor(session, levels = c("pre", "post")),
         x       = radius * sin(theta * pi / 180),
         y       = radius * cos(theta * pi / 180))

HEAD_RADIUS <- max(means$radius) * 1.05

grid <- expand_grid(x = seq(-HEAD_RADIUS, HEAD_RADIUS, GRID_STEP),
                    y = seq(-HEAD_RADIUS, HEAD_RADIUS, GRID_STEP)) |>
  filter(sqrt(x^2 + y^2) <= HEAD_RADIUS)

outline <- tibble(angle = seq(0, 2 * pi, length.out = 200),
                  x = HEAD_RADIUS * sin(angle),
                  y = HEAD_RADIUS * cos(angle))

nose <- tibble(x = c(-0.08, 0, 0.08) * HEAD_RADIUS,
               y = c(0.997, 1.12, 0.997) * HEAD_RADIUS)

# helper functions
distance <- function(a, b) {
  sqrt(outer(a$x, b$x, "-")^2 + outer(a$y, b$y, "-")^2)
}

green <- function(d) {
  ifelse(d == 0, 0, d^2 * (log(d) - 1))
}

interpolate <- function(channels, ...) {
  weights <- solve(green(distance(channels, channels)), channels$amplitude)
  mutate(grid, amplitude = as.vector(green(distance(grid, channels)) %*% weights))
}

plot_topography <- function(name) {
  channels <- filter(means, component == name) |>
    mutate(roi = channel %in% ROIS[[name]])
  surface  <- channels |>
    group_by(group, session) |>
    group_modify(interpolate) |>
    ungroup()
  limit <- max(abs(surface$amplitude))

  ggplot(surface, aes(x, y)) +
    geom_raster(aes(fill = amplitude)) +
    geom_contour(aes(z = amplitude), colour = "grey40", linewidth = 0.15, bins = 7) +
    geom_path(data = outline, linewidth = 0.5) +
    geom_path(data = nose, linewidth = 0.5) +
    geom_point(data = channels, aes(size = roi)) +
    scale_size_manual(values = c(`FALSE` = 0.5, `TRUE` = 2.2), guide = "none") +
    facet_grid(group ~ session, switch = "y") +
    scale_fill_distiller(palette = "RdBu", limits = c(-limit, limit)) +
    coord_equal(clip = "off") +
    labs(title = TITLES[[name]], fill = "µV",
         caption = paste("ROI (large dots):", paste(ROIS[[name]], collapse = ", "))) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(size = 11, margin = margin(b = 10)),
          plot.caption = element_text(size = 10, hjust = 0.5),
          strip.text = element_text(size = 11, margin = margin(4, 4, 8, 4)),
          panel.spacing = unit(1.2, "lines"),
          plot.margin = margin(10, 10, 10, 10))
}

# plots
FIGURES <- map(set_names(names(TITLES)), plot_topography)

for (figure_name in names(FIGURES)) {
  print(FIGURES[[figure_name]])
  ggsave(path(FIGURE_DIR, paste0("erp_topography_", figure_name, ".png")),
         FIGURES[[figure_name]], width = 5.5, height = 5, dpi = 300, bg = "white")
}

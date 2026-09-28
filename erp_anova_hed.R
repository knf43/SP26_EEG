# ============================================================
# Hed-style repeated-measures ANOVAs on the SP26 ERP data
#
# Hed et al. (2019, p. 109): "the same procedure as in Roll et al. (2015)
# was used. Thus, repeated measure ANOVAs were conducted with factors
# session (pre and post), tone (accent 1 and accent 2) and validity (valid,
# invalid), using the average value for each participant and condition for
# the ANOVA."
#
# So: average within participant x condition FIRST, then analyse the means.
# This averages away the trial-level noise that makes the single-trial LMMs
# look so weak (marginal R2 ~ .002), and is the analysis the source paper
# would recognise.
#
# Group (L1 vs L2) is added as a between-subjects factor; Hed had only L2
# learners.
#
# EXPLORATORY. The confirmatory analysis is erp_statistical_analysis.R.
# This is here to see the shape of the data, not to replace it. Note that a
# group x session x tone x validity ANOVA has 15 main effects and
# interactions, so its familywise error rate is around .54 (Luck & Gaspelin
# 2017, p. 10). Treat anything here as descriptive.
# ============================================================

library(tidyverse)
library(fs)

OUTPUT_DIR <- "~/Desktop/SP26_data/erp_output"
ANOVA_DIR  <- path(OUTPUT_DIR, "anova")
dir_create(ANOVA_DIR)

has_afex <- requireNamespace("afex", quietly = TRUE)
if (!has_afex) {
  message("afex is not installed. Falling back to base aov(), which gives no ",
          "Greenhouse-Geisser correction.\n  install.packages('afex') for the ",
          "standard treatment.")
}

# reuse `erp` if erp_script.R has already been sourced this session
if (!exists("erp")) source("erp_script.R")
dat <- erp$erp

# ---- the two 2x2 designs -------------------------------------
# verb validity:  a = parox valid, c = parox invalid, f = oxy valid, e = oxy invalid
# noun validity:  a = parox valid, d = parox invalid, f = oxy valid, b = oxy invalid
#
# A plain list, not a tibble: storing character vectors in a tribble makes
# list columns, which dplyr then refuses to print or bind.
design_spec <- function(which_design) {
  switch(which_design,
    verb = list(conds = c("a","c","e","f"), valid = c("a","f"), invalid = c("c","e")),
    noun = list(conds = c("a","b","d","f"), valid = c("a","f"), invalid = c("d","b")),
    stop("design_spec: unknown design '", which_design, "'")
  )
}

COMPONENT_DESIGN <- tibble::tribble(
  ~component,  ~design,
  "LAN_noun",  "noun",
  "N400_noun", "noun",
  "P600_noun", "noun",
  "LAN_verb",  "verb",
  "N400_verb", "verb",
  "P600_verb", "verb",
  "PrAN_verb", "noun",   # stress is the effect of interest; validity carried along
  "PrAN_noun", "cue"     # handled separately: splice a vs c, not a 2x2
)

# ---- participant x cell means, exactly what Hed analysed ------
cell_means <- function(comp, which_design) {
  spec <- design_spec(which_design)
  d <- dat |>
    filter(component == comp, erp_phase == "experimental",
           cond %in% spec$conds)

  d |>
    mutate(
      session  = factor(as.character(phase), levels = c("exp_pre", "exp_post"),
                        labels = c("pre", "post")),
      tone     = factor(verb_stress, levels = c("paroxytone", "oxytone")),
      validity = factor(if_else(cond %in% spec$valid, "valid", "invalid"),
                        levels = c("valid", "invalid")),
      group    = factor(group, levels = c("L1", "L2")),
      participant = factor(participant)
    ) |>
    group_by(participant, group, session, tone, validity) |>
    summarise(amp = mean(amplitude, na.rm = TRUE),
              n_trials = n(), .groups = "drop")
}

run_anova <- function(cm, label) {
  if (has_afex) {
    fit <- afex::aov_ez(id = "participant", dv = "amp", data = cm,
                        between = "group",
                        within  = c("session", "tone", "validity"),
                        type = 3, include_aov = FALSE)
    as.data.frame(fit$anova_table) |>
      rownames_to_column("effect") |>
      transmute(model = label, effect,
                df1 = `num Df`, df2 = `den Df`,
                F = round(F, 2), p = signif(`Pr(>F)`, 3),
                ges = round(ges, 3))
  } else {
    fit <- aov(amp ~ group * session * tone * validity +
                 Error(participant / (session * tone * validity)), data = cm)
    s <- summary(fit)
    map_dfr(names(s), function(stratum) {
      tb <- as.data.frame(s[[stratum]][[1]]) |> rownames_to_column("effect")
      tb |> filter(!str_detect(effect, "Residuals")) |>
        transmute(model = label, effect = str_trim(effect),
                  df1 = Df, df2 = NA_real_,
                  F = round(`F value`, 2), p = signif(`Pr(>F)`, 3),
                  ges = NA_real_)
    })
  }
}

# ---- run --------------------------------------------------------
res_l <- list(); means_l <- list()

for (i in seq_len(nrow(COMPONENT_DESIGN))) {
  comp <- COMPONENT_DESIGN$component[i]
  des  <- COMPONENT_DESIGN$design[i]
  if (des == "cue") next          # PrAN_noun handled below

  message("\n=== ", comp, "  (", des, " validity design) ===")
  cm <- cell_means(comp, des)
  means_l[[comp]] <- cm |> mutate(model = comp)
  res_l[[comp]]   <- run_anova(cm, comp)
  print(as_tibble(res_l[[comp]]), n = 20)
}

# ---- PrAN_noun: the noun cue, splice a vs c ---------------------
message("\n=== PrAN_noun  (noun cue: splice a vs c) ===")
cm_cue <- dat |>
  filter(component == "PrAN_noun", erp_phase == "experimental",
         !is.na(splice)) |>
  mutate(session = factor(as.character(phase),
                          levels = c("exp_pre", "exp_post"),
                          labels = c("pre", "post")),
         cue     = factor(splice, levels = c("a", "c")),
         group   = factor(group, levels = c("L1", "L2")),
         participant = factor(participant)) |>
  group_by(participant, group, session, cue) |>
  summarise(amp = mean(amplitude, na.rm = TRUE),
            n_trials = n(), .groups = "drop")

res_cue <- if (has_afex) {
  fit <- afex::aov_ez(id = "participant", dv = "amp", data = cm_cue,
                      between = "group", within = c("session", "cue"),
                      type = 3, include_aov = FALSE)
  as.data.frame(fit$anova_table) |>
    rownames_to_column("effect") |>
    transmute(model = "PrAN_noun", effect,
              df1 = `num Df`, df2 = `den Df`,
              F = round(F, 2), p = signif(`Pr(>F)`, 3),
              ges = round(ges, 3))
} else {
  fit <- aov(amp ~ group * session * cue +
               Error(participant / (session * cue)), data = cm_cue)
  s <- summary(fit)
  map_dfr(names(s), function(stratum) {
    tb <- as.data.frame(s[[stratum]][[1]]) |> rownames_to_column("effect")
    tb |> filter(!str_detect(effect, "Residuals")) |>
      transmute(model = "PrAN_noun", effect = str_trim(effect),
                df1 = Df, df2 = NA_real_,
                F = round(`F value`, 2), p = signif(`Pr(>F)`, 3),
                ges = NA_real_)
  })
}
print(as_tibble(res_cue), n = 20)

# ---- training growth, simple correlation -----------------------
# Hed correlated PrAN with behaviour rather than modelling trials; here the
# simplest look is whether mean training PrAN differs by group.
cm_train <- dat |>
  filter(component == "PrAN_verb", erp_phase == "training") |>
  mutate(group = factor(group, levels = c("L1", "L2")),
         half  = if_else(trial <= median(trial), "first", "second")) |>
  group_by(participant, group, half) |>
  summarise(amp = mean(amplitude, na.rm = TRUE), .groups = "drop")

message("\n=== training PrAN, first vs second half ===")
print(cm_train |> group_by(group, half) |>
        summarise(mean_amp = round(mean(amp), 2),
                  se = round(sd(amp) / sqrt(n()), 2), .groups = "drop"))

# ---- collect ----------------------------------------------------
res_l[["PrAN_noun"]] <- res_cue
anova_all <- bind_rows(res_l)
means_all <- bind_rows(means_l)

write_csv(anova_all, path(ANOVA_DIR, "hed_style_anova.csv"))
write_csv(means_all, path(ANOVA_DIR, "cell_means.csv"))
write_csv(cm_cue,    path(ANOVA_DIR, "cell_means_cue.csv"))

message("\n--- everything at p < .05 ---")
print(as_tibble(anova_all) |> filter(p < .05), n = 50)

message("\nSaved to ", ANOVA_DIR)

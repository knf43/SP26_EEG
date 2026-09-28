# anovas on participant x cell means

# libraries
library(tidyverse)
library(fs)

# configuration
OUTPUT_DIR <- "~/Desktop/SP26_data/erp_output"
ANOVA_DIR  <- path(OUTPUT_DIR, "anova")
dir_create(ANOVA_DIR)
has_afex <- requireNamespace("afex", quietly = TRUE)
if (!has_afex) message("afex not installed; falling back to aov() with no GG correction")

# data
if (!exists("erp")) source("erp_script.R")
dat <- erp$erp

# component definitions
# noun validity: a,f valid / b,d invalid.  verb validity: a,f valid / c,e invalid.
# stress: a,c,d paroxytone / b,e,f oxytone.
spec_for <- function(comp) {
  switch(comp,
    PrAN_noun = list(conds = c("a","b","c","d","e","f"), within = c("session","cue")),
    PrAN_verb = list(conds = c("a","b","d","f"),         within = c("session","stress")),
    LAN_noun  = , N400_noun = ,
    P600_noun = list(conds = c("a","b","d","f"),         within = c("session","stress","validity")),
    LAN_verb  = , N400_verb = ,
    P600_verb = list(conds = c("a","c","e","f"),         within = c("session","stress","validity"))
  )
}

COMPONENTS <- c("PrAN_noun","PrAN_verb",
                "LAN_noun","N400_noun","P600_noun",
                "LAN_verb","N400_verb","P600_verb")

# helper functions
code_factors <- function(d) {
  d |> mutate(session  = factor(as.character(phase),
                                levels = c("exp_pre","exp_post"),
                                labels = c("pre","post")),
              stress   = factor(verb_stress, levels = c("paroxytone","oxytone")),
              cue      = factor(splice, levels = c("a","c")),
              validity = factor(if_else(cond %in% c("a","f"), "valid", "invalid"),
                                levels = c("valid","invalid")),
              group    = factor(group, levels = c("L1","L2")),
              participant = factor(participant))
}

cell_means <- function(comp) {
  sp <- spec_for(comp)
  dat |>
    filter(component == comp, erp_phase == "experimental", cond %in% sp$conds) |>
    code_factors() |>
    group_by(across(all_of(c("participant","group", sp$within)))) |>
    summarise(amp = mean(amplitude, na.rm = TRUE), n = n(), .groups = "drop")
}

run_anova <- function(cm, within, label, between = "group") {
  if (has_afex) {
    fit <- afex::aov_ez(id = "participant", dv = "amp", data = cm,
                        between = between, within = within,
                        type = 3, include_aov = FALSE)
    as.data.frame(fit$anova_table) |>
      rownames_to_column("effect") |>
      transmute(model = label, effect,
                df1 = `num Df`, df2 = `den Df`,
                F = round(F, 2), p = signif(`Pr(>F)`, 3), ges = round(ges, 3))
  } else {
    w <- paste(within, collapse = " * ")
    f <- as.formula(paste("amp ~", between, "*", w,
                          "+ Error(participant/(", w, "))"))
    s <- summary(aov(f, data = cm))
    map_dfr(names(s), function(st) {
      as.data.frame(s[[st]][[1]]) |> rownames_to_column("effect") |>
        filter(!str_detect(effect, "Residuals")) |>
        transmute(model = label, effect = str_trim(effect),
                  df1 = Df, df2 = NA_real_,
                  F = round(`F value`, 2), p = signif(`Pr(>F)`, 3), ges = NA_real_)
    })
  }
}

show_means <- function(cm, within) {
  cm |> group_by(across(all_of(c("group", within)))) |>
    summarise(mean_uV = round(mean(amp), 2),
              se = round(sd(amp)/sqrt(n()), 2), .groups = "drop") |>
    print(n = 40)
}

# experimental
res_l <- list(); means_l <- list()

for (comp in COMPONENTS) {
  sp <- spec_for(comp)
  cm <- cell_means(comp)
  means_l[[comp]] <- cm |> mutate(model = comp)
  res_l[[comp]]   <- run_anova(cm, sp$within, comp)

  cat("\n=== ", comp, "  group x ", paste(sp$within, collapse = " x "), " ===\n", sep = "")
  print(as_tibble(res_l[[comp]]), n = 20)
  show_means(cm, sp$within)
}

# training
# no validity manipulation in training: 60 stressed->present, 60 unstressed->past,
# all valid. so the factors are stress and half of training, plus group.
cm_train <- dat |>
  filter(component == "PrAN_verb", erp_phase == "training") |>
  code_factors() |>
  group_by(participant) |>
  mutate(half = factor(if_else(trial <= median(trial), "first", "second"),
                       levels = c("first","second"))) |>
  ungroup() |>
  group_by(participant, group, half, stress) |>
  summarise(amp = mean(amplitude, na.rm = TRUE), n = n(), .groups = "drop")

cat("\n=== PrAN_train  group x half x stress ===\n")
res_train <- run_anova(cm_train, c("half","stress"), "PrAN_train")
print(as_tibble(res_train), n = 20)
show_means(cm_train, c("half","stress"))

# results
res_l[["PrAN_train"]] <- res_train
anova_all <- bind_rows(res_l)
means_all <- bind_rows(means_l)

cat("\n--- p < .05 ---\n")
print(as_tibble(anova_all) |> filter(p < .05), n = 60)

cat("\n--- largest effect sizes ---\n")
print(as_tibble(anova_all) |> filter(!is.na(ges)) |> arrange(desc(ges)) |> head(12))

# save
write_csv(anova_all, path(ANOVA_DIR, "anova_results.csv"))
write_csv(means_all, path(ANOVA_DIR, "cell_means.csv"))
write_csv(cm_train,  path(ANOVA_DIR, "cell_means_training.csv"))
message("\nsaved to ", ANOVA_DIR)

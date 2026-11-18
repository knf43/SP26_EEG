# erp models, built up in levels

# libraries
library(tidyverse)
library(lme4)
library(lmerTest)
library(DHARMa)
library(MuMIn)
library(broom.mixed)
library(ggeffects)
library(openxlsx)
library(fs)

# configuration
# fixed effects are always group * phase * validity * stress.
# LEVEL controls the random structure, which is what causes convergence trouble.
# LEVEL 1  (1 | participant) + (1 | item)
# LEVEL 2  (1 + validity || participant) + (1 | item)
# LEVEL 3  (1 + validity || participant) + (1 + validity || item)
# LEVEL 4  (1 + phase * validity || participant) + (1 + phase * validity || item)
LEVEL      <- 1
OUTPUT_DIR <- "~/Desktop/SP26_data/erp_output"
LMM_DIR    <- path(OUTPUT_DIR, paste0("lmm_L", LEVEL))
dir_create(LMM_DIR)
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

# data
source("erp_script.R")
dat <- erp$erp

# component definitions
# noun validity: a,f valid / b,d invalid.  verb validity: a,f valid / c,e invalid.
# stress: a,c,d paroxytone / b,e,f oxytone.
spec_for <- function(comp) {
  switch(comp,
    PrAN_noun = list(conds = c("a","b","c","d","e","f"), preds = "cue"),
    PrAN_verb = list(conds = c("a","b","d","f"),         preds = "stress"),
    LAN_noun  = , N400_noun = ,
    P600_noun = list(conds = c("a","b","d","f"),         preds = c("validity","stress")),
    LAN_verb  = , N400_verb = ,
    P600_verb = list(conds = c("a","c","e","f"),         preds = c("validity","stress"))
  )
}

COMPONENTS <- c("PrAN_noun","PrAN_verb",
                "LAN_noun","N400_noun","P600_noun",
                "LAN_verb","N400_verb","P600_verb")

# helper functions
prep <- function(comp) {
  sp <- spec_for(comp)
  d <- dat |>
    filter(component == comp, erp_phase == "experimental", cond %in% sp$conds) |>
    mutate(group    = factor(group, levels = c("L1","L2")),
           phase    = factor(as.character(phase), levels = c("exp_pre","exp_post")),
           stress   = factor(verb_stress, levels = c("paroxytone","oxytone")),
           cue      = factor(splice, levels = c("a","c")),
           validity = factor(if_else(cond %in% c("a","f"), "valid", "invalid"),
                             levels = c("valid","invalid")),
           item = factor(item), participant = factor(participant)) |>
    drop_na(any_of(c("phase","item","participant", sp$preds)))
  contrasts(d$group) <- contr.sum(2)
  contrasts(d$phase) <- contr.sum(2)
  for (p in sp$preds) contrasts(d[[p]]) <- contr.sum(nlevels(d[[p]]))
  d
}

factors_at <- function(preds, level) c("group", "phase", preds)

formula_at <- function(preds, level) {
  fx  <- paste(factors_at(preds, level), collapse = " * ")
  key <- preds[1]
  re  <- switch(as.character(level),
    "1" = "(1 | participant) + (1 | item)",
    "2" = sprintf("(1 + %s || participant) + (1 | item)", key),
    "3" = sprintf("(1 + %s || participant) + (1 + %s || item)", key, key),
    "4" = sprintf("(1 + phase * %s || participant) + (1 + phase * %s || item)", key, key))
  as.formula(paste("amplitude ~", fx, "+", re))
}

drop_test <- function(bigger, term, label) {
  red <- update(bigger, as.formula(paste(". ~ . -", term)))
  cmp <- as.data.frame(anova(red, bigger, refit = TRUE))
  nm  <- names(cmp)
  c_df <- if ("Chi Df" %in% nm) "Chi Df" else "Df"
  data.frame(model = label, effect = term,
             chisq = round(cmp[[nm[grepl("Chisq", nm)][1]]][2], 2),
             df    = cmp[[c_df]][2],
             p     = signif(cmp[[nm[grepl("^Pr", nm)][1]]][2], 3))
}

nmc <- function(full, label, factors) {
  orders <- lapply(rev(seq_along(factors)), function(k)
    combn(factors, k, FUN = \(x) paste(x, collapse = ":")))
  out <- list(); cur <- full
  for (i in seq_along(orders)) {
    tms <- orders[[i]]
    for (tm in tms) out[[length(out) + 1]] <- drop_test(cur, tm, label)
    if (i < length(orders))
      cur <- update(cur, as.formula(paste(". ~ . -", paste(tms, collapse = " - "))))
  }
  bind_rows(out)
}

diagnose <- function(m, label) {
  sim <- simulateResiduals(m, n = 1000)
  png(path(LMM_DIR, paste0("dharma_", label, ".png")), width = 1200, height = 600, res = 120)
  plot(sim); dev.off()
  disp <- testDispersion(sim, plot = FALSE)
  data.frame(model = label,
             KS_pvalue    = signif(testUniformity(sim, plot = FALSE)$p.value, 3),
             dispersion   = signif(disp$statistic, 3),
             dispersion_p = signif(disp$p.value, 3))
}

get_fixed <- function(m, label) {
  tidy(m, effects = "fixed", conf.int = TRUE) |>
    mutate(model = label,
           across(c(estimate, std.error, conf.low, conf.high), \(x) round(x, 3)),
           p.value = signif(p.value, 3)) |>
    select(model, term, estimate, std.error, conf.low, conf.high, p.value)
}

# models
models <- list(); fixed_l <- list(); nmc_l <- list()
diag_l <- list(); r2_l <- list(); pred_l <- list(); struct_l <- list()

message("\nLEVEL ", LEVEL)

for (comp in COMPONENTS) {
  sp <- spec_for(comp)
  fx <- factors_at(sp$preds, LEVEL)
  message("\n=== ", comp, "  ", paste(fx, collapse = " * "), " ===")
  t0 <- Sys.time()

  d <- prep(comp)
  f <- formula_at(sp$preds, LEVEL)
  m <- suppressWarnings(lmer(f, data = d, control = ctrl))
  if (isSingular(m, tol = 1e-5)) message("  singular")

  models[[comp]]  <- m
  fixed_l[[comp]] <- get_fixed(m, comp)
  nmc_l[[comp]]   <- nmc(m, comp, fx)
  diag_l[[comp]]  <- diagnose(m, comp)

  r2 <- r.squaredGLMM(m)
  r2_l[[comp]] <- data.frame(model = comp,
                             R2_marginal    = round(r2[1,"R2m"], 3),
                             R2_conditional = round(r2[1,"R2c"], 3),
                             row.names = NULL)
  struct_l[[comp]] <- data.frame(model = comp, level = LEVEL,
                                 formula = deparse1(f),
                                 n_trials = nrow(d),
                                 singular = isSingular(m, tol = 1e-5),
                                 row.names = NULL)

  pred_l[[comp]] <- as.data.frame(ggpredict(m, terms = rev(sp$preds))) |>
    mutate(model = comp)

  message("  ", nrow(d), " trials, ",
          round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s")
}

# results
fixed_all           <- bind_rows(fixed_l)
nmc_all             <- bind_rows(nmc_l)
r2_summary          <- bind_rows(r2_l)
diagnostics_summary <- bind_rows(diag_l)
model_structure     <- bind_rows(struct_l)
pred_all            <- bind_rows(pred_l)

cat("\n--- all effects ---\n")
print(as_tibble(nmc_all), n = 100)
cat("\n--- p < .05 ---\n")
print(as_tibble(nmc_all) |> filter(p < .05), n = 60)

# plots
theme_set(theme_minimal(base_size = 12))

p_pred <- ggplot(pred_all, aes(x, predicted,
                               colour = if ("group" %in% names(pred_all)) group else NULL)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                  position = position_dodge(0.25), size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  facet_wrap(~ model, scales = "free_y") +
  labs(x = NULL, y = "predicted amplitude (uV)", colour = NULL)
print(p_pred)
ggsave(path(LMM_DIR, "pred_amplitude.png"), p_pred, width = 11, height = 8, dpi = 150)

# save
write.xlsx(list(fixed_effects = fixed_all, nested_comparisons = nmc_all,
                model_fit_R2 = r2_summary, diagnostics = diagnostics_summary,
                model_structure = model_structure, predicted = pred_all),
           path(LMM_DIR, "erp_results.xlsx"), overwrite = TRUE)

saveRDS(list(level = LEVEL, models = models, fixed_all = fixed_all,
             nmc_all = nmc_all, r2_summary = r2_summary,
             diagnostics_summary = diagnostics_summary,
             model_structure = model_structure, pred_all = pred_all,
             amp_limit = erp$amp_limit),
        path(LMM_DIR, "erp_results.rds"))

message("\nsaved to ", LMM_DIR)

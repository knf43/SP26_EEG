# ============================================================
# RT family check: StressLearning2025
# Fits the primary RT model under three response distributions and
# compares DHARMa diagnostics and AIC, so the family is chosen from the
# data rather than from a citation.
#
# The Gamma-log fit currently in statistical_analysis.R is underdispersed
# (DHARMa dispersion 0.591, p < .001). Underdispersion inflates standard
# errors, so it makes tests conservative rather than anti-conservative,
# but it still means the variance function is wrong.
#
# Run this BEFORE followup_analysis.R, since the four follow-up RT models
# should use whichever family wins here.
# ============================================================

# libraries
library(tidyverse)
library(lme4)
library(DHARMa)
library(fs)

# configuration
OUTPUT_DIR <- "~/Desktop/SP26_data/behavioral_output"
GLMM_DIR   <- path(OUTPUT_DIR, "glmm")

dir_create(GLMM_DIR)
theme_set(theme_minimal(base_size = 12))

ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

# data preparation
source("script.R")
long <- res$long

long <- long |>
  mutate(
    sound_for_item = coalesce(sounds, sounds_stem),
    item = str_extract(sound_for_item, "(?<=[xt]_)\\d{2}")
  )

COND_ORDER <- c("parox_nounmatch_verbmatch", "parox_nounmismatch_verbmatch",
                "parox_nounmatch_verbmismatch", "oxy_nounmatch_verbmatch",
                "oxy_nounmismatch_verbmatch", "oxy_nounmatch_verbmismatch")

long_testing <- long |>
  filter(phase %in% c("pre", "post"), !is.na(cond_name), !is.na(item)) |>
  mutate(
    group       = factor(group, levels = c("L1", "L2")),
    phase       = factor(phase, levels = c("pre", "post")),
    cond_name   = factor(cond_name, levels = COND_ORDER),
    item        = factor(item),
    participant = factor(participant)
  )

contrasts(long_testing$group)     <- contr.sum(2)
contrasts(long_testing$phase)     <- contr.sum(2)
contrasts(long_testing$cond_name) <- contr.sum(6)

rt_dat <- long_testing |>
  filter(!is.na(rt_correct), rt_correct > 0) |>
  mutate(log_rt = log(rt_correct))

message("RT data: ", nrow(rt_dat), " trials, ",
        n_distinct(rt_dat$participant), " participants")

# models
FORM <- ~ group * phase + cond_name + (1 + phase | participant) + (1 | item)

message("\n=== Gamma, log link ===")
m_gamma <- glmer(update(FORM, rt_correct ~ .), data = rt_dat,
                 family = Gamma(link = "log"), control = ctrl)

message("\n=== inverse Gaussian, log link ===")
m_invgauss <- glmer(update(FORM, rt_correct ~ .), data = rt_dat,
                    family = inverse.gaussian(link = "log"), control = ctrl)

message("\n=== lognormal (Gaussian on log RT) ===")
m_lognorm <- lmer(update(FORM, log_rt ~ .), data = rt_dat,
                  control = lmerControl(optimizer = "bobyqa",
                                        optCtrl = list(maxfun = 2e5)))

models <- list(gamma_log = m_gamma,
               invgauss_log = m_invgauss,
               lognormal = m_lognorm)

# model checks
check_one <- function(m, nm) {
  sim <- simulateResiduals(m, n = 500, seed = 1)
  ks   <- testUniformity(sim, plot = FALSE)
  disp <- testDispersion(sim, plot = FALSE)
  png(path(GLMM_DIR, paste0("rtfamily_dharma_", nm, ".png")),
      width = 1000, height = 500)
  plot(sim)
  dev.off()
  plot(sim)
  data.frame(
    family        = nm,
    AIC           = round(AIC(m), 1),
    KS_p          = signif(ks$p.value, 3),
    dispersion    = round(as.numeric(disp$statistic), 3),
    dispersion_p  = signif(disp$p.value, 3),
    singular      = isSingular(m)
  )
}

family_table <- bind_rows(
  lapply(names(models), \(nm) check_one(models[[nm]], nm))
)

cat("\n=== RT family comparison ===\n")
cat("dispersion near 1 is good; KS_p above .05 is good;",
    "AIC is only comparable within the same response scale\n\n")
print(family_table)

cat("\nNote: AIC for lognormal is on log(RT) and is NOT comparable",
    "to the two GLMMs, which are on raw RT. Compare it on the",
    "DHARMa columns only.\n")

# results
# Do the three families agree about the effects? If they do, the family
# choice does not change any conclusion and the simplest defensible one wins.
lrt <- function(m, drop_term, nm) {
  reduced <- update(m, as.formula(paste(". ~ . -", drop_term)))
  cmp <- anova(reduced, m)
  df_col <- if ("Chi Df" %in% names(cmp)) "Chi Df" else "Df"
  data.frame(family = nm, effect = drop_term,
             chisq = round(cmp$Chisq[2], 2),
             df = cmp[[df_col]][2],
             p = signif(cmp$`Pr(>Chisq)`[2], 3))
}

effects_table <- bind_rows(
  lapply(names(models), \(nm)
         bind_rows(lrt(models[[nm]], "group:phase", nm),
                   lrt(models[[nm]], "cond_name",   nm)))
)

cat("\n=== do the families agree about the effects? ===\n")
print(effects_table)

# save
saveRDS(list(models = models,
             family_table = family_table,
             effects_table = effects_table),
        path(GLMM_DIR, "rt_family_check.rds"))

message("\nRT family check saved to: ", GLMM_DIR)

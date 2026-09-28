# ============================================================
# Add confidence intervals to the cost tables
#
# The eight follow-up models are already fitted and saved. This computes the
# match/mismatch contrasts from those stored objects and writes the cost
# tables back into the RDS. NOTHING IS REFITTED.
#
# Run once, then behavioral_plots.R will draw the intervals.
# ============================================================

# libraries
library(tidyverse)
library(lme4)
library(emmeans)
library(fs)

# configuration
GLMM_DIR <- "~/Desktop/SP26_data/behavioral_output/glmm"
RDS      <- path(GLMM_DIR, "glmm_followup_results.rds")

emm_options(lmer.df = "asymptotic")

# data
# emmeans re-evaluates the model's data argument, so dat and dat_rt have to
# exist under those names, built exactly as followup_analysis.R built them.
source("script.R")

long <- res$long |>
  mutate(
    sound_for_item = coalesce(sounds, sounds_stem),
    item = str_extract(sound_for_item, "(?<=[xt]_)\\d{2}")
  )

dat <- long |>
  filter(phase %in% c("pre", "post"),
         !is.na(item), !is.na(verb_match), !is.na(noun_match),
         !is.na(verb_stress)) |>
  mutate(
    group       = factor(group,       levels = c("L1", "L2")),
    phase       = factor(phase,       levels = c("pre", "post")),
    verb_match  = factor(verb_match,  levels = c("Match", "Mismatch")),
    noun_match  = factor(noun_match,  levels = c("Match", "Mismatch")),
    verb_stress = factor(verb_stress, levels = c("paroxytone", "oxytone")),
    item        = factor(item),
    participant = factor(participant)
  )

contrasts(dat$group)       <- contr.sum(2)
contrasts(dat$phase)       <- contr.sum(2)
contrasts(dat$verb_match)  <- contr.sum(2)
contrasts(dat$noun_match)  <- contr.sum(2)
contrasts(dat$verb_stress) <- contr.sum(2)

dat_rt <- filter(dat, !is.na(rt_correct), rt_correct > 0) |>
  mutate(log_rt = log(rt_correct))

followup <- readRDS(RDS)

# helper functions
# regrid() moves the estimates onto the response scale before contrasting, so
# the difference is in probability or seconds rather than log odds or log rt.
# Cost convention: RT = Mismatch - Match (ms); accuracy = Match - Mismatch (pp).
cost_emm <- function(model, factor_name, by_terms, measure) {
  args <- list(object = model, specs = factor_name, by = by_terms)
  if (measure == "rt") args$tran <- "log"
  rg <- regrid(do.call(emmeans, args))

  scale <- if (measure == "rt") 1000 else 100

  # emmeans names the estimate column by response type: prob for a binomial
  # model, response for the log-transformed one. The interval columns are
  # lower.CL/upper.CL with df and asymp.LCL/asymp.UCL without. Both spellings
  # are looked up rather than assumed.
  pick <- function(d, cands) intersect(cands, names(d))[1]

  md  <- as.data.frame(rg)
  est <- pick(md, c("response", "prob", "emmean", "rate"))
  means <- md[, c(factor_name, by_terms, est)]
  names(means)[ncol(means)] <- ".est"
  means$.est <- means$.est * scale
  means <- tidyr::pivot_wider(means, names_from = all_of(factor_name),
                              values_from = ".est")

  ct <- contrast(rg, if (measure == "rt") "revpairwise" else "pairwise")
  cf <- as.data.frame(confint(ct))
  lo <- pick(cf, c("lower.CL", "asymp.LCL"))
  hi <- pick(cf, c("upper.CL", "asymp.UCL"))
  ci <- cf[, c(by_terms, "estimate", lo, hi)]
  names(ci) <- c(by_terms, "cost", "cost_low", "cost_high")
  ci[c("cost", "cost_low", "cost_high")] <-
    ci[c("cost", "cost_low", "cost_high")] * scale

  left_join(means, ci, by = by_terms) |> mutate(Measure = measure)
}

# cost tables
message("computing contrasts, no models are being refitted")

followup$cost_q1q2 <- bind_rows(
  cost_emm(followup$acc_verbmatch, "verb_match", c("group", "phase"), "accuracy") |>
    mutate(Factor = "verb_match"),
  cost_emm(followup$rt_verbmatch,  "verb_match", c("group", "phase"), "rt") |>
    mutate(Factor = "verb_match"),
  cost_emm(followup$acc_nounmatch, "noun_match", c("group", "phase"), "accuracy") |>
    mutate(Factor = "noun_match"),
  cost_emm(followup$rt_nounmatch,  "noun_match", c("group", "phase"), "rt") |>
    mutate(Factor = "noun_match")
) |>
  rename(Group = group, Phase = phase) |>
  select(Factor, Measure, Group, Phase, Match, Mismatch,
         cost, cost_low, cost_high)

followup$cost_q3 <- bind_rows(
  cost_emm(followup$acc_stress_vm, "verb_match", "verb_stress", "accuracy") |>
    mutate(Factor = "verb_match"),
  cost_emm(followup$rt_stress_vm,  "verb_match", "verb_stress", "rt") |>
    mutate(Factor = "verb_match")
) |>
  rename(Stress = verb_stress) |>
  select(Factor, Measure, Stress, Match, Mismatch, cost, cost_low, cost_high)

followup$cost_q4 <- bind_rows(
  cost_emm(followup$acc_g_s_vm, "verb_match", c("group", "verb_stress"), "accuracy") |>
    mutate(Factor = "verb_match"),
  cost_emm(followup$rt_g_s_vm,  "verb_match", c("group", "verb_stress"), "rt") |>
    mutate(Factor = "verb_match")
) |>
  rename(Group = group, Stress = verb_stress) |>
  select(Factor, Measure, Group, Stress, Match, Mismatch,
         cost, cost_low, cost_high)

print(followup$cost_q1q2)
print(followup$cost_q4)

# save
saveRDS(followup, RDS)
message("\ncost tables now carry 95% confidence intervals")

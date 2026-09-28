# ============================================================
# Follow-up analysis: StressLearning2025
# Verb-mismatch and noun-mismatch effects on accuracy and RT
# Same conventions as statistical_analysis.R: lme4, sum contrasts,
# nested model comparisons, lognormal for RT, binomial-logit for accuracy
#
# Run AFTER statistical_analysis.R. Writes glmm_followup_results.rds
# and the followup_*.png figures that behavioral.Rmd reads.
# ============================================================

# libraries
library(tidyverse)
library(lme4)
library(ggeffects)
library(emmeans)
library(fs)

# configuration
DATA_DIR   <- "~/Desktop/SP26_data/2_psychopy_data"
OUTPUT_DIR <- "~/Desktop/SP26_data/behavioral_output"
GLMM_DIR   <- path(OUTPUT_DIR, "glmm")

dir_create(GLMM_DIR)

ctrl_bin <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

theme_set(theme_minimal(base_size = 12))

# Kenward-Roger df on 20,000 rows is impractically slow and the asymptotic
# approximation is what the nested model comparisons already assume.
emm_options(lmer.df = "asymptotic")

# data preparation
# Mirrors section 1 of statistical_analysis.R so this script runs on its own.
source("script.R")
long <- res$long

long <- long |>
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

# Lognormal for RT, chosen by rt_family_check.R. exp() of a fitted value is
# therefore a geometric mean, which is what the cost tables report.
dat_rt <- filter(dat, !is.na(rt_correct), rt_correct > 0) |>
  mutate(log_rt = log(rt_correct))
ctrl_lmm <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

message("Follow-up data: ", nrow(dat), " accuracy trials, ",
        nrow(dat_rt), " RT trials, ",
        n_distinct(dat$participant), " participants, ",
        n_distinct(dat$item), " items")

# helper functions
# Every reduced model is built with update() rather than a fresh glmer() call.
# anova.merMod compares models by the deparsed name of their data argument, so
# refitting with `data = <local variable>` makes lme4 refuse with
# "all models must be fit to the same data object". update() keeps the original
# symbol, which is why each model below is fit with `dat` or `dat_rt` written
# out literally rather than passed into a wrapper function.
# Reduced models are built by rewriting the formula from its own term labels
# rather than with `. ~ . - term` arithmetic. Subtracting a two-way term from a
# model that also contains the three-way interaction silently left the model
# unchanged, which showed up as chisq 0, df 0, p NA. Rebuilding from term
# labels removes exactly the named term, and stops with an error if that term
# is not in the model rather than reporting a null result.
drop_term_formula <- function(model, drop_term) {
  f     <- formula(model)
  tl    <- attr(terms(lme4::nobars(f)), "term.labels")
  if (!drop_term %in% tl) {
    stop("term not in model: ", drop_term, "\n  available: ",
         paste(tl, collapse = ", "))
  }
  bars  <- vapply(lme4::findbars(f),
                  function(b) paste0("(", deparse(b), ")"), character(1))
  stats::as.formula(
    paste(deparse(f[[2]]), "~",
          paste(c(setdiff(tl, drop_term), bars), collapse = " + ")),
    env = environment(f)
  )
}

nmc_compare <- function(model, drop_term, label) {
  reduced <- update(model, drop_term_formula(model, drop_term))
  cmp <- anova(reduced, model)
  df_col <- if ("Chi Df" %in% names(cmp)) "Chi Df" else "Df"
  if (isTRUE(cmp[[df_col]][2] == 0)) {
    stop("dropping ", drop_term, " from ", label, " changed nothing")
  }
  data.frame(
    model  = label,
    effect = drop_term,
    chisq  = round(cmp$Chisq[2], 2),
    df     = cmp[[df_col]][2],
    p      = signif(cmp$`Pr(>Chisq)`[2], 3)
  )
}

# Main effects involved in an interaction are tested against a no-interaction
# model, so marginality is preserved.
nmc_main <- function(model, drop_term, label, interactions) {
  base <- model
  for (ix in interactions) base <- update(base, drop_term_formula(base, ix))
  nmc_compare(base, drop_term, label)
}

# Model-derived match/mismatch values and the cost between them, with a 95%
# confidence interval on the cost itself. The interval comes from an emmeans
# contrast rather than from subtracting two predictions, because the difference
# of two intervals is not the interval of the difference.
#
# regrid() moves the estimates onto the response scale first, so the contrast
# is in probability or seconds rather than log odds or log rt.
#
# Cost convention, unchanged: RT cost = Mismatch - Match (ms);
# accuracy cost = Match - Mismatch (percentage points).
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

save_and_print <- function(p, filename, w = 8, h = 4.5) {
  print(p)
  ggsave(path(GLMM_DIR, filename), p, width = w, height = h, dpi = 150)
  invisible(p)
}

# models
message("\n=== verb match x group x phase ===")
acc_verbmatch <- glmer(corr ~ group * phase * verb_match +
                         (1 + phase | participant) + (1 | item),
                       data = dat, family = binomial("logit"), control = ctrl_bin)
rt_verbmatch <- lmer(log_rt ~ group * phase * verb_match +
                        (1 + phase | participant) + (1 | item),
                 data = dat_rt, control = ctrl_lmm)

message("\n=== noun match x group x phase ===")
acc_nounmatch <- glmer(corr ~ group * phase * noun_match +
                         (1 + phase | participant) + (1 | item),
                       data = dat, family = binomial("logit"), control = ctrl_bin)
rt_nounmatch <- lmer(log_rt ~ group * phase * noun_match +
                        (1 + phase | participant) + (1 | item),
                 data = dat_rt, control = ctrl_lmm)

message("\n=== verb stress x verb match ===")
acc_stress_vm <- glmer(corr ~ verb_stress * verb_match +
                         (1 + verb_match || participant) + (1 | item),
                       data = dat, family = binomial("logit"), control = ctrl_bin)
rt_stress_vm <- lmer(log_rt ~ verb_stress * verb_match +
                        (1 + verb_match || participant) + (1 | item),
                 data = dat_rt, control = ctrl_lmm)

message("\n=== group x verb stress x verb match ===")
acc_g_s_vm <- glmer(corr ~ group * verb_stress * verb_match +
                      (1 + verb_match || participant) + (1 | item),
                    data = dat, family = binomial("logit"), control = ctrl_bin)
rt_g_s_vm <- lmer(log_rt ~ group * verb_stress * verb_match +
                     (1 + verb_match || participant) + (1 | item),
                 data = dat_rt, control = ctrl_lmm)

# model checks
singular_table <- data.frame(
  model = c("acc_verbmatch", "rt_verbmatch", "acc_nounmatch", "rt_nounmatch",
            "acc_stress_vm", "rt_stress_vm", "acc_g_s_vm", "rt_g_s_vm"),
  singular = vapply(list(acc_verbmatch, rt_verbmatch, acc_nounmatch, rt_nounmatch,
                         acc_stress_vm, rt_stress_vm, acc_g_s_vm, rt_g_s_vm),
                    isSingular, logical(1))
)
cat("\n=== singular fits ===\n")
print(singular_table)

# nested model comparisons
message("\n=== nested model comparisons ===")

# A two-way term cannot be dropped while the three-way that contains it is
# still in the model: model.matrix re-absorbs it into the higher-order term and
# df comes back 0. Each two-way is therefore tested against the model with the
# three-way removed, the same marginality rule already applied to main effects.
nmc_followup <- bind_rows(
  nmc_compare(acc_verbmatch, "group:phase:verb_match", "acc_verbmatch"),
  nmc_main(acc_verbmatch, "group:verb_match", "acc_verbmatch", "group:phase:verb_match"),
  nmc_main(acc_verbmatch, "phase:verb_match", "acc_verbmatch", "group:phase:verb_match"),

  nmc_compare(rt_verbmatch, "group:phase:verb_match", "rt_verbmatch"),
  nmc_main(rt_verbmatch, "group:verb_match", "rt_verbmatch", "group:phase:verb_match"),
  nmc_main(rt_verbmatch, "phase:verb_match", "rt_verbmatch", "group:phase:verb_match"),

  nmc_compare(acc_nounmatch, "group:phase:noun_match", "acc_nounmatch"),
  nmc_main(acc_nounmatch, "group:noun_match", "acc_nounmatch", "group:phase:noun_match"),
  nmc_main(acc_nounmatch, "phase:noun_match", "acc_nounmatch", "group:phase:noun_match"),

  nmc_compare(rt_nounmatch, "group:phase:noun_match", "rt_nounmatch"),
  nmc_main(rt_nounmatch, "group:noun_match", "rt_nounmatch", "group:phase:noun_match"),
  nmc_main(rt_nounmatch, "phase:noun_match", "rt_nounmatch", "group:phase:noun_match"),

  nmc_compare(acc_stress_vm, "verb_stress:verb_match", "acc_stress_vm"),
  nmc_main(acc_stress_vm, "verb_stress", "acc_stress_vm", "verb_stress:verb_match"),
  nmc_main(acc_stress_vm, "verb_match",  "acc_stress_vm", "verb_stress:verb_match"),

  nmc_compare(rt_stress_vm, "verb_stress:verb_match", "rt_stress_vm"),
  nmc_main(rt_stress_vm, "verb_stress", "rt_stress_vm", "verb_stress:verb_match"),
  nmc_main(rt_stress_vm, "verb_match",  "rt_stress_vm", "verb_stress:verb_match"),

  nmc_compare(acc_g_s_vm, "group:verb_stress:verb_match", "acc_g_s_vm"),
  nmc_main(acc_g_s_vm, "group:verb_stress",      "acc_g_s_vm", "group:verb_stress:verb_match"),
  nmc_main(acc_g_s_vm, "group:verb_match",       "acc_g_s_vm", "group:verb_stress:verb_match"),
  nmc_main(acc_g_s_vm, "verb_stress:verb_match", "acc_g_s_vm", "group:verb_stress:verb_match"),

  nmc_compare(rt_g_s_vm, "group:verb_stress:verb_match", "rt_g_s_vm"),
  nmc_main(rt_g_s_vm, "group:verb_stress",      "rt_g_s_vm", "group:verb_stress:verb_match"),
  nmc_main(rt_g_s_vm, "group:verb_match",       "rt_g_s_vm", "group:verb_stress:verb_match"),
  nmc_main(rt_g_s_vm, "verb_stress:verb_match", "rt_g_s_vm", "group:verb_stress:verb_match")
)

print(nmc_followup)

# model-derived predictions
pred_acc_vm <- ggpredict(acc_verbmatch, terms = c("verb_match", "group", "phase"), bias_correction = TRUE)
pred_rt_vm  <- ggpredict(rt_verbmatch,  terms = c("verb_match", "group", "phase"))
pred_acc_nm <- ggpredict(acc_nounmatch, terms = c("noun_match", "group", "phase"), bias_correction = TRUE)
pred_rt_nm  <- ggpredict(rt_nounmatch,  terms = c("noun_match", "group", "phase"))
pred_acc_sv <- ggpredict(acc_stress_vm, terms = c("verb_match", "verb_stress"), bias_correction = TRUE)
pred_rt_sv  <- ggpredict(rt_stress_vm,  terms = c("verb_match", "verb_stress"))
pred_acc_gsv <- ggpredict(acc_g_s_vm, terms = c("verb_match", "group", "verb_stress"), bias_correction = TRUE)
pred_rt_gsv  <- ggpredict(rt_g_s_vm,  terms = c("verb_match", "group", "verb_stress"))

# cost tables
cost_q1q2 <- bind_rows(
  cost_emm(acc_verbmatch, "verb_match", c("group", "phase"), "accuracy") |>
    mutate(Factor = "verb_match"),
  cost_emm(rt_verbmatch,  "verb_match", c("group", "phase"), "rt") |>
    mutate(Factor = "verb_match"),
  cost_emm(acc_nounmatch, "noun_match", c("group", "phase"), "accuracy") |>
    mutate(Factor = "noun_match"),
  cost_emm(rt_nounmatch,  "noun_match", c("group", "phase"), "rt") |>
    mutate(Factor = "noun_match")
) |>
  rename(Group = group, Phase = phase) |>
  select(Factor, Measure, Group, Phase, Match, Mismatch,
         cost, cost_low, cost_high)

cost_q3 <- bind_rows(
  cost_emm(acc_stress_vm, "verb_match", "verb_stress", "accuracy") |>
    mutate(Factor = "verb_match"),
  cost_emm(rt_stress_vm,  "verb_match", "verb_stress", "rt") |>
    mutate(Factor = "verb_match")
) |>
  rename(Stress = verb_stress) |>
  select(Factor, Measure, Stress, Match, Mismatch, cost, cost_low, cost_high)

cost_q4 <- bind_rows(
  cost_emm(acc_g_s_vm, "verb_match", c("group", "verb_stress"), "accuracy") |>
    mutate(Factor = "verb_match"),
  cost_emm(rt_g_s_vm,  "verb_match", c("group", "verb_stress"), "rt") |>
    mutate(Factor = "verb_match")
) |>
  rename(Group = group, Stress = verb_stress) |>
  select(Factor, Measure, Group, Stress, Match, Mismatch,
         cost, cost_low, cost_high)

print(cost_q1q2)
print(cost_q3)
print(cost_q4)

# figures
pred_plot <- function(pred, ylab, pct = FALSE, expo = FALSE) {
  d <- as.data.frame(pred)
  if (expo) d <- mutate(d, across(c(predicted, conf.low, conf.high), exp))
  p <- ggplot(d, aes(x = x, y = predicted, colour = group, group = group)) +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                    position = position_dodge(0.25), size = 0.6) +
    geom_line(position = position_dodge(0.25)) +
    labs(x = NULL, y = ylab, colour = "Group")
  if ("facet" %in% names(d)) p <- p + facet_wrap(~ facet)
  if (pct) p <- p + scale_y_continuous(labels = scales::percent)
  p
}

cost_plot <- function(d, xvar, fillvar, ylab, title) {
  ggplot(d, aes(x = .data[[xvar]], y = cost, fill = .data[[fillvar]])) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    geom_hline(yintercept = 0, colour = "grey40") +
    labs(x = NULL, y = ylab, fill = NULL, title = title)
}

save_and_print(pred_plot(pred_rt_vm,  "Predicted RT (s, geometric mean)", expo = TRUE),
               "followup_rt_verbmatch.png")
save_and_print(pred_plot(pred_acc_vm, "Predicted accuracy", pct = TRUE),
               "followup_acc_verbmatch.png")
save_and_print(pred_plot(pred_rt_nm,  "Predicted RT (s, geometric mean)", expo = TRUE),
               "followup_rt_nounmatch.png")
save_and_print(pred_plot(pred_acc_nm, "Predicted accuracy", pct = TRUE),
               "followup_acc_nounmatch.png")

vm_cost <- filter(cost_q1q2, Factor == "verb_match")
save_and_print(cost_plot(filter(vm_cost, Measure == "rt"), "Phase", "Group",
                         "RT cost (ms)", "Verb-mismatch RT cost"),
               "followup_cost_rt_phase.png", w = 6)
save_and_print(cost_plot(filter(vm_cost, Measure == "accuracy"), "Phase", "Group",
                         "Accuracy cost (pp)", "Verb-mismatch accuracy cost"),
               "followup_cost_acc_phase.png", w = 6)

save_and_print(cost_plot(filter(cost_q3, Measure == "rt"), "Stress", "Stress",
                         "RT cost (ms)", "Verb-mismatch RT cost by stress"),
               "followup_cost_rt_stress.png", w = 6)
save_and_print(cost_plot(filter(cost_q3, Measure == "accuracy"), "Stress", "Stress",
                         "Accuracy cost (pp)", "Verb-mismatch accuracy cost by stress"),
               "followup_cost_acc_stress.png", w = 6)

save_and_print(cost_plot(filter(cost_q4, Measure == "rt"), "Stress", "Group",
                         "RT cost (ms)", "Verb-mismatch RT cost by group and stress"),
               "followup_cost_rt_group_stress.png", w = 6)
save_and_print(cost_plot(filter(cost_q4, Measure == "accuracy"), "Stress", "Group",
                         "Accuracy cost (pp)",
                         "Verb-mismatch accuracy cost by group and stress"),
               "followup_cost_acc_group_stress.png", w = 6)

# save
saveRDS(list(
  acc_verbmatch = acc_verbmatch, rt_verbmatch = rt_verbmatch,
  acc_nounmatch = acc_nounmatch, rt_nounmatch = rt_nounmatch,
  acc_stress_vm = acc_stress_vm, rt_stress_vm = rt_stress_vm,
  acc_g_s_vm    = acc_g_s_vm,    rt_g_s_vm    = rt_g_s_vm,
  nmc_followup  = nmc_followup,
  singular_table = singular_table,
  cost_q1q2     = cost_q1q2,
  cost_q3       = cost_q3,
  cost_q4       = cost_q4
), path(GLMM_DIR, "glmm_followup_results.rds"))

message("\nFollow-up outputs saved to: ", GLMM_DIR)

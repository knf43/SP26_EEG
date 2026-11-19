# erp permutation tests on difference scores

# libraries
library(tidyverse)
library(fs)

# configuration
ANOVA_DIR    <- "~/Desktop/SP26_data/erp_output/anova"
PERMUTATIONS <- 10000
SEED         <- 26

# the difference score is first minus second. where the anova has a second
# condition factor, the interaction score is the difference score at
# other_first minus the difference score at other_second.
EFFECTS <- tribble(
  ~model,      ~effect,       ~first,     ~second,      ~other,        ~other_first, ~other_second,
  "LAN_verb",  "verb_match",  "Mismatch", "Match",      "verb_stress", "oxytone",    "paroxytone",
  "N400_verb", "verb_match",  "Mismatch", "Match",      "verb_stress", "oxytone",    "paroxytone",
  "P600_verb", "verb_match",  "Mismatch", "Match",      "verb_stress", "oxytone",    "paroxytone",
  "PrAN_verb", "verb_stress", "oxytone",  "paroxytone", "noun_match",  "Match",      "Mismatch",
  "PrAN_noun", "splice",      "c",        "a",          NA,            NA,           NA,
  "LAN_noun",  "noun_match",  "Mismatch", "Match",      "verb_stress", "oxytone",    "paroxytone",
  "N400_noun", "noun_match",  "Mismatch", "Match",      "verb_stress", "oxytone",    "paroxytone",
  "P600_noun", "noun_match",  "Mismatch", "Match",      "verb_stress", "oxytone",    "paroxytone"
)

WITHIN_ORDER <- c("verb_stress", "verb_match", "noun_match", "splice")

ANOVA_TERMS <- c("all participants"                        = "%s",
                 "L1 against L2"                           = "group:%s",
                 "post against pre"                        = "phase:%s",
                 "change from pre to post, L1 against L2"  = "group:phase:%s")

MEAN_TERMS <- c("L1 against L2"                           = "group",
                "post against pre"                        = "phase",
                "change from pre to post, L1 against L2"  = "group:phase")

# data
cell_means <- read_csv(path(ANOVA_DIR, "anova_cell_means.csv"), show_col_types = FALSE)
anova      <- read_csv(path(ANOVA_DIR, "anova_table.csv"), show_col_types = FALSE)

# one difference score, one mean and one interaction score per participant
# and phase
scores <- pmap_dfr(EFFECTS, \(model, effect, first, second, other, other_first, other_second) {
  rows <- filter(cell_means, .data$model == .env$model)

  main <- rows |>
    group_by(participant, group, phase, level = .data[[effect]]) |>
    summarise(amplitude = mean(amplitude), .groups = "drop") |>
    pivot_wider(names_from = level, values_from = amplitude) |>
    transmute(participant, group, phase,
              difference = .data[[first]] - .data[[second]],
              mean       = (.data[[first]] + .data[[second]]) / 2)

  if (!is.na(other)) {
    main <- rows |>
      transmute(participant, phase, level = .data[[effect]],
                other_level = .data[[other]], amplitude) |>
      pivot_wider(names_from = level, values_from = amplitude) |>
      mutate(difference = .data[[first]] - .data[[second]]) |>
      select(participant, phase, other_level, difference) |>
      pivot_wider(names_from = other_level, values_from = difference) |>
      transmute(participant, phase,
                interaction = .data[[other_first]] - .data[[other_second]]) |>
      right_join(main, by = c("participant", "phase"))
  }

  factors <- intersect(WITHIN_ORDER, c(effect, other))
  mutate(main, model, effect, interaction_term = paste(factors, collapse = ":"),
         .before = 1)
})

# helper functions
t_one <- function(x) mean(x) / (sd(x) / sqrt(length(x)))

# the mean of the two group means against zero, as in the anova
t_both <- function(x, group) {
  n      <- table(group)
  means  <- tapply(x, group, mean)
  pooled <- sum(tapply(x, group, \(v) sum((v - mean(v))^2))) / (length(x) - 2)
  mean(means) / (sqrt(pooled * sum(1 / n)) / 2)
}

t_two <- function(x, group) {
  unname(t.test(x[group == "L1"], x[group == "L2"], var.equal = TRUE)$statistic)
}

permutation_p <- function(observed, permuted) {
  (sum(abs(permuted) >= abs(observed)) + 1) / (length(permuted) + 1)
}

# is the score different from zero? the sign of each participant's score is
# flipped at random
test_zero <- function(x) {
  observed <- t_one(x)
  permuted <- replicate(PERMUTATIONS, t_one(x * sample(c(-1, 1), length(x), replace = TRUE)))
  tibble(estimate = mean(x), t = observed, df = length(x) - 1,
         p_parametric = 2 * pt(-abs(observed), length(x) - 1),
         p_permutation = permutation_p(observed, permuted))
}

# the same question for both groups together, with the groups weighted equally
test_both <- function(x, group) {
  observed <- t_both(x, group)
  permuted <- replicate(PERMUTATIONS,
                        t_both(x * sample(c(-1, 1), length(x), replace = TRUE), group))
  tibble(estimate = mean(tapply(x, group, mean)), t = observed, df = length(x) - 2,
         p_parametric = 2 * pt(-abs(observed), length(x) - 2),
         p_permutation = permutation_p(observed, permuted))
}

# does the score differ between the groups? the group labels are shuffled
test_groups <- function(x, group) {
  observed <- t_two(x, group)
  permuted <- replicate(PERMUTATIONS, t_two(x, sample(group)))
  tibble(estimate = mean(x[group == "L1"]) - mean(x[group == "L2"]),
         t = observed, df = length(x) - 2,
         p_parametric = 2 * pt(-abs(observed), length(x) - 2),
         p_permutation = permutation_p(observed, permuted))
}

test_score <- function(data, score) {
  wide <- data |>
    select(participant, group, phase, value = all_of(score)) |>
    pivot_wider(names_from = phase, values_from = value) |>
    mutate(average = (pre + post) / 2, change = post - pre)

  cells <- data |>
    group_by(group, phase) |>
    group_modify(\(rows, key) test_zero(rows[[score]])) |>
    ungroup() |>
    mutate(test = paste(group, "at", phase),
           p_parametric  = p.adjust(p_parametric, "fdr"),
           p_permutation = p.adjust(p_permutation, "fdr")) |>
    arrange(group, desc(phase)) |>
    select(-group, -phase)

  results <- bind_rows(
    mutate(test_both(wide$average, wide$group),   test = "all participants"),
    mutate(test_groups(wide$average, wide$group), test = "L1 against L2"),
    mutate(test_both(wide$change, wide$group),    test = "post against pre"),
    mutate(test_groups(wide$change, wide$group),  test = "change from pre to post, L1 against L2"),
    cells) |>
    mutate(score = score, .before = 1)

  # the mean amplitude against zero is not a question of the study
  if (score == "mean") results <- filter(results, test %in% names(MEAN_TERMS))
  results
}

# tests
set.seed(SEED)

SCORES <- c("difference", "interaction", "mean")

results <- scores |>
  group_by(model, effect, interaction_term) |>
  group_modify(\(data, key) {
    available <- SCORES[map_lgl(SCORES, \(score) score %in% names(data) &&
                                  !all(is.na(data[[score]])))]
    map_dfr(available, \(score) test_score(data, score))
  }) |>
  ungroup() |>
  mutate(model = factor(model, levels = EFFECTS$model),
         score = factor(score, levels = SCORES)) |>
  arrange(model, score) |>
  relocate(test, .after = score)

# the same effect in the anova, for comparison
results <- results |>
  mutate(anova_effect = case_when(
    score == "difference"  ~ sprintf(ANOVA_TERMS[test], effect),
    score == "interaction" ~ sprintf(ANOVA_TERMS[test], interaction_term),
    score == "mean"        ~ unname(MEAN_TERMS[test]))) |>
  left_join(select(anova, model, anova_effect = effect, anova_F = F, anova_p = p),
            by = c("model", "anova_effect")) |>
  mutate(across(c(estimate, t), \(x) round(x, 2)),
         across(c(p_parametric, p_permutation), \(x) round(x, 4)))

# results
cat("\npermutation tests on difference scores\n")
print(filter(results, score == "difference"), n = Inf, width = Inf)

cat("\npermutation tests on the interaction scores\n")
print(filter(results, score == "interaction"), n = Inf, width = Inf)

cat("\npermutation tests on the mean of the two conditions\n")
print(filter(results, score == "mean"), n = Inf, width = Inf)

# save
write_csv(results, path(ANOVA_DIR, "permutation_tests.csv"))
write_csv(scores,  path(ANOVA_DIR, "permutation_scores.csv"))

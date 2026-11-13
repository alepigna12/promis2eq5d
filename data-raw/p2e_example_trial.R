# Builds data/p2e_example_trial.rda.
# SYNTHETIC data: respondents resampled from the CC0 p2e_training data with a
# made-up treatment effect. This is not a real trial.
# Usage, from the package directory: Rscript data-raw/p2e_example_trial.R
#
# Draw order (deterministic, set.seed(20261113)):
#   1. sample() of 300 distinct UK rows (first 150 control, rest treatment);
#   2. month 6: one runif() per person x item, items in ITEMS order then P9
#      (a 300 x 29 matrix, filled column by column);
#   3. month 12: the same, from the baseline values (not from month 6).
# Per item: U < p_better -> one step better; U < p_better + 0.05 -> one step
# worse; otherwise unchanged. Clamped to 1..5 (items) and 0..10 (P9).
# p_better (made up): control 0.05 at months 6 and 12; treatment 0.35 at
# month 6 and 0.50 at month 12.
# Better direction: +1 for P1_* and P7_*; -1 for all other items and for P9.

if (!file.exists("DESCRIPTION") || !dir.exists("data-raw")) stop("run from the package directory")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
set.seed(20261113)

uk <- p2e_training[p2e_training$country == "UK", ]
base <- uk[sample.int(nrow(uk), 300L), ]
arm <- rep(c("control", "treatment"), each = 150L)
vars <- c(ITEMS, "P9")
up <- ifelse(grepl("^P1_|^P7_", vars), 1L, -1L)
lo <- c(rep(1L, length(ITEMS)), 0L); hi <- c(rep(5L, length(ITEMS)), 10L)
p_better <- list(`6` = c(control = 0.05, treatment = 0.35),
                 `12` = c(control = 0.05, treatment = 0.50))

B <- as.matrix(base[vars])
follow <- function(month) {
  U <- matrix(runif(nrow(B) * length(vars)), nrow(B))
  pb <- matrix(p_better[[as.character(month)]][arm], nrow(B), length(vars))
  step <- ifelse(U < pb, 1L, ifelse(U < pb + 0.05, -1L, 0L))
  X <- B + sweep(step, 2, up, `*`)
  X <- pmin(pmax(X, rep(lo, each = nrow(B))), rep(hi, each = nrow(B)))
  storage.mode(X) <- "integer"; X
}
waves <- list(`0` = B, `6` = follow(6), `12` = follow(12))

one <- function(month) {
  X <- waves[[as.character(month)]]
  d <- data.frame(id = 1:300, arm = arm, month = as.integer(month), stringsAsFactors = FALSE)
  for (j in seq_along(ITEMS)) d[[ITEM_IDS[[ITEMS[j]]]]] <- X[, ITEMS[j]]
  d[[PAIN_ID]] <- X[, "P9"]
  d$age <- base$AGE + month / 12
  d$female <- as.integer(base$female)
  d
}
p2e_example_trial <- do.call(rbind, lapply(c(0L, 6L, 12L), one))
p2e_example_trial <- p2e_example_trial[order(p2e_example_trial$id, p2e_example_trial$month), ]
rownames(p2e_example_trial) <- NULL

save(p2e_example_trial, file = "data/p2e_example_trial.rda", compress = "xz")

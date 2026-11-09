# Ported from iHEA_2027 R/functions/analysis_data_dictionary.R at freeze-ihea2027-abstract-v3 (a79e102); changes: dropped DIAG, S1D/S3D, MODEL_DATA_COLS and VT_COUNTRIES; added ITEM_IDS, PAIN_ID, METHOD_INFO and the Monte Carlo constants

# EQ-5D-5L dimensions (levels 1-5): E1 = MO, E2 = SC, E3 = UA, E4 = PD, E5 = AD
DIMS <- c("E1", "E2", "E3", "E4", "E5")

# The 28 PROMIS-29 items (PROMIS scoring orientation, integers 1-5). P9 (pain
# intensity, 0-10) is kept separate from ITEMS.
ITEMS <- c(paste0("P1_", 1:4), paste0("P2_", 1:4), paste0("P3_", 1:4),
           paste0("P4_", 1:4), "P5", paste0("P6_", 1:3),
           paste0("P7_", 1:4), paste0("P8_", 1:4))

# HealthMeasures item IDs (the user-facing item column names), in ITEMS order.
# Verified against the iHEA_2027 data dictionary and the promis.R header.
ITEM_IDS <- stats::setNames(
  c("PFA11", "PFA21", "PFA23", "PFA53",
    "EDANX01", "EDANX40", "EDANX41", "EDANX53",
    "EDDEP04", "EDDEP06", "EDDEP29", "EDDEP41",
    "HI7", "AN3", "FATEXP41", "FATEXP40",
    "Sleep109", "Sleep116", "Sleep20", "Sleep44",
    "SRPPER11_CaPS", "SRPPER18_CaPS", "SRPPER23_CaPS", "SRPPER46_CaPS",
    "PAININ9", "PAININ22", "PAININ31", "PAININ34"),
  ITEMS)

# User-facing name of the pain intensity item (alias P9)
PAIN_ID <- "Global07"

# PROMIS T-scores: physical function, anxiety, depression, fatigue,
# sleep disturbance, social roles, pain interference
TSCORES <- c("T_PF", "T_ANX", "T_DEP", "T_FAT", "T_SLP", "T_SR", "T_PI")

# Covariates
COVS <- c("AGE", "female")

# Predictor sets. S3K: Klapproth et al. (2020), T-scores, age and sex without
# P9 (D12 only).
PREDICTOR_SETS <- list(
  S1 = c(ITEMS, "P9", COVS),
  S2 = c(ITEMS, "P9"),
  S3 = c(TSCORES, "P9", COVS),
  S4 = c(TSCORES, "P9"),
  S3K = c(TSCORES, COVS)
)

# Conceptually motivated predictors per EQ-5D dimension, as item-based and
# T-score-based versions. For S1/S3-type sets, COVS are added to every
# dimension.
CONCEPTUAL_SUBSETS <- list(
  E1 = list(items   = c(paste0("P1_", 1:4), paste0("P8_", 1:4), "P9"),
            tscores = c("T_PF", "T_PI", "P9")),
  E2 = list(items   = paste0("P1_", 1:4),
            tscores = "T_PF"),
  E3 = list(items   = c(paste0("P7_", 1:4), paste0("P1_", 1:4),
                        paste0("P4_", 1:4), paste0("P8_", 1:4)),
            tscores = c("T_SR", "T_PF", "T_FAT", "T_PI")),
  E4 = list(items   = c(paste0("P8_", 1:4), "P9"),
            tscores = c("T_PI", "P9")),
  E5 = list(items   = c(paste0("P2_", 1:4), paste0("P3_", 1:4)),
            tscores = c("T_ANX", "T_DEP"))
)

# Methods: family and preferred predictor-set order
METHOD_INFO <- list(
  R06_mvop_dim  = list(family = "response", sets = c("S1", "S3", "S2", "S4")),
  R09_chain     = list(family = "response", sets = c("S3", "S1", "S4", "S2")),
  D03_clad      = list(family = "direct",   sets = c("S3", "S4", "S1", "S2")),
  D12_klapproth = list(family = "direct",   sets = "S3K")
)

# Monte Carlo settings, defined once
N_MC_MEAN <- 2000L          # draws per respondent for Monte Carlo means
MC_SEED <- 20261108L        # fixed seed of every Monte Carlo mean
COHORT_LEVEL <- 0.95        # default cohort summary interval level
R09_EXACT_MEAN <- FALSE     # R09 mean by Monte Carlo (FALSE) or exact sum

# Runs expr under the fixed Monte Carlo seed with pinned RNG kinds, so a Monte
# Carlo mean does not depend on the caller's RNGkind() and the caller's stream
# is left alone. Used by every Monte Carlo mean.
.with_mc_seed <- function(expr) {
  withr::with_seed(MC_SEED, expr, .rng_kind = "Mersenne-Twister",
                   .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
}

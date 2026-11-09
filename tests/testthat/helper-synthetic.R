# Ported from iHEA_2027 R/functions/synthetic.R at freeze-ihea2027-abstract-v3 (a79e102); changes: test helper; scoring with .score() and p2e_value_set() per country instead of score_profiles(); countries checked against FR/DE/UK; DIAG and MODEL_DATA_COLS defined locally; no pipeline sourcing or library() calls
#
# Synthetic data with the columns of the pipeline analysis table. Data come from
# a latent 2-factor model (physical and mental health), so the items, T-scores,
# P9 and E1-E5 are positively dependent. The relationships are plausible but not
# calibrated to the real survey. The T-scores are generated here from the items
# (standardised domain sums); no conversion table is read or hard-coded.

.SYN_COUNTRIES <- c("FR", "DE", "UK")
.SYN_DIAG <- c(paste0("Q12_", 1:12), "Q12_15", "Q12_16", "Q12_other")
.SYN_COLS <- c("ID", "country", "AGE", "female", "WEIGHT", "util_band",
               ITEMS, "P9", TSCORES, .SYN_DIAG, DIMS, "U_own")

make_synthetic_data <- function(n = 400, seed = 1, countries = c("DE", "UK")) {
  if (!all(countries %in% .SYN_COUNTRIES)) {
    stop("countries must be a subset of: ", paste(.SYN_COUNTRIES, collapse = ", "))
  }

  # the RNG kind is pinned, so the data never depend on the caller's RNGkind()
  withr::with_seed(seed, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
                   .rng_sample_kind = "Rejection", {
    country <- sample(countries, n, replace = TRUE)
    AGE     <- as.integer(pmin(pmax(round(rnorm(n, 48, 16)), 18), 85))
    female  <- rbinom(n, 1, 0.5)

    # Latent severity (higher = worse health): physical and mental, r = 0.4
    z <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, 0.4, 0.4, 1), 2))
    phys <- z[, 1] + 0.01 * (AGE - 48)
    ment <- z[, 2]

    # Ordinal 1-5 from a latent score; thresholds fixed so that the lower
    # levels dominate, as in general-population samples
    to_level <- function(x, cuts) as.integer(findInterval(x, cuts) + 1L)
    item_cuts <- c(0.2, 0.9, 1.6, 2.3)
    make_item <- function(lat, higher_better = FALSE) {
      lvl <- to_level(lat + rnorm(n, 0, 0.6), item_cuts)
      if (higher_better) 6L - lvl else lvl
    }
    # domain -> latent severity
    lat <- list(
      P1 = phys, P2 = ment, P3 = ment, P4 = 0.5 * phys + 0.5 * ment,
      P5 = 0.5 * phys + 0.5 * ment, P6 = 0.5 * phys + 0.5 * ment,
      P7 = 0.7 * phys + 0.3 * ment, P8 = 0.8 * phys + 0.2 * ment
    )
    # P1 (physical function) and P7 (social roles) are scored higher = better
    better <- c(P1 = TRUE, P7 = TRUE)

    items <- sapply(ITEMS, function(v) {
      dom <- sub("_.*$", "", v)
      make_item(lat[[dom]], isTRUE(better[dom]))
    })
    items <- as.data.frame(items)

    # P9: pain intensity 0-10
    P9 <- as.integer(pmin(pmax(round(2.2 + 1.8 * (0.7 * phys + 0.3 * ment) +
                                       rnorm(n, 0, 1.5)), 0), 10))

    # T-scores: standardised domain item sums, mean 50 and SD 10; high = more
    # of the concept named (for T_PF and T_SR, high = better, as for P1 and P7)
    tscore <- function(prefixes) {
      cols <- ITEMS[grepl(paste0("^(", prefixes, ")(_|$)"), ITEMS)]
      stopifnot(length(cols) == 4)  # every domain has 4 items (sleep: P5, P6_1-3)
      s <- rowSums(items[cols])
      50 + 10 * as.numeric(scale(s))
    }
    TS <- data.frame(
      T_PF  = tscore("P1"), T_ANX = tscore("P2"), T_DEP = tscore("P3"),
      T_FAT = tscore("P4"), T_SLP = tscore("P5|P6"), T_SR  = tscore("P7"),
      T_PI  = tscore("P8")
    )
    TS <- round(TS, 1)

    # EQ-5D-5L dimensions from the same latent factors
    eq_lat <- cbind(
      E1 = 0.8 * phys + 0.2 * ment,
      E2 = 0.7 * phys + 0.2 * ment,
      E3 = 0.6 * phys + 0.4 * ment,
      E4 = 0.7 * phys + 0.3 * ment,
      E5 = 0.2 * phys + 0.8 * ment
    )
    eq_cuts <- c(0.25, 0.95, 1.65, 2.35)
    E <- as.data.frame(sapply(DIMS, function(v) {
      to_level(eq_lat[, v] + rnorm(n, 0, 0.35), eq_cuts)
    }))

    # Diagnoses (0/1); prevalence rises with latent severity
    prev <- c(0.12, 0.10, 0.08, 0.15, 0.06, 0.07, 0.05, 0.08, 0.04, 0.06,
              0.10, 0.09, 0.05, 0.04, 0.20)
    names(prev) <- .SYN_DIAG
    D <- as.data.frame(sapply(.SYN_DIAG, function(v) {
      rbinom(n, 1, plogis(qlogis(prev[[v]]) + 0.4 * (phys + ment)))
    }))

    # Survey weight: mean 1, coefficient of variation about 0.3
    shape <- 1 / 0.3^2
    WEIGHT <- rgamma(n, shape, rate = shape)
    WEIGHT <- WEIGHT / mean(WEIGHT)
  })

  # U_own: utility on the value set of the respondent's own country
  U_own <- numeric(n)
  for (cc in unique(country)) {
    idx <- country == cc
    U_own[idx] <- .score(as.matrix(E[idx, DIMS]), p2e_value_set(cc))
  }

  # util_band: FH = full health (11111); Q1-Q4 = utility quartiles among the
  # remaining rows within each country (as in R/functions/splits.R)
  full_health <- rowSums(E == 1) == 5
  util_band <- rep("FH", n)
  for (cc in countries) {
    idx <- which(country == cc & !full_health)
    if (length(idx) == 0) next
    q <- ceiling(4 * rank(U_own[idx], ties.method = "first") / length(idx))
    util_band[idx] <- paste0("Q", q)
  }

  ID <- sprintf("%s%04d", country, seq_len(n))
  out <- data.frame(ID = ID, country = country, AGE = AGE, female = female,
                    WEIGHT = WEIGHT, util_band = util_band, items, P9 = P9,
                    TS, D, E, U_own = U_own, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out[.SYN_COLS]
}

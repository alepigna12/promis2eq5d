# Ported from iHEA_2027 R/functions/models/R06_mvop_dim.R at freeze-ihea2027-abstract-v3 (a79e102); changes: constructor .method_R06() instead of register_mvop_response(); fit() adds $set and $variant; simulate() takes a key; mean_utility() takes a p2e_vs value set; profile_prob not ported

# R06: multivariate ordered probit with a different predictor set for each
# dimension (spec §5.2, R6): R05 with the conceptual columns of each
# dimension (CONCEPTUAL_SUBSETS, item or T-score version, plus the set's
# covariates), as R02's conceptual variant. One variant, conceptual, in
# every set. Estimation, predictions and the mean: mvop.R. The fit holds no
# formula or environment (clm's coefficients are extracted as numbers), so
# saving a bundle does not serialise the training data.
.method_R06 <- function() {
  list(
    id = "R06_mvop_dim", family = "response",
    sets = METHOD_INFO$R06_mvop_dim$sets,
    variants = list(list(form = "conceptual")),
    aliases = c("R06", "mvop", "mvop_dim"),
    needs = c("ordinal", "pbivnorm", "Matrix"),
    predict_needs = character(0),
    fit = function(train, set, variant, oof = NULL) {
      f <- mvop_fit(train, list(set = set, variant = variant), link = "probit", label = "R06")
      f$set <- set
      f$variant <- variant
      f
    },
    simulate = function(fit, newdata, M, key) mvop_simulate_profiles(fit, newdata, M, key),
    mean_utility = function(fit, newdata, vs) mvop_mean_utility(fit, newdata, vs),
    marginals = function(fit, newdata) ordinal_marginals(fit, newdata)
  )
}

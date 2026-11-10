# Pipeline-equivalence fixtures (Task 12). Fits the four ported methods with the
# FROZEN research pipeline (iHEA_2027, tag freeze-ihea2027-abstract-v3) and stores
# the numbers the package must reproduce. The fixtures hold plain numbers only:
# no pipeline fit objects, no lookup tables, no item or level wording.
#
# Run from the PROMIS_to_EQ5D directory, on the normal R library (not iHEA's renv):
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript data-raw/pipeline_fixtures.R <iHEA repo> tests/testthat/fixtures
#
# The pipeline code is extracted with `git archive` into a tempdir (read-only on the
# iHEA repo) and sourced into a dedicated environment. Files sourced, in order:
# analysis_data_dictionary.R, the other R/functions/*.R alphabetically
# (beta, cli, design, evaluate, klapproth2020, load_split, metrics, ml_configs,
# model_interface, mvop, promis, report, results, selection, splits, synthetic,
# task_grid, tasks, twin_failures, twopart, value_sets), then the model files
# R/functions/models/*.R (alphabetical; each calls register_method()).
# Nothing extracted is copied into the package.
#
# Selection fixture (Plan 2, Task 4): pipeline_selection.rds holds the pipeline's
# ml_configs("xgb"), the row .select_cell() selects from three made-up metric
# tables, crps_draws() on a fixed 5 x 200 draw matrix, and the R09 profile
# probabilities of the first 5 newdata rows of the pipeline_R09.rds fit (refitted
# from that fixture's train, set, variant and seed). The pipeline's CRPS needs the
# scoringRules package on the normal library (1.1.3, as in iHEA's renv.lock).
# A third argument `selection` builds this fixture only and leaves the four
# method fixtures untouched:
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript data-raw/pipeline_fixtures.R <iHEA repo> tests/testthat/fixtures selection

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 2:3 || (length(args) == 3 && args[3] != "selection")) {
  stop("usage: Rscript pipeline_fixtures.R <iHEA repo> <output dir> [selection]")
}
only_selection <- length(args) == 3
repo <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
outdir <- normalizePath(args[2], winslash = "/", mustWork = TRUE)
tag <- "freeze-ihea2027-abstract-v3"

tmp <- tempfile("ihea_fn_")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
tarfile <- tempfile(fileext = ".tar")
on.exit(unlink(tarfile), add = TRUE)
rc <- system2("git", c("-C", shQuote(repo), "archive", "--format=tar", "-o", shQuote(tarfile), tag, "R/functions"))
if (rc != 0) stop("git archive failed")
utils::untar(tarfile, exdir = tmp)

if (capabilities("NLS")) Sys.setLanguage("en")
env <- new.env()                                   # pipeline objects live here, not in globalenv
fdir <- file.path(tmp, "R", "functions")
dict <- file.path(fdir, "analysis_data_dictionary.R")
sys.source(dict, envir = env)
fl <- sort(list.files(fdir, pattern = "[.]R$", full.names = TRUE))
for (f in fl[normalizePath(fl) != normalizePath(dict)]) sys.source(f, envir = env)
for (f in sort(list.files(file.path(fdir, "models"), pattern = "[.]R$", full.names = TRUE))) sys.source(f, envir = env)
env$ML_THREADS <- 1L                               # single-threaded xgboost

train <- env$make_synthetic_data(400, seed = 1, countries = c("DE", "UK"))
train$fold <- env$make_cv_folds(train, v = 3, seed = 1)$fold
train$y <- env$score_profiles(train[env$DIMS], "UK")        # direct methods' outcome: E1-E5 rescored on UK
newdata <- env$make_synthetic_data(10, seed = 2)            # 10 rows, default countries DE/UK
vt <- "UK"; fit_seed <- 1L; seed <- 11L; M <- 20L

mkspec <- function(id, set, variant, arm = "none", train_vt = NA_character_)
  list(method = id, set = set, arm = arm, train_vt = train_vt, variant = variant, weighted = FALSE, cw_weights = NA)
fit_pipeline <- function(id, spec) {
  set.seed(fit_seed)
  f <- env$METHODS[[id]]$fit(train, spec)
  f$spec <- spec
  class(f) <- c(id, "mapfit")
  f
}
num <- function(x) lapply(x, function(v) if (is.matrix(v)) unclass(v) else as.numeric(v))
save_fx <- function(id, set, variant, core, mean, u) {
  fx <- list(id = id, set = set, variant = variant, vt = vt, fit_seed = fit_seed, seed = seed, M = M,
             train = train, newdata = newdata, core = core, mean = mean, draws = list(u = u))
  saveRDS(fx, file.path(outdir, paste0("pipeline_", substr(id, 1, 3), ".rds")), version = 3)
}

# The four method fixtures (skipped with `selection`; the block is not indented,
# so the Plan 1 code below reads unchanged)
if (!only_selection) {
## R06: mean = predict_model(mean) under set.seed(12); draws = predict_model(draws, M) under set.seed(11)
v <- list(form = "conceptual")
f <- fit_pipeline("R06_mvop_dim", mkspec("R06_mvop_dim", "S1", v))
core <- list(alpha = lapply(f$dims, function(d) as.numeric(d$alpha)), beta = lapply(f$dims, function(d) as.numeric(d$beta)),
             coef = unclass(f$coef), R = unclass(f$R))
set.seed(12); mu <- as.numeric(env$predict_model(f, newdata, vt, "mean"))
set.seed(seed); u <- unname(env$predict_model(f, newdata, vt, "draws", M = M)$u)
save_fx("R06_mvop_dim", "S1", v, core, mu, u)

## R09: registered functions directly (the stand-in variant is not in ml_configs("xgb"))
v <- list(eta = 0.3, max_depth = 2, subsample = 1, colsample_bytree = 1, min_child_weight = 1, nrounds = 10)
f <- fit_pipeline("R09_chain", mkspec("R09_chain", "S3", v))
set.seed(seed); u <- unname(env$score_profiles(env$METHODS$R09_chain$simulate_profiles(f, newdata, M), vt))
set.seed(12); mu <- as.numeric(rowMeans(env$score_profiles(env$METHODS$R09_chain$simulate_profiles(f, newdata, 2000L), vt)))
save_fx("R09_chain", "S3", v, NULL, mu, u)

## D03 and D12: predict_model() (mean = clamped point; draws = clamped as_published draws)
v <- list(form = "linear")
f <- fit_pipeline("D03_clad", mkspec("D03_clad", "S3", v, "rescored", "UK"))
core <- list(beta = as.numeric(f$beta), oof_mu = as.numeric(f$oof_mu), oof_resid = as.numeric(f$oof_resid))
mu <- as.numeric(env$predict_model(f, newdata, vt, "mean"))
set.seed(seed); u <- unname(env$predict_model(f, newdata, vt, "draws", M = M)$u)
save_fx("D03_clad", "S3", v, core, mu, u)

v <- list(form = "klapproth")
f <- fit_pipeline("D12_klapproth", mkspec("D12_klapproth", "S3K", v, "rescored", "UK"))
core <- list(terms_kept = f$terms_kept, beta = as.numeric(f$beta), beta_names = names(f$beta), sigma = as.numeric(f$sigma))
mu <- as.numeric(env$predict_model(f, newdata, vt, "mean"))
set.seed(seed); u <- unname(env$predict_model(f, newdata, vt, "draws", M = M)$u)
save_fx("D12_klapproth", "S3K", v, core, mu, u)
cat("fixtures written to", outdir, "\n")
}

## Selection fixture (see the header)
# (a) the R09 grid
ml <- env$ml_configs("xgb")

# (b) three made-up metric tables in the package's form (variant_id, ok, mae,
# tie_metric), selected by the pipeline's .select_cell() on the ok rows (the
# pipeline's candidates are the usable units). tie_metric is cv_neglogp for a
# response cell and cv_crps for a direct cell.
tables <- list(
  # no tie: the lowest-MAE row (id 2) failed; id 3 wins over id 4 (0.005 apart)
  data.frame(variant_id = 1:4, ok = c(TRUE, FALSE, TRUE, TRUE), mae = c(0.120, 0.100, 0.110, 0.115),
             tie_metric = c(2.0, 1.0, 3.0, 1.5)),
  # two-way tie at the 0.0005 boundary (ids 1, 2); id 3 (0.0006 above) is not tied
  data.frame(variant_id = 1:5, ok = TRUE, mae = c(0.1000, 0.1005, 0.1006, 0.1200, 0.1300),
             tie_metric = c(0.30, 0.25, 0.10, 0.05, 0.01)),
  # three-way tie (ids 7, 4, 9): 7 and 4 share tie_metric and MAE, so variant_id decides
  data.frame(variant_id = c(7L, 4L, 9L, 2L, 5L), ok = TRUE, mae = c(0.2002, 0.2002, 0.2000, 0.2100, 0.2500),
             tie_metric = c(1.0, 1.0, 1.2, 0.1, 0.5))
)
families <- c("response", "direct", "response")
pipeline_select <- function(tab, family) {
  d <- tab[tab$ok, , drop = FALSE]
  resp <- family == "response"
  cells <- data.frame(scenario = "holdout_UK", method = "m", set = "S3", arm = if (resp) "none" else "rescored",
                      train_vt = if (resp) NA_character_ else "UK", cw_weights = NA_character_, family = family,
                      variant_id = d$variant_id, cv_mae = d$mae,
                      cv_neglogp = if (resp) d$tie_metric else NA_real_,
                      cv_crps = if (resp) NA_real_ else d$tie_metric, stringsAsFactors = FALSE)
  env$.select_cell(cells)
}
picked <- Map(pipeline_select, tables, families)

# (c) CRPS on a fixed 5 x 200 draw matrix, with mass at full health and a draw equal to y
set.seed(20261110)
crps_draws <- matrix(round(runif(1000, -0.2, 1), 3), 5, 200)
crps_draws[1:2, 1:60] <- 1
crps_y <- c(1, 0.9, 0.5, -0.1, crps_draws[5, 17])
crps <- as.numeric(env$crps_draws(crps_y, crps_draws))

# (d) R09 profile probabilities of the pipeline_R09.rds fit on its first 5 newdata rows
r09 <- readRDS(file.path(outdir, "pipeline_R09.rds"))
if (!identical(r09$train, train)) stop("pipeline_R09.rds was built from another train; rebuild it first")
spec <- mkspec("R09_chain", r09$set, r09$variant)
set.seed(r09$fit_seed)
f <- env$METHODS$R09_chain$fit(r09$train, spec)
f$spec <- spec
class(f) <- c("R09_chain", "mapfit")
nd5 <- r09$newdata[1:5, , drop = FALSE]
obs <- matrix(as.integer(as.matrix(nd5[env$DIMS])), ncol = 5)
r09_prob <- as.numeric(env$METHODS$R09_chain$profile_prob(f, nd5, obs))
r09_logp <- as.numeric(env$log_profile_prob(f, nd5))

fx <- list(ml_configs = ml, tables = tables, families = families,
           selected = lapply(picked, function(s) s$variant_id),
           n_tied = vapply(picked, function(s) as.integer(s$n_tied), 1L),
           n_candidates = vapply(picked, function(s) as.integer(s$n_candidates), 1L),
           crps_y = crps_y, crps_draws = crps_draws, crps = crps,
           r09_rows = 1:5, r09_prob = r09_prob, r09_logp = r09_logp,
           scoringRules = as.character(utils::packageVersion("scoringRules")))
saveRDS(fx, file.path(outdir, "pipeline_selection.rds"), version = 3)
cat("selection fixture written to", outdir, "\n")

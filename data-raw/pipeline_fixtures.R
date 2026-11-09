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

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: Rscript pipeline_fixtures.R <iHEA repo> <output dir>")
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

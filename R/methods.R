# Method registry. Each method file defines a constructor function only
# (.method_R06(), .method_R09(), .method_D03(), .method_D12()); .onLoad() calls
# .register_all_methods(), which registers every namespace function named
# `.method_<id letter><digits>` (the constructors; helpers such as .method_variants() do not match).
#
# A method is a list:
#   id, family            "R06_mvop_dim", "response" or "direct"
#   sets                  supported predictor sets, preferred order first
#   variants              list of variant lists, or function(set) returning one; the first
#                         is the default; resolve with .method_variants(m, set)
#   aliases               alternative names, e.g. c("R06", "mvop")
#   needs                 Suggests packages needed to FIT
#   predict_needs         Suggests packages needed to PREDICT
#   uses_oof              optional logical (absent = FALSE). When TRUE,
#                         .fit_bundle() computes out-of-fold pairs
#                         oof = data.frame(oof_mu, oof_resid) from the fold
#                         fits, passes them to fit() of the point fit, and
#                         passes oof[idx, ] to bootstrap replicate b, with the
#                         indices idx that resample `train`. Fold fits get
#                         oof = FALSE. Direct methods only (D03).
#                         The oof argument of fit() means: a data.frame
#                         (oof_mu, oof_resid) = use these residual pairs;
#                         NULL = compute your own pairs (standalone fit, as in
#                         the pipeline); FALSE = no pairs needed, the fit is
#                         used for point predictions only (skip the nested
#                         out-of-fold step).
#   fit(train, set, variant, oof = NULL)    the fit list (+ $set, $variant);
#                         may set $converged = FALSE. Direct methods find the
#                         outcome in train$y; train has a `fold` column.
#   response only:  simulate(fit, newdata, M, key), mean_utility(fit, newdata, vs),
#                   marginals(fit, newdata)
#   direct only:    point(fit, newdata), draw(fit, newdata, M, key, draw_type, lo)

.p2e_methods <- new.env(parent = emptyenv())

.register_p2e_method <- function(m) {
  need <- c("id", "family", "sets", "variants", "fit")
  miss <- setdiff(need, names(m))
  if (length(miss)) stop("method is missing: ", paste(miss, collapse = ", "))
  if (!m$family %in% c("response", "direct")) stop("method family must be 'response' or 'direct'")
  if (!is.list(m$variants) && !is.function(m$variants)) stop("method variants must be a list or a function of the set")
  assign(m$id, m, envir = .p2e_methods)
  invisible(m)
}

.register_all_methods <- function() {
  ns <- topenv(environment())
  for (nm in ls(ns, pattern = "^\\.method_[A-Z][0-9]", all.names = TRUE)) {
    f <- get(nm, envir = ns)
    if (is.function(f)) .register_p2e_method(f())
  }
  invisible(NULL)
}

# A registered method by id or alias
.p2e_method <- function(id_or_alias) {
  ms <- mget(ls(.p2e_methods), envir = .p2e_methods)
  if (is.character(id_or_alias) && length(id_or_alias) == 1L && !is.na(id_or_alias)) {
    if (!is.null(ms[[id_or_alias]])) return(ms[[id_or_alias]])
    for (m in ms) if (id_or_alias %in% m$aliases) return(m)
  }
  stop("unknown method '", paste(format(id_or_alias), collapse = " "), "'; available methods: ",
       paste(names(ms), collapse = ", "), call. = FALSE)
}

# The variants of method m for predictor set `set`: m$variants(set) if it is a
# function, else m$variants. Must be a non-empty list.
.method_variants <- function(m, set) {
  v <- if (is.function(m$variants)) m$variants(set) else m$variants
  if (!is.list(v) || !length(v)) {
    stop("method ", if (is.null(m$id)) "" else paste0("'", m$id, "' "),
         "variants must be a non-empty list (or a function of the set returning one)", call. = FALSE)
  }
  v
}

# Tests for value-set objects, the 3,125-state index and scoring.
# Ported in part from iHEA_2027 tests/testthat/test-value-sets.R (state_index,
# score_profiles, vt_min, vt_gap, clamp_utility, vt_additive), rewritten for p2e_vs.

test_that("lookup by name, code and alias gives the eq5d value set", {
  uk <- p2e_value_set("UK")
  expect_equal(length(uk$u), 3125); expect_equal(uk$u[1], 1)
  expect_equal(uk$u, eq5d::eq5d(.all_states_df(), version = "5L", type = "VT", country = "UK"))
  expect_identical(p2e_value_set("gbr")$u, uk$u)
  expect_identical(p2e_value_set("United States")$name, "USA")
  expect_message(p2e_value_set("Japan"), "Japan_DCE")
  expect_identical(p2e_value_set("Japan", version = "DCE")$name, "Japan_DCE")
  expect_error(p2e_value_set("Narnia"), "no EQ-5D-5L value set")
  expect_identical(p2e_value_set("usa")$u, p2e_value_set("United States")$u)
  expect_identical(p2e_value_set("GB")$u, uk$u); expect_identical(p2e_value_set("England")$name, "England")
  expect_identical(p2e_value_set("Japan", version = "DCE")$version, "DCE")
  expect_error(p2e_value_set("Germny"), "Germany")                 # close matches are listed
})
test_that("codes, defaults for multi-version countries and version checks", {
  expect_identical(p2e_value_set("UK")$code, "UK")
  expect_identical(p2e_value_set("France")$code, "FR")
  expect_identical(p2e_value_set("deu")$name, "Germany")
  expect_identical(p2e_value_set("Japan", version = "DCE")$code, "JP_DCE")
  expect_identical(suppressMessages(p2e_value_set("japan"))$name, "Japan_cTTO")
  expect_identical(suppressMessages(p2e_value_set("peru"))$name, "Peru_cTTO")
  expect_identical(suppressMessages(p2e_value_set("sweden"))$name, "Sweden_2022")
  expect_identical(p2e_value_set("Sweden", version = "2020")$name, "Sweden_2020")
  expect_identical(p2e_value_set("UK")$version, NA_character_)
  expect_identical(p2e_value_set("UK")$source, "eq5d")
  expect_error(p2e_value_set("UK", version = "DCE"), "version")
  expect_error(p2e_value_set("Japan", version = "XYZ"), "version")
  expect_s3_class(p2e_value_set("UK"), "p2e_vs")
  vs <- p2e_value_set("DE"); expect_identical(p2e_value_set(vs), vs)
  expect_error(p2e_value_set(1), "country")
})
test_that("alias table covers every eq5d 5L country", {
  vt <- eq5d::valuesets(type = "VT", version = "5L")$Country
  expect_true(all(vt %in% .vs_env$aliases$eq5d_name))
  expect_identical(.vs_env$aliases$alias, tolower(.vs_env$aliases$alias))
  expect_false(anyDuplicated(.vs_env$aliases$alias) > 0)
})
# Verified at the freeze tag: the pipeline builds its table with eq5d::eq5d() for
# "France", "Germany" and "UK", and its tests assert these constants.
test_that("FR, DE, UK minima and gaps match the pipeline table", {
  expect_equal(.vs_min(p2e_value_set("FR")), -0.525, tolerance = 1e-9)
  expect_equal(.vs_gap(p2e_value_set("DE")), 0.026, tolerance = 1e-9)
  expect_equal(.vs_gap(p2e_value_set("UK")), 0.032, tolerance = 1e-9)
  expect_equal(c(.vs_min(p2e_value_set("DE")), .vs_min(p2e_value_set("UK"))), c(-0.661, -0.567), tolerance = 1e-9)
  expect_equal(.vs_gap(p2e_value_set("FR")), 0.020, tolerance = 1e-9)
})
test_that("user tables must cover all states once", {
  t <- cbind(.all_states_df(), utility = p2e_value_set("DE")$u)
  expect_identical(p2e_value_set(t)$u, p2e_value_set("DE")$u)
  expect_identical(p2e_value_set(t)$code, "user")
  expect_identical(p2e_value_set(t)$source, "user")
  expect_error(p2e_value_set(t[-1, ]), "3,125")
  expect_error(p2e_value_set(rbind(t, t[1, ])), "3,125")
  expect_error(p2e_value_set(t[, -6]), "utility")
  tt <- t; tt$utility[5] <- NA; expect_error(p2e_value_set(tt), "utility")
  shuffled <- t[sample(3125), ]                                     # row order is free
  expect_identical(p2e_value_set(shuffled)$u, p2e_value_set("DE")$u)
  t2 <- t; t2$utility <- rev(t$utility)                            # same code "user", other utilities: never cached
  expect_false(isTRUE(all.equal(.vs_additive(p2e_value_set(t))$intercept, .vs_additive(p2e_value_set(t2))$intercept)))
})
test_that("additive decomposition reproduces the table", {
  for (cc in c("FR", "DE", "UK", "USA")) {
    vs <- p2e_value_set(cc); a <- .vs_additive(vs); S <- .all_states_df()
    add <- a$intercept - rowSums(sapply(1:5, function(j) a$decrement[j, S[[j]]]))
    expect_equal(add + a$remainder, vs$u, tolerance = 1e-12)
  }
  expect_identical(.vs_additive(p2e_value_set("DE")), .vs_additive(p2e_value_set("DE")))  # cached
  expect_lt(max(abs(.vs_additive(p2e_value_set("UK"))$remainder)), 1e-12)   # UK is additive (raw residual 2.7e-13, stored as 0)
  expect_equal(dim(.vs_additive(p2e_value_set("UK"))$decrement), c(5, 5))
  expect_true(all(.vs_additive(p2e_value_set("UK"))$decrement[, 1] == 0))
  expect_lte(max(abs(.vs_additive(p2e_value_set("FR"))$remainder)), 0.00056)  # FR: additive model rounded to 3 decimals
})
test_that("state order and scoring match eq5d", {
  S <- .all_states_df()
  expect_equal(dim(S), c(3125L, 5L)); expect_identical(names(S), c("MO", "SC", "UA", "PD", "AD"))
  st <- do.call(paste0, S)
  expect_identical(st[c(1, 2, 6, 3125)], c("11111", "11112", "11121", "55555"))
  expect_identical(.score(as.matrix(S), p2e_value_set("UK")), p2e_value_set("UK")$u)
  set.seed(3); p <- matrix(sample(1:5, 100, TRUE), 20)
  for (cc in c("FR", "DE", "UK")) expect_equal(.score(p, p2e_value_set(cc)),
    eq5d::eq5d(setNames(as.data.frame(p), c("MO", "SC", "UA", "PD", "AD")), version = "5L", type = "VT",
               country = c(FR = "France", DE = "Germany", UK = "UK")[[cc]]), tolerance = 1e-10)
})
test_that(".score handles value-set lists, draw arrays and bad input", {
  fr <- p2e_value_set("FR"); de <- p2e_value_set("DE"); uk <- p2e_value_set("UK")
  p <- rbind(c(1, 1, 1, 1, 2), c(2, 1, 1, 1, 1), c(5, 5, 5, 5, 5))
  vl <- list(fr, de, uk)
  expect_identical(.score(p, vl), c(fr$u[2], de$u[626], uk$u[3125]))
  expect_identical(.score(as.data.frame(p), vl), .score(p, vl))
  a <- array(0L, c(3, 2, 5)); a[, 1, ] <- p; a[, 2, ] <- p[c(3, 1, 2), ]
  u <- .score(a, vl)
  expect_equal(dim(u), c(3, 2))
  for (i in 1:3) for (j in 1:2) expect_identical(u[i, j], .score(rbind(a[i, j, ]), vl[[i]]))
  expect_error(.score(cbind(p[, 1:4], 6), uk), "1-5")
  expect_error(.score(rbind(c(1, 1, NA, 1, 1)), uk), "1-5")
  expect_error(.score(data.frame(a = "1", b = "1", c = "1", d = "1", e = "1"), uk), "1-5")
  expect_error(.score(p, list(uk, de)), "value set")
  expect_error(.score(p, "UK"), "value set")
})
test_that(".state_index gives the row of the state table for each profile shape", {
  S <- .all_states_df(); st <- do.call(paste0, S)
  expect_identical(.state_index(rbind(rep(1, 5))), 1L)
  expect_identical(.state_index(rbind(rep(5, 5))), 3125L)
  set.seed(5); p <- matrix(sample(1:5, 500, TRUE), 100)
  i <- .state_index(p)
  expect_identical(length(i), 100L)
  expect_identical(st[i], do.call(paste0, as.data.frame(p)))
  expect_identical(.state_index(as.data.frame(p)), i)
  a <- array(sample(1:5, 100 * 3 * 5, TRUE), c(100, 3, 5))
  ia <- .state_index(a)
  expect_equal(dim(ia), c(100, 3))
  for (j in 1:3) expect_identical(ia[, j], .state_index(a[, j, ]))
  expect_error(.state_index(cbind(p[, 1:4], 6)), "1-5")
  expect_error(.state_index(matrix(1, 3, 4)), "n x 5")
})
test_that("value-set minimum and gap", {
  uk <- p2e_value_set("UK"); fr <- p2e_value_set("FR"); de <- p2e_value_set("DE")
  expect_identical(.vs_min(uk), min(uk$u))
  expect_identical(.vs_min(de), de$u[3125])                          # 55555 is the worst state
  expect_equal(.vs_gap(fr), 1 - max(fr$u[fr$u < 1]))
})
test_that(".clamp keeps attributes and refuses missing values", {
  uk <- p2e_value_set("UK"); fr <- p2e_value_set("FR")
  expect_identical(.clamp(c(a = -5, b = 0.5), uk), c(a = .vs_min(uk), b = 0.5))
  m <- matrix(c(-5, 2, .5, .4), 2, dimnames = list(c("r1", "r2"), c("d1", "d2")))
  out <- .clamp(m, list(uk, fr))
  expect_identical(dimnames(out), dimnames(m))
  expect_identical(unname(out), matrix(c(.vs_min(uk), 1, .5, .4), 2))
  expect_error(.clamp(c(.5, NA), uk), "missing")
  expect_error(.clamp(c(.5, .5, .5), list(uk, fr)), "length")
})

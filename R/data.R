#' Training data of the shipped fits
#'
#' The respondents the shipped mapping fits were trained on, and the data used
#' to retrain a direct method on another value set (`retrain = TRUE`).
#' Items are in PROMIS scoring orientation: `P7_1` to `P7_4` are stored with
#' the opposite direction in the source file and are reversed (`6 - x`); every
#' other item is used as stored. There are no T-scores and no respondent
#' identifier.
#'
#' @format A data frame with 4,512 rows and 37 columns, in this order:
#' \itemize{
#'   \item `P1_1` to `P8_4`: the 28 PROMIS-29 items, integer 1 to 5, named as in the Dataverse file (`P1_1` ... `P8_4`); [p2e_columns()] maps them (`internal` to `column`) for renaming before [predict_eq5d()].
#'   \item `P9`: pain intensity, integer 0 to 10.
#'   \item `AGE`: age in years, numeric.
#'   \item `female`: integer, 1 for female and 0 otherwise.
#'   \item `country`: character, `"FR"`, `"DE"` or `"UK"`.
#'   \item `E1` to `E5`: the five EQ-5D-5L dimension responses, integer 1 to 5.
#' }
#' @source Dataverse, doi:10.7910/DVN/FNYYQR (CC0). Built by `data-raw/p2e_training.R`.
"p2e_training"

#' Synthetic example trial
#'
#' A small synthetic two-arm trial with three time points, for the vignettes
#' and examples. The respondents are resampled from the UK rows of
#' [p2e_training] and followed up with a made-up treatment effect; it is not a
#' real trial, and the treatment effect is made up. Items are in PROMIS
#' scoring orientation and carry the user-facing item column names, ready for
#' [predict_eq5d()]. There are no EQ-5D columns, no T-scores and no country
#' column.
#'
#' @format A data frame with 900 rows (300 persons at 3 times) and 33 columns, in this order:
#' \itemize{
#'   \item `id`: person identifier, integer 1 to 300.
#'   \item `arm`: character, `"control"` or `"treatment"` (150 persons each).
#'   \item `month`: integer 0, 6 or 12.
#'   \item 28 item columns, named by item ID (see [p2e_columns()]): integer 1 to 5.
#'   \item `Global07`: pain intensity, integer 0 to 10.
#'   \item `age`: age in years at that time, numeric.
#'   \item `female`: integer, 1 for female and 0 otherwise.
#' }
#' @source Synthetic; built by `data-raw/p2e_example_trial.R` from [p2e_training] (Dataverse doi:10.7910/DVN/FNYYQR, CC0).
"p2e_example_trial"

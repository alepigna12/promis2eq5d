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
#' \describe{
#'   \item{`P1_1` to `P8_4`}{The 28 PROMIS-29 items, integer 1 to 5, named by item ID.}
#'   \item{`P9`}{Pain intensity, integer 0 to 10.}
#'   \item{`AGE`}{Age in years, numeric.}
#'   \item{`female`}{Integer, 1 for female and 0 otherwise.}
#'   \item{`country`}{Character, `"FR"`, `"DE"` or `"UK"`.}
#'   \item{`E1` to `E5`}{The five EQ-5D-5L dimension responses, integer 1 to 5.}
#' }
#' @source Dataverse, doi:10.7910/DVN/FNYYQR (CC0). Built by `data-raw/p2e_training.R`.
"p2e_training"

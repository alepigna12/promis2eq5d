# promis2eq5d

promis2eq5d maps PROMIS-29 profiles to EQ-5D-5L utilities, with uncertainty. For each person it predicts either the five EQ-5D-5L dimension responses, which are then scored under the value set of any country in the 'eq5d' package, or the utility directly on the scale of a shipped value set. It can return a mean prediction with a conformal interval, individual draws, or draws of the cohort mean for use in cost-effectiveness and QALY analyses. Prediction models are chosen automatically from the columns in your data (items or T-scores, pain intensity, age and sex).

## Installation

The package is not on CRAN. Once the repository is public, install it from GitHub (placeholder until then):

```r
# install.packages("remotes")
remotes::install_github("alepigna12/promis2eq5d")
```

The fitting packages (ordinal, pbivnorm, Matrix, quantreg, glmnet, xgboost) are in Suggests and are needed only for the methods that use them.

## Example

```r
library(promis2eq5d)
predict_eq5d(my_data, country = "UK")
predict_eq5d(my_data, country = "Japan", mode = "individual", draws = 200, seed = 1)
```

`p2e_columns()` lists the columns the function reads.

## Three modes

* `mode = "mean"`: one expected utility per row, optionally with a conformal interval (`interval = TRUE`) and the most likely dimension levels (`dimensions = TRUE`).
* `mode = "individual"`: `draws` draws per row, for individual-level uncertainty.
* `mode = "cohort"`: for each draw, the mean utility of each group (`group = ...`), with draws aligned across groups and time points for QALY calculations.

## Large downloads (R09)

The R09_chain response method has large fits that are not part of the package. Download them once with `p2e_download("R09_chain")`. The function asks before downloading in an interactive session, and in a non-interactive session needs `options(promis2eq5d.allow_download = TRUE)`. Files are checked against a published checksum and installed in a user cache folder.

## T-score tables

PROMIS T-score conversion tables are not included in this package. If your data have raw domain scores, supply a conversion table obtained from the HealthMeasures Scoring Service (<https://www.healthmeasures.net/score-and-interpret/calculate-scores>) through the `tscore_table` argument (see `?p2e_tscores`), or supply T-score columns directly. Item-level predictor sets need no table.

## Trademarks

EQ-5D™ is a trade mark of the EuroQol Research Foundation.

PROMIS® is a registered trademark of the U.S. Department of Health and Human Services. This package is not affiliated with or endorsed by HHS, the PROMIS Health Organization, HealthMeasures or the EuroQol Research Foundation.

## Licence

MIT.

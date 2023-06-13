---
title: "NEWS.md"
author: "Simon Dedman"
date: "2023-05-25"
output: html_document
---
# v2023.06.13
* CRAN release

# v2023.05.23
* added stop condition when all resvar are zero

# v2023.05.22
* dismo1.3.14 pushed to CRAN, updated dependency

# v2023.05.18
* readme - gbm.step.sd added, changed devtools to remotes
* gbm.auto - dismo 1.3.10 remotes install note added, dev.print fix so it saves residual deviance lineplots for each model combo, ablines added at y=0 for pdp lineplots

# v2023.03.16
* added categorical_pdp_plotter to extras
* added gbm.plots and lmplot
* random var added
* gbm.factorplot added
* add gbm.step.sd and model evaluation metrics like RMSE etc

# v1.5.5
* gbm.basemap added returnsf param, allows returning sf object

# v1.5.2
* gbm.basemap added sf::sf_use_s2(FALSE)

# v1.5.1
* More reporting info - Self_CV_Statistics.csv created with all outputs from all models for those 2 gbm.step objects
* gbm.loop automatically checks for the presence of Report.csv's in numbered folders and doesn't run for that folder if so

# v1.5.0
* gbm.loop fixed for running gbm.auto & post-run-results-gathering loops separately.
* gbm.loop added all params for internal call to gbm.auto.

# v1.4.1, 2021-02-23
* match.arg removed for ZI and res in basemap and elsewhere. Function can't have mixed inputs i.e. character and logical/numeric (#​65)

# v1.4, 2021-01-14: CRAN release
* Reduced large datasets to 75000 rows for CRAN size constraints.
* gbm.utils split into roc calibration & gbm.predict.grids. g.p.g. merged into gbm.auto and removed.


<!-- If an item is related to an issue in GitHub, include the issue number in parentheses, e.g. (#​10).
If an item is related to a pull request, include the pull request number and the author, e.g. (#​101, @hadley).
Doing this makes it easy to navigate to the relevant issues on GitHub.-->

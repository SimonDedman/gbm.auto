---
title: "NEWS.md"
author: "Simon Dedman"
date: "2021-01-14"
output: html_document
---

# v1.5.0
* gbm.loop fixed for running gbm.auto & post-run-results-gatehring loops separately.
* gbm.loop added all params for internal call to gbm.auto.

# v1.4.1, 2021-02-23
* match.arg removed for ZI and res in basemap and elsewhere. Function can't have mixed inputs i.e. character and logical/numeric (#​65)

# v1.4, 2021-01-14: CRAN release
* Reduced large datasets to 75000 rows for CRAN size constraints.
* gbm.utils split into roc calibration & gbm.predict.grids. g.p.g. merged into gbm.auto and removed.


<!-- If an item is related to an issue in GitHub, include the issue number in parentheses, e.g. (#​10).
If an item is related to a pull request, include the pull request number and the author, e.g. (#​101, @hadley).
Doing this makes it easy to navigate to the relevant issues on GitHub.-->

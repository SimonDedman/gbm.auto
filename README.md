
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gbm.auto

<!-- badges: start -->

[![R-CMD-check](https://github.com/SimonDedman/gbm.auto/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/SimonDedman/gbm.auto/actions/workflows/R-CMD-check.yaml)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/gbm.auto)](https://cran.r-project.org/package=gbm.auto)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/gbm.auto)](https://cran.r-project.org/package=gbm.auto)
<!-- badges: end -->
<!-- badgeplacer(location = ".", status = "active", githubaccount = SimonDedman, githubrepo = gbm.auto, branch = master, name = "README.Rmd") -->

Automatically runs numerous processes from R packages ‘gbm’ and ‘dismo’
and script ‘gbm.utils.R’ which contains Elith et al.’s functions: roc,
calibration, and gbm.predict.grids, as well as running my packages
gbm.bfcheck, gbm.basemap, gbm.map, gbm.rsb, gbm.cons, gbm.valuemap, and
gbm.loop.

Also see each script’s Details section in the manual pages, as these
frequently contain tips or common bugfixes.

I strongly recommend that you download papers 1 to 5 (or just the
doctoral thesis) on <https://www.simondedman.com>, with emphasis on P4
(the guide) and P1 (statistical background). Elith et al 2008
(<https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.1365-2656.2008.01390.x>)
is also strongly recommended. Also it’s imperative you read the R help
files for each function before you use them. In RStudio: Packages tab,
scroll to gbm.auto, click its name, the click the function to see its
man (manual) page. Read the whole thing. Function man pages can also be
accessed from the console by typing

``` r
?function
```

Just because you CAN try every conceivable combination of tc, lr, bf,
all, at once doesn’t mean you should. Try a range of lr in shrinking
orders of magnitude from 0.1 to 0.000001, find the best, THEN try tc
c(2, n.expvars), find the best THEN bf c(0.5, 0.75, 0.9) and then in
between if either outperform 0.5.

------------------------------------------------------------------------

### gbm.auto

Automated Boosted Regression Tree modelling and mapping suite

Automates delta log normal boosted regression trees abundance
prediction. Loops through all permutations of parameters provided
(learning rate, tree complexity, bag fraction), chooses the best, then
simplifies it. Generates line, dot and bar plots, and outputs these and
the predictions and a report of all variables used, statistics for
tests, variable interactions, predictors used and dropped, etc. If
selected, generates predicted abundance maps, and Unrepresentativeness
surfaces.

------------------------------------------------------------------------

### gbm.bfcheck

Calculates minimum Bag Fraction size for gbm.auto

Provides minimum bag fractions for gbm.auto, preventing failure due to
bf & samples rows limit.

------------------------------------------------------------------------

### gbm.basemap

Creates Basemaps for Gbm.auto mapping from your data range

Downloads unzips crops & saves NOAAs global coastline shapefiles to
user-set box. Use for ‘shape’ in gbm.map. If downloading in RStudio
uncheck “Use secure download method for HTTP” in Tools \> Global Options
\> Packages.

------------------------------------------------------------------------

### gbm.map

Maps of predicted abundance from Boosted Regression Tree modelling

Generates maps from the outputs of gbm.step then gbm.predict.grids,
handled automatically within gbm.auto but can be run alone, and
generates representativeness surfaces from the output of gbm.rsb.

------------------------------------------------------------------------

### gbm.rsb

Representativeness Surface Builder

Loops through explanatory variables comparing their histogram in
‘samples’ to their histogram in ‘grids’ to see how well the explanatory
variable range in samples represents the range being predicted to in
grids. Assigns a representativeness score per variable per site in
grids, and takes the average score per site if there’s more than 1
expvar. Saves this to a CSV; it’s plotted by gbm.map if called in
gbm.auto. This shows you which areas have the most and least
representative coverage by samples, therefore where you can have the
most/least confidence in the predictions from gbm.predict.grids. Can be
called directly, and choosing a subset of expvars allows one to see
their individual / collective representativeness.

------------------------------------------------------------------------

### gbm.cons

Conservation Area Mapping

Runs gbm.auto for multiple subsets of the same overall dataset and
scales the combined results, leading to maps which highlight areas of
high conservation importance for multiple species in the same study area
e.g. using juvenile and adult female subsets to locate candidate nursery
grounds and spawning areas respectively.

------------------------------------------------------------------------

### gbm.valuemap

Decision Support Tool that generates (Marine) Protected Area options
using species predicted abundance maps

Scales response variable data, maps a user-defined explanatory variable
to be avoided, e.g. fishing effort, combines them into a map showing
areas to preferentially close. Bpa, the precautionary biomass required
to protect the spawning stock, is used to calculate MPA size. MPA is
then grown to add subsequent species starting from the most
conservationally at-risk species, resulting in one MPA map per species,
and a multicolour MPA map of all. All maps list the percentage of the
avoid-variables total that is overlapped by the MPA in the map legend.

------------------------------------------------------------------------

### gbm.loop

Calculate Coefficient Of Variation surfaces for gbm.auto predictions

Processes a user-specified number of loops through the same gbm.auto
parameter combinations and calculates the Coefficient Of Variation in
the predicted abundance scores for each site aka cell. This can be
mapped to spatially demonstrate the output variance range.

------------------------------------------------------------------------

### gbm.factorplot

ggplot-based update to PDP for factorial/categorical/character
variables, allows changing order of categorical variables, and changing
angle of x-axis labels to avoid them being cut off.

------------------------------------------------------------------------

### lmplot

Linear plot of two variables.

------------------------------------------------------------------------

### gbm.lmplots

Loops through lmplots for all expvars (x) against the same resvar (y).

------------------------------------------------------------------------

### roc & calibration

Internal functions authored by Elith & Leathwick, used by gbm.auto.R

------------------------------------------------------------------------

### gbm.step.sd

Local copy of dismo’s gbm.step, with added functions to generate model
evaluation metrics such as root mean squared error and amount of
deviance explained relative to null.

------------------------------------------------------------------------

## Installation

You can install the released version of gbm.auto from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("gbm.auto")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
remotes::install_github("SimonDedman/gbm.auto")
```

------------------------------------------------------------------------

## Example

(See each function’s help file for specific examples, and the documents
listed above)

------------------------------------------------------------------------

## ToDo List

See GitHub issues section
<https://github.com/SimonDedman/gbm.auto/issues> Feel free to contribute
to this!

<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->
<!-- ```{r cars} -->
<!-- summary(cars) -->
<!-- ``` -->
<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/master/examples>. -->
<!-- You can also embed plots, for example: -->
<!-- ```{r pressure, echo = FALSE} -->
<!-- plot(pressure) -->
<!-- ``` -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->

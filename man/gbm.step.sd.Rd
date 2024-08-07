% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gbm.step.sd.R
\name{gbm.step.sd}
\alias{gbm.step.sd}
\title{Function to assess optimal no of boosting trees using k-fold cross validation}
\usage{
gbm.step.sd(
  data,
  gbm.x,
  gbm.y,
  offset = NULL,
  fold.vector = NULL,
  tree.complexity = 1,
  learning.rate = 0.01,
  bag.fraction = 0.75,
  site.weights = rep(1, nrow(data)),
  var.monotone = rep(0, length(gbm.x)),
  n.folds = 10,
  prev.stratify = TRUE,
  family = "bernoulli",
  n.trees = 50,
  step.size = n.trees,
  max.trees = 10000,
  tolerance.method = "auto",
  tolerance = 0.001,
  plot.main = TRUE,
  plot.folds = FALSE,
  verbose = TRUE,
  silent = FALSE,
  keep.fold.models = FALSE,
  keep.fold.vector = FALSE,
  keep.fold.fit = FALSE,
  ...
)
}
\arguments{
\item{data}{The input dataframe.}

\item{gbm.x}{The predictors.}

\item{gbm.y}{The response.}

\item{offset}{Allows an offset to be specified.}

\item{fold.vector}{Allows a fold vector to be read in for CV with offsets,.}

\item{tree.complexity}{Sets the complexity of individual trees.}

\item{learning.rate}{Sets the weight applied to inidivudal trees.}

\item{bag.fraction}{Sets the proportion of observations used in selecting variables.}

\item{site.weights}{Allows varying weighting for sites.}

\item{var.monotone}{Restricts responses to individual predictors to monotone.}

\item{n.folds}{Number of folds.}

\item{prev.stratify}{Prevalence stratify the folds - only for p/a data.}

\item{family}{Family - bernoulli (=binomial), poisson, laplace or gaussian.}

\item{n.trees}{Number of initial trees to fit.}

\item{step.size}{Numbers of trees to add at each cycle.}

\item{max.trees}{Max number of trees to fit before stopping.}

\item{tolerance.method}{Method to use in deciding to stop - "fixed" or "auto".}

\item{tolerance}{Tolerance value to use - if method == fixed is absolute, if auto is multiplier * total mean deviance.}

\item{plot.main}{Plot hold-out deviance curve.}

\item{plot.folds}{Plot the individual folds as well.}

\item{verbose}{Control amount of screen reporting.}

\item{silent}{To allow running with no output for simplifying model).}

\item{keep.fold.models}{Keep the fold models from cross valiation.}

\item{keep.fold.vector}{Allows the vector defining fold membership to be kept.}

\item{keep.fold.fit}{Allows the predicted values for observations from CV to be kept.}

\item{...}{Allows for any additional plotting parameters.}
}
\value{
GBM models using gbm as the engine.
}
\description{
SD fork of dismo's gbm.step to add evaluation metrics like d.squared and rmse. J. Leathwick and
J. Elith - 19th September 2005, version 2.9. Function to assess optimal no of boosting trees
using k-fold cross validation. Implements the cross-validation procedure described on page 215 of
Hastie T, Tibshirani R, Friedman JH (2001) The Elements of Statistical Learning: Data Mining,
Inference, and Prediction Springer-Verlag, New York.
}
\details{
Divides the data into 10 subsets, with stratification by prevalence if required for pa data then
fits a gbm model of increasing complexity along the sequence from n.trees to n.trees + (n.steps *
step.size) calculating the residual deviance at each step along the way after each fold
processed, calculates the average holdout residual deviance and its standard error then
identifies the optimal number of trees as that at which the holdout deviance is minimised and
fits a model with this number of trees, returning it as a gbm model along with additional
information from the cv selection process.

D squared is 1 - (cv.dev / total.deviance). Abeare thesis: For each of the
fitted models, the pseudo-R2, or D2, or Explained Deviance, was calculated
for comparison, where: D2 = 1 – (residual deviance/total deviance).

requires gbm library from Cran
requires roc and calibration scripts of J Elith
requires calc.deviance script of J Elith/J Leathwick
}

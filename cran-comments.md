---
title: "cran-comments"
author: "Simon Dedman"
date: "23 August 2023"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

***

## Test environments
* local R installation, R 4.2.2, xubuntu 23.04
* win-builder (devel and release)

***

## R CMD check results

0 errors | 1 warnings | 2 notes

Warnings:

* "Rd files with duplicated alias 'gbm.auto': ‘gbm.auto-package.Rd’ ‘gbm.auto.Rd’". Have found no way to remove this, caused by lifecycle and usethis autogeneration of gbm.auto-package.R and Rd.

Notes:

* Possible code problems: no visible binding for global variables: named variables are column names in a csv exported by another function.

* Non-standard things, png file created by examples, doesn't exist permanently.

***

## Downstream dependencies

There are currently no downstream dependencies for this package

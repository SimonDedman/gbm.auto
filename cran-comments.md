---
title: "cran-comments"
author: "Simon Dedman"
date: "14 January 2021"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

***

## Test environments
* local R installation, R 4.0.3, xubuntu 20.10
* win-builder (devel and release)

***

## R CMD check results

0 errors | 0 warnings | 1 note

"checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    ‘COPYING.LESSERv3’ ‘GSHHG.zip’ ‘GSHHS_shp’ ‘LICENSE.TXT’ ‘README.TXT’
    ‘SHAPEFILES.TXT’ ‘WDBII_shp’"
These are shapefiles downloaded for the basemap function, into a tmp folder.
They are not present in the check directory after the check has finished.

* This is a new release.

***

## Downstream dependencies

There are currently no downstream dependencies for this package

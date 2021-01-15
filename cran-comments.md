---
title: "cran-comments"
author: "Simon Dedman"
date: "14 January 2021"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Test environments
* local linux install, xubuntu 20.10, R 4.0.3
* win-builder (devel and release)
* 

***

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking top-level files ... NOTE
  Non-standard file/directory found at top level: ‘cran-comments.md’
  ./Rbuildignore contains "^cran-comments\.Rmd$" but this still appears and I can't figure out why.

## Downstream dependencies

There are currently no downstream dependencies for this package

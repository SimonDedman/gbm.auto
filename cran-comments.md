---
title: "cran-comments"
author: "Simon Dedman"
date: "6 August 2018"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# From Hadley's blog:

## Test environments
* local OS X install, R 3.1.2
* ubuntu 12.04 (on travis-ci), R 3.1.2
* win-builder (devel and release)

It’s painful to manage multiple R versions, especially since you’ll need to reinstall all your packages. Instead, you can run R CMD check on CRAN’s servers with devtools::build_win(). This builds your package and submits it to the CRAN win-builder. 10-20 minutes after submission, you’ll receive an e-mail telling you the check results.

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking dependencies in R code ... NOTE
  Namespace in Imports field not imported from: 'R6'. R6 is a build-time dependency.
* Any NOTEs go in a bulleted list. For each NOTE, I include the message from R CMD check and a brief description of why I think it’s OK. If there were no NOTEs, I’d say “There were no ERRORs, WARNINGs or NOTEs”


## Downstream dependencies
I have also run R CMD check on downstream dependencies of httr 
(https://github.com/wch/checkresults/blob/master/httr/r-release). 
All packages that I could install passed except:

* Ecoengine: this appears to be a failure related to config on 
  that machine. I couldn't reproduce it locally, and it doesn't 
  seem to be related to changes in httr (the same problem exists 
  with httr 0.4).
  
If there are downstream dependencies, I run R CMD check on each package and summarise the results. If there are no downstream dependencies, keep this section, but say: “There are currently no downstream dependencies for this package”.

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gbm.factorplot.R
\name{gbm.factorplot}
\alias{gbm.factorplot}
\title{Creates ggplots of marginal effect for factorial variables from plot.gbm in gbm.auto.}
\usage{
gbm.factorplot(
  x,
  factorplotlevels = NULL,
  ggplot2guideaxisangle = 0,
  ggplot2labsx = "",
  ggplot2labsy = "Marginal Effect",
  ggplot2axistext = 1.5,
  ggplot2axistitle = 2,
  ggplot2legendtext = 1,
  ggplot2legendtitle = 1.5,
  ggplot2legendtitlealign = 0,
  ggplot2plotbackgroundfill = "white",
  ggplot2plotbackgroundcolour = "grey50",
  ggplot2striptextx = 2,
  ggplot2panelbordercolour = "black",
  ggplot2panelborderfill = NA,
  ggplot2panelborderlinewidth = 1,
  ggplot2legendspacingx = grid::unit(0, "cm"),
  ggplot2legendbackground = ggplot2::element_blank(),
  ggplot2panelbackgroundfill = "white",
  ggplot2panelbackgroundcolour = "grey50",
  ggplot2panelgridcolour = "grey90",
  ggplot2legendkey = ggplot2::element_blank(),
  ggsavefilename = paste0(lubridate::today(), "_Categorical-variable.png"),
  ggsaveplot = ggplot2::last_plot(),
  ggsavedevice = "png",
  ggsavepath = "",
  ggsavescale = 2,
  ggsavewidth = 10,
  ggsaveheight = 4,
  ggsaveunits = "in",
  ggsavedpi = 300,
  ggsavelimitsize = TRUE,
  ...
)
}
\arguments{
\item{x}{Input data.frame or tibble or csv (full file address including .csv) to read, must be a
categorical variable.}

\item{factorplotlevels}{Character vector of the variable's levels to reorder the x axis by, all
must match those in the first column of the csv exactly. Default NULL orders from high to low Y
value.}

\item{ggplot2guideaxisangle}{Default 0. Set at e.g. 90 to rotate.}

\item{ggplot2labsx}{Default: "".}

\item{ggplot2labsy}{Default: "Marginal Effect".}

\item{ggplot2axistext}{Default: 1.5.}

\item{ggplot2axistitle}{Default: 2.}

\item{ggplot2legendtext}{Default: 1.}

\item{ggplot2legendtitle}{Default: 1.5.}

\item{ggplot2legendtitlealign}{Default: 0, # otherwise effect type title centre aligned for some reason.}

\item{ggplot2plotbackgroundfill}{Default: "white", white background.}

\item{ggplot2plotbackgroundcolour}{Default: "grey50", background lines.}

\item{ggplot2striptextx}{Default: 2.}

\item{ggplot2panelbordercolour}{Default: "black".}

\item{ggplot2panelborderfill}{Default: NA.}

\item{ggplot2panelborderlinewidth}{Default: 1.}

\item{ggplot2legendspacingx}{Default: unit(0, "cm"), # compress spacing between legend items, this is min.}

\item{ggplot2legendbackground}{Default: ggplot2::element_blank().}

\item{ggplot2panelbackgroundfill}{Default: "white".}

\item{ggplot2panelbackgroundcolour}{Default: "grey50".}

\item{ggplot2panelgridcolour}{Default: "grey90".}

\item{ggplot2legendkey}{Default: ggplot2::element_blank().}

\item{ggsavefilename}{Default: paste0(saveloc, lubridate::today(), "_SankeyAlluvial_EMT.SoEv-EfTyp_Col-EfSz.png").}

\item{ggsaveplot}{Default: last_plot().}

\item{ggsavedevice}{Default: "png".}

\item{ggsavepath}{Default: "".}

\item{ggsavescale}{Default: 2.}

\item{ggsavewidth}{Default: 10.}

\item{ggsaveheight}{Default: 4.}

\item{ggsaveunits}{Default: "in".}

\item{ggsavedpi}{Default: 300.}

\item{ggsavelimitsize}{Default: TRUE.}

\item{...}{Allow params to be called from higher function esp gbm.auto.}
}
\value{
Factorial ggplot saved with users preferred location and name.
}
\description{
Creates an additional plot to those created by gbm.plot within gbm.auto. Can also take
Bin/Gaus_Best_line.csv or similar csvs directly. Allows changing of x axis levels and all ggplot
and ggsave params.
}
\details{
`r lifecycle::badge("experimental")
}
\author{
Simon Dedman, \email{simondedman@gmail.com}
}

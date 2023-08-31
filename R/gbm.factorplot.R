#' Creates ggplots of marginal effect for factorial variables from plot.gbm in gbm.auto.
#'
#' Creates an additional plot to those created by gbm.plot within gbm.auto. Can also take
#' Bin/Gaus_Best_line.csv or similar csvs directly. Allows changing of x axis levels and all ggplot
#' and ggsave params.
#'
#' `r lifecycle::badge("experimental")

#' @param x Input data.frame or tibble or csv (full file address including .csv) to read, must be a
#' categorical variable.
#' @param factorplotlevels Character vector of the variable's levels to reorder the x axis by, all
#' must match those in the first column of the csv exactly. Default NULL orders from high to low Y
#' value.
#' @param ggplot2guideaxisangle Default 0. Set at e.g. 90 to rotate.
#' @param ggplot2labsx Default: "".
#' @param ggplot2labsy Default: "Marginal Effect".
#' @param ggplot2axistext Default: 1.5.
#' @param ggplot2axistitle Default: 2.
#' @param ggplot2legendtext Default: 1.
#' @param ggplot2legendtitle Default: 1.5.
#' @param ggplot2legendtitlealign Default: 0, # otherwise effect type title centre aligned for some reason.
#' @param ggplot2plotbackgroundfill Default: "white", white background.
#' @param ggplot2plotbackgroundcolour Default: "grey50", background lines.
#' @param ggplot2striptextx Default: 2.
#' @param ggplot2panelbordercolour Default: "black".
#' @param ggplot2panelborderfill Default: NA.
#' @param ggplot2panelborderlinewidth Default: 1.
#' @param ggplot2legendspacingx Default: unit(0, "cm"), # compress spacing between legend items, this is min.
#' @param ggplot2legendbackground Default: ggplot2::element_blank().
#' @param ggplot2panelbackgroundfill Default: "white".
#' @param ggplot2panelbackgroundcolour Default: "grey50".
#' @param ggplot2panelgridcolour Default: "grey90".
#' @param ggplot2legendkey Default: ggplot2::element_blank().
#' @param ggsavefilename Default: paste0(saveloc, lubridate::today(), "_SankeyAlluvial_EMT.SoEv-EfTyp_Col-EfSz.png").
#' @param ggsaveplot Default: last_plot().
#' @param ggsavedevice Default: "png".
#' @param ggsavepath Default: "".
#' @param ggsavescale Default: 2.
#' @param ggsavewidth Default: 10.
#' @param ggsaveheight Default: 4.
#' @param ggsaveunits Default: "in".
#' @param ggsavedpi Default: 300.
#' @param ggsavelimitsize Default: TRUE.
#' @param ... Allow params to be called from higher function esp gbm.auto.
#'
#' @return Factorial ggplot saved with users preferred location and name.
#'
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#'
#' @export
#' @importFrom dplyr arrange mutate rename
#' @importFrom ggplot2 aes element_blank element_line element_rect element_text geom_col ggplot ggsave guide_axis guides labs last_plot rel theme theme_minimal %+replace%
#' @importFrom grid unit
#' @importFrom lubridate today
#' @importFrom readr read_csv
#' @importFrom tidyselect last_col

# x <- "/home/simon/Documents/Si Work/PostDoc Work/Gbmauto help/2022-03_Frances_Naomi/2023-02-21 categorical variables issue/count_elasmo/Gaus_Best_line_isl_grp.csv"
# factorplotlevels <- c("west tuamotu",
#                        "east tuamotu",
#                        "leeward",
#                        "windward",
#                        "australes",
#                        "marquesas")

gbm.factorplot <- function(x,
                           factorplotlevels = NULL,
                           ggplot2guideaxisangle = 0,
                           ggplot2labsx = "",
                           ggplot2labsy = "Marginal Effect",
                           ggplot2axistext = 1.5,
                           ggplot2axistitle = 2,
                           ggplot2legendtext = 1,
                           ggplot2legendtitle = 1.5,
                           ggplot2legendtitlealign = 0, # otherwise effect type title centre aligned for some reason
                           ggplot2plotbackgroundfill = "white", # white background
                           ggplot2plotbackgroundcolour = "grey50",
                           ggplot2striptextx = 2,
                           ggplot2panelbordercolour = "black",
                           ggplot2panelborderfill = NA,
                           ggplot2panelborderlinewidth = 1,
                           ggplot2legendspacingx = grid::unit(0, "cm"), # compress spacing between legend items, this is min
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
) {
  if (is.character(x)) {
    x <- readr::read_csv(x) # readr::read_csv("/home/simon/Documents/Si Work/PostDoc Work/Gbmauto help/2022-03_Frances_Naomi/2023-02-21 categorical variables issue/count_elasmo/Gaus_Best_line_isl_grp.csv")
  } else {
    x <- x
  }

  # add check to see whether csv is categorical: first column
  # if (class(x[, 1][[1]]) != "character") stop("csv is not a categorical/factorial variable")
  # wasn't working, should debug but since this is already checked by gbm.auto before this is run,
  # and gbm.factorplot is going to be used by gbm.auto 99.99% of it's life, just remove for now.

  x <- x |>
    # dplyr::mutate(ycentred = y - mean(y)) |> # make mean column
    # The above exists already in gbm.auto L843
    # Also it uses preset colname y produced by gbm.auto.
    dplyr::rename("Category" = tidyselect::last_col(offset = 2), # no first_col option
                  "ycentred" = tidyselect::last_col()) |> # attempt to address no visible binding for global variable ‘Category’
    # 2023-08-30 quoted Category to hopefully address gbm.factorplot: no visible binding for global variable ‘Category’
    # also need to do for ycentred
    dplyr::arrange(ycentred) # re-order the x axis for categorical variables in order from high to low value

  # re-reorder by factorplotlevels if present
  if (!is.null(factorplotlevels)) {
    # check level number matches x categories
    if (length(factorplotlevels) != length(x[, "Category"])) stop("number of factorplotlevels doesn't match number of categories in x") # x$Category
    # check names match
    if (!all(factorplotlevels %in% x[, "Category"])) stop(paste0("The following level names not present in categories in x: ", # x$Category
                                                            factorplotlevels[!factorplotlevels %in% x[, "Category"]])) # x$Category
    x <- x |>
      dplyr::mutate(Category = ordered(Category, levels = factorplotlevels)) |> # recreate Category as ordered factor with factorplotlevels. Can't quote "Category" in ordered
      dplyr::arrange(Category) # arrange by labels (implicit). Can't quote "Category"
  } else {
    if (is.integer(x$Category)) x <- x |> dplyr::mutate(Category = ordered(Category)) # make integers ordered factor to avoid defaulting to continuous
  }

  ggplot2::ggplot(x) +
    ggplot2::geom_col(ggplot2::aes(x = Category, # x[, "Category"]
                                   y = ycentred)) + # x[, "ycentred"]
    # rotate x axis labels
    ggplot2::guides(x =  ggplot2::guide_axis(angle = ggplot2guideaxisangle)) +
    # alter axis labels
    ggplot2::labs(x = "",
                  y = "Marginal Effect") +
    # change theme elements
    ggplot2::theme_minimal() %+replace% ggplot2::theme(
      axis.text = ggplot2::element_text(size = ggplot2::rel(ggplot2axistext)),
      axis.title = ggplot2::element_text(size = ggplot2::rel(ggplot2axistitle)),
      legend.text = ggplot2::element_text(size = ggplot2::rel(ggplot2legendtext)),
      legend.title = ggplot2::element_text(size = ggplot2::rel(ggplot2legendtitle)),
      legend.title.align = ggplot2legendtitlealign, # otherwise effect type title centre aligned for some reason
      plot.background = ggplot2::element_rect(fill = ggplot2plotbackgroundfill,
                                              colour = ggplot2plotbackgroundcolour), # white background
      strip.text.x = ggplot2::element_text(size = ggplot2::rel(ggplot2striptextx)),
      panel.border = ggplot2::element_rect(colour = ggplot2panelbordercolour,
                                           fill = ggplot2panelborderfill,
                                           linewidth = ggplot2panelborderlinewidth),
      legend.spacing.x = ggplot2legendspacingx, # compress spacing between legend items, this is min
      legend.background = ggplot2legendbackground,
      panel.background = ggplot2::element_rect(fill = ggplot2panelbackgroundfill,
                                               colour = ggplot2panelbackgroundcolour),
      panel.grid = ggplot2::element_line(colour = ggplot2panelgridcolour),
      legend.key = ggplot2legendkey)

  # save
  ggplot2::ggsave(filename = ggsavefilename,
                  plot = ggsaveplot,
                  device = ggsavedevice,
                  path = ggsavepath,
                  scale = ggsavescale,
                  width = ggsavewidth,
                  height = ggsaveheight,
                  units = ggsaveunits,
                  dpi = ggsavedpi,
                  limitsize = ggsavelimitsize)
}

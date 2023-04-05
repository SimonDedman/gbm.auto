# 2023-03-10 code to plot categorical variables with ggplot
# Todo: could automate the process of scanning through the output folder,
# looking for csvs with certain naming conventions,
# screening out those which aren't factorial variables,
# then auto-plotting the rest?
# add abline @ Y=0

# library(readr)
# library(utils)
library(ggplot2)
# library(gplots)
# library(rstatix)
library(ggpubr)
# library(viridis)
remotes::install_github("wilkelab/ungeviz")
library(ungeviz)
library(lubridate)
library(stringi)

options(scipen = 5)
# readloc = "/home/simon/Documents/Si Work/PostDoc Work/Saving The Blue/Projects/2021-10_Drumline_Reefshark/BRT/CPUE_Loops/1/CaribbeanReef_CPUE/"
# readfile = "Gaus_Best_line_Habitat.csv"
# saveloc = "/home/simon/Documents/Si Work/PostDoc Work/Saving The Blue/Projects/2021-10_Drumline_Reefshark/BRT/CPUE_Loops/1/CaribbeanReef_CPUE/"
readloc = NULL
readfile = NULL
saveloc = NULL
savename = stringi::stri_replace_all_fixed(str = readfile,
                        ".csv",
                        ".png")

readr::read_csv(paste0(readloc, readfile)) |>
  dplyr::mutate(ycentred = y - mean(y)) |> # centre y conforming to plot.gbm format
  dplyr::rename(Category = last_col(offset = 2)) |> # Rename first column to a generic name. No first_col option
  # dplyr::case_match(Category, # uncomment this to change the category values to new names if desired.
  #                   "OLD NAME 1" ~ "NEW NAME 1",
  #                   "OLD NAME 2" ~ "NEW NAME 2",
  #                   "OLD NAME 3" ~ "NEW NAME 3",
  #                   .default = Category) |>
  # make category an ordered factor. Default alphabetical, can change as desired.
  dplyr::mutate(Category = ordered(Category, levels = Category)) |>
  # re-order the x axis for categorical variables so they're in order from high to low or can add factor levels
  dplyr::arrange(ycentred) |>
  ggplot2::ggplot(ggplot2::aes(x = Category,
                               y = ycentred)) +
  ungeviz::geom_hpline() +
  ggplot2::labs(x = "",
       y = "Marginal Effect") +
  ggpubr::theme_pubclean(base_size = 16) +
  ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 45))

ggplot2::ggsave(paste0(saveloc, lubridate::today(), "_", savename),
                plot = ggplot2::last_plot(), device = "png", scale = 1.75, width = 7,
                height = 4, units = "in", dpi = 300, limitsize = TRUE)

# library(tidyverse)

# 2023-03-08: work in progress, commented out to avoid problems in build

# as a function:
# levels param: default Alphabetical, other options Rising (see below), Falling (reverse, arrange(desc)), then user entry - detect a character vector & check it matches the Category entries
# save options?
# basically people could just edit this themselves for the most part?
# do the basic plot then save as an rds object?
# option for turning x axis labels at an angle

# tmp <-
# read_csv("/home/simon/Documents/Si Work/PostDoc Work/Gbmauto help/2022-03_Frances_Naomi/2023-02-21 categorical variables issue/count_elasmo/Gaus_Best_line_isl_grp.csv") |>
#   dplyr::mutate(ycentred = y - mean(y)) |>
#   dplyr::rename(Category = last_col(offset = 2)) |> # no first_col option
#   # alter axis titles
#   dplyr::mutate(Category = case_match(Category,
#                                       "australes" ~ "Australes",
#                                       # add more here. Or do this in csv
#                                       .default = Category)) |>
#   # re-order the x axis for categorical variables so they are either
#   # 1) in order from high to low or
#   dplyr::arrange(ycentred) |>
#   dplyr::mutate(Category = ordered(Category, levels = Category)) |>
#   # 2) in some other logical order I specify such as factor level, rather than alphabetical.
#   # Do this by setting the levels of Category. I did this by arranging by ycentred (=high to low Y value)
#   # Default is alphabetical, so if you want that the previous 2 lines can be an if statement
#   # And if you want something else then you can manually set the category levels
#
#   ggplot() +
#   geom_col(aes(x = Category,
#                y = ycentred)) +
#   # alter axis labels
#   labs(x = "",
#        y = "Marginal Effect") +
#   theme_minimal() %+replace% theme(
#     axis.text = element_text(size = rel(1.5)),
#     axis.title = element_text(size = rel(2)),
#     legend.text = element_text(size = rel(1)),
#     legend.title = element_text(size = rel(1.5)),
#     legend.title.align = 0, # otherwise effect type title centre aligned for some reason
#     plot.background = element_rect(fill = "white", colour = "grey50"), # white background
#     strip.text.x = element_text(size = rel(2)),
#     panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
#     legend.spacing.x = unit(0, "cm"), # compress spacing between legend items, this is min
#     legend.background = element_blank(),
#     panel.background = element_rect(fill = "white", colour = "grey50"),
#     panel.grid = element_line(colour = "grey90"),
#     legend.key = element_blank())


# ggsave(paste0(saveloc, today(), "_SankeyAlluvial_EMT.SoEv-EfTyp_Col-EfSz.png"),
#        plot = last_plot(), device = "png", path = "",
#        scale = 2, width = 10, height = 4, units = "in", dpi = 300, limitsize = TRUE)

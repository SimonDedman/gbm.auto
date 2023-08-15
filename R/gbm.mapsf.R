#' Maps of predicted abundance from Boosted Regression Tree modelling
#'
#' Generates maps from the outputs of gbm.step then Gbm.predict.grids, handled automatically within
#' gbm.auto but can be run alone, and generates representativeness surfaces from the output of
#' gbm.rsb.
#'
#' @author Simon Dedman, \email{simondedman@@gmail.com}

#' @import ggplot2
#' @importFrom ggmap get_map ggmap register_google
#' @importFrom ggspatial layer_spatial
#' @importFrom lubridate today
#' @importFrom stars geom_stars st_as_stars
#' @importFrom sf st_set_crs st_bbox st_transform st_as_sfc st_as_sf st_buffer
#' @importFrom starsExtra trim2
#' @importFrom viridis scale_fill_viridis
#' @export
#'
#' @param predabund Predicted abundance data frame produced by gbm.auto (Abundance_Preds_only.csv),
#' with Latitude, Longitude, and Predicted Abundance columns. Default NULL. You need to read the csv
#'  in R if not already present as an object in the environment.
#' @param predabundlon Longitude column number. Default 2.
#' @param predabundlat Latitude column number. Default 1.
#' @param predabundpreds Predicted abundance column number, default 3.
#' @param myLocation Location for extents, format c(xmin, ymin, xmax, ymax). Default NULL, extents
#' autocreated from data.
#' @param trim Remove NA & <=0 values and crop to remaining date extents? Default TRUE.
#' @param scale100 Scale Predicted Abundance to 100? Default FALSE.
#' @param gmapsAPI Enter your Google maps API here, quoted character string. Default NULL.
#' @param mapsource Source for ggmap::get_map; uses Stamen as fallback if no Google Maps API present
#' . Options: "google", "stamen", "gbm.basemap". Default "google". Using "gbm.basemap" requires one
#' to have run that functiuon already, and enter its location using the shape paramater below.
#' @param googlemap If pulling basemap from Google maps, this sets expansion factors since Google
#' Maps tiling zoom setup doesn't align to myLocation extents. Default TRUE.
#' @param maptype = "satellite", # Type of map for ggmap::get_map. Options: "terrain",
#' "terrain-background", "satellite", "roadmap", "hybrid", "toner", "terrain-labels",
#' "terrain-lines", "toner-2010", "toner-2011", "toner-background", "toner-hybrid", "toner-labels",
#'  "toner-lines", "toner-lite". Google options: “terrain”, “satellite”, “roadmap”, “hybrid”.
#' @param darkenproportion Amount to darken the google/stamen basemap, 0-1. Default 0.
#' @param mapzoom Highest number = zoomed in. Google: 3 (continent) - 21 (building). stamen: 0-18.
#' Default 9.
#' @param shape If mapsource is "gbm.basemap", enter the full path to gbm.basemaps downloaded map,
#' typically Crop_Map.shp, including the .shp. Default NULL.
#' @param colourscale Scale fill colour scheme to use, default "viridis", other option is
#' "gradient".
#' @param colorscale Scale fill colour scheme to use, default NULL, populating this will overwrite
#' colourscale.
#' @param heatcolours Vector of colours if gradient selected for colourscale, defaults to heatmap
#' theme.
#' @param colournumber Number of colours to spread heatcolours over, if gradient selected for
#' colourscale. Default 8.
#' @param colourscalelimits Colour scale limits, default NULL, vector of 2, e.g. c(0, 0).
#' @param colourscalebreaks Colour scale breaks, default NULL.
#' @param colourscalelabels Colour scale labels, default NULL, must match number of breaks.
#' @param colourscaleexpand Colour scale expand, default NULL, vector of 2, e.g. c(0, 0).
#' @param studyspecies Name of your study species, appears in plot title and savename. Default
#' "MySpecies".
#' @param plottitle Title of the resultant plot, default paste0("Predicted abundance of ",
#' studyspecies).
#' @param plotsubtitle Plot subtitle, default ""CPUE". Can add the n of your individuals.
#' @param legendtitle Legend title, default "CPUE".
#' @param plotcaption Plot caption, default "gbm.auto::gbm.map" + today's date.
#' @param axisxlabel Default "Longitude".
#' @param axisylabel Default "Latitude".
#' @param legendposition Vector of 2, format c(1,2), Proportional distance of (middle?) of legend
#' box from L to R, percent distance from Bottom to Top. Values 0 to 1. Default c(0.05, 0.15).
#' @param fontsize Font size, default 12.
#' @param fontfamily = Font family, default "Times New Roman".
#' @param filesavename File savename, default today's date + studyspecies + legendtitle.
#' @param savedir Save outputs to a temporary directory (default) else change to current directory
#' e.g. "/home/me/folder". Do not use getwd() here. No terminal slash. E.g.
#' paste0(movegroupsavedir, "Plot/") .
#' @param receiverlats Vector of latitudes for receivers to be plotted.
#' @param receiverlons Vector of longitudes for receivers to be plotted. Same length as receiverlats.
#' @param receivernames Vector of names for receivers to be plotted. Same length as receiverlats.
#' @param receiverrange Single (will be recycled), or vector (same length as receiverlats) of
#' detection ranges in metres for receivers to be plotted. If you have a max and a (e.g.) 90 percent
#'  detection range, probably use max.
#' @param recpointscol Colour of receiver centrepoint outlines. Default "black".
#' @param recpointsfill Colour of receiver centrepoint fills. Default "white".
#' @param recpointsalpha Alpha value of receiver centrepoint fills, 0 (invisible) to 1 (fully
#' visible). Default 0.5.
#' @param recpointssize Size of receiver points. Default 1.
#' @param recpointsshape Shape of receiver points, default 21, circle with outline and fill.
#' @param recbufcol Colour of the receiver buffer circle outlines. Default "grey75"
#' @param recbuffill Colour of the receiver buffer circle fills. Default "grey".
#' @param recbufalpha Alpha value of receiver buffer fills, 0 (invisible) to 1 (fully visible).
#' Default 0.5.
#' @param reclabcol Receiver label text colour. Default "black".
#' @param reclabfill Receiver label fill colour, NA for no fill. Default NA.
#' @param reclabnudgex Receiver label offset nudge in X dimension. Default 0.
#' @param reclabnudgey Receiver label offset nudge in Y dimension. Default -200.
#' @param reclabpad Receiver label padding in lines. Default 0.
#' @param reclabrad Receiver label radius in lines. Default 0.15.
#' @param reclabbord Receiver label border in mm. Default 0.
#'
#' @return Species abundance maps using data provided by gbm.auto, and
#' Representativeness Surface Builder maps using data provided by gbm.rsb, to be
#' run in a png/par/gbm.map/dev.off sequence.
#'
#' @details
#'
#' Error in seq.default(xlim[1], xlim[2], by = byx):wrong sign in 'by' argument
#' Check that your lat & long columns are the right way around. Ensure grids (predabund) data are
#' gridded, i.e. they are in a regular pattern of same/similar lines of lat/lon, even if they're
#' missing sections.
#'
#' Suggested parameter values:
#' z = rsbdf[,"Unrepresentativeness"]
#'
#' mapmain = "Unrepresentativeness: "
#'
#' legendtitle = "UnRep 0-1"
#'
#' How to get Google map basemaps:
#'
#' @examples
#' \donttest{
#' # Not run
#' }

gbm.mapsf <- function(
    predabund = NULL, # predicted abundance data frame produced by gbm.auto (Abundance_Preds_only.csv), with Latitude, Longitude, and Predicted Abundance columns.
    predabundlon = 2, # Longitude column number.
    predabundlat = 1, # Latitude column number.
    predabundpreds = 3, # Predicted abundance column number.
    myLocation = NULL, # location for extents, format c(xmin, ymin, xmax, ymax).
    # Default NULL, extents autocreated from data.
    # c(-79.3, 25.68331, -79.24, 25.78)
    trim = TRUE, # remove NA & 0 values and crop to remaining date extents? Default TRUE.
    scale100 = FALSE, # scale Predicted Abundance to 100? Default FALSE.
    gmapsAPI = NULL, # enter your Google maps API here, quoted character string
    mapsource = "google", # Source for ggmap::get_map; uses Stamen as fallback if no Google Maps API present. Options: "google", "stamen", "gbm.basemap".
    googlemap = TRUE, # If pulling basemap from Google maps, this sets expansion factors since
    # Google Maps tiling zoom setup doesn't align to myLocation extents.
    maptype = "satellite", # Type of map for ggmap::get_map. Options: "terrain", "terrain-background", "satellite", "roadmap", "hybrid", "toner",
    #  "terrain-labels", "terrain-lines", "toner-2010", "toner-2011",
    # "toner-background", "toner-hybrid", "toner-labels", "toner-lines", "toner-lite". Google options: “terrain”, “satellite”, “roadmap”, “hybrid”.
    darkenproportion = 0, # amount to darken the basemap, 0-1.
    mapzoom = 9, # google: 3 (continent) - 21 (building). stamen: 0-18
    shape = NULL, # If mapsource is "gbm.basemap", enter the full path to gbm.basemaps downloaded map, typically Crop_Map.shp, including the .shp.
    expandfactor = 0, # extents expansion factor for basemap. default was 1.6
    colourscale = "viridis", # Scale fill colour scheme to use, default "viridis", other option is "gradient".
    colorscale = NULL, # Scale fill colour scheme to use, default NULL, populating this will overwrite colourscale.
    heatcolours = c("white", "yellow", "orange","red", "brown4"), # Vector of colours if gradient selected for colourscale, defaults to heatmap theme.
    colournumber = 8, # Number of colours to spread heatcolours over, if gradient selected for colourscale. Default 8.
    colourscalelimits = NULL, # Colour scale limits, default NULL, vector of 2, e.g. c(0, 0).
    colourscalebreaks = NULL, # Colour scale breaks, default NULL.
    colourscalelabels = NULL, # Colour scale labels, default NULL, must match number of breaks.
    colourscaleexpand = NULL, # Colour scale expand, default NULL, vector of 2, e.g. c(0, 0).
    studyspecies = "MySpecies", # The name of your study species, appears in plot title and savename.
    plottitle = paste0("Predicted abundance of ", studyspecies),
    plotsubtitle  = "CPUE", # Plot subtitle. Can add the n of your individuals.
    legendtitle = "CPUE",
    plotcaption = paste0("gbm.auto::gbm.map, ", lubridate::today()),
    axisxlabel = "Longitude",
    axisylabel = "Latitude",
    legendposition = c(0.05, 0.15), # Percent distance (of middle? of legend box) from L to R, percent distance from Bottom to Top.
    fontsize = 12,
    fontfamily = "Times New Roman",
    filesavename = paste0(lubridate::today(), "_", studyspecies, "_", legendtitle, ".png"),
    savedir = tempdir(), # file.path(work.dir, out.dir, "Scaled")
    receiverlats = NULL, # vector of latitudes for receivers to be plotted
    receiverlons = NULL, # vector of longitudes for receivers to be plotted
    receivernames = NULL, # vector of names for receivers to be plotted
    receiverrange = NULL, # single (will be recycled), or vector of detection ranges in metres for receivers to be plotted
    recpointscol = "black", # Colour of receiver centrepoint outlines.
    recpointsfill = "white", # Colour of receiver centrepoint fills.
    recpointsalpha = 0.5, # Alpha value of receiver centrepoint fills, 0 (invisible) to 1 (fully visible).
    recpointssize = 1, # Size of receiver points.
    recpointsshape = 21, # Shape of receiver points, default 21, circle with outline and fill.
    recbufcol = "grey75", # Colour of the receiver buffer circle outlines.
    recbuffill = "grey", # Colour of the receiver buffer circle fills.
    recbufalpha = 0.5,  # Alpha value of receiver buffer fills, 0 (invisible) to 1 (fully visible).
    reclabcol = "black", # Receiver label text colour.
    reclabfill = NA, # Receiver label fill colour, NA for no fill.
    reclabnudgex = 0, # Receiver label offset nudge in X dimension.
    reclabnudgey = -200, # Receiver label offset nudge in Y dimension.
    reclabpad = 0, # Receiver label padding in lines.
    reclabrad = 0.15, # Receiver label radius in lines.
    reclabbord = 0 # Receiver label border in mm.
){
  # Todo:
  # unrep exponential breaks for RSB legend: breaks = exp01seq
  # gbm.auto L1847
  # linear01seq <- seq(from = 0, to = 1, length.out = 9) #linear sequence from 0:1, 9 bins
  # exp01seq <- expm1(4*linear01seq)/expm1(4) # exponentiate to change shape then scale back to 1

  # gbm.auto allow entry for user preference elements, google, etc.

  # remove Hans as gbm.auto author since the stuff he contributed to is gone?

  # check receiver inputs are the correct lengths, if present.
  if (!is.null(receiverlats) & !is.null(receiverlons)) if (length(receiverlats) != length(receiverlons)) stop("length of receiverlats must equal length of receiverlons")
  if (!is.null(receiverlats) & !is.null(receivernames)) if (length(receiverlats) != length(receivernames)) stop("length of receivernames must equal length of receiverlats/lons")
  if (!is.null(receiverlats) & !is.null(receiverrange)) if (length(receiverrange) != length(receiverlons)) if (length(receiverrange) != 1) stop("length of receiverrange must equal length of receiverlats/lons, or 1")

  # Create receiver objects
  if (!is.null(receiverlats) & !is.null(receiverlons)) {
    receiver <- data.frame(lon = receiverlons,
                           lat = receiverlats)
    receiver <- sf::st_as_sf(receiver, coords = c("lon","lat"))  |>
      sf::st_set_crs(4326) |>
      sf::st_transform(3857)
    if (!is.null(receivernames)) {
      receiver <- cbind(receiver, receivernames)
    }
    if (!is.null(receiverrange)) {
      receiver <- cbind(receiver, receiverrange)
    }
  }

  # location for extents, format c(xmin, ymin, xmax, ymax). Default NULL, extents autocreated from data.
  myLocation <- c(min(predabund[,predabundlon], na.rm = TRUE),
                  min(predabund[,predabundlat], na.rm = TRUE),
                  max(predabund[,predabundlon], na.rm = TRUE),
                  max(predabund[,predabundlat], na.rm = TRUE))


  if (!is.null(gmapsAPI)) ggmap::register_google(key = gmapsAPI, # an api key
                                                 account_type = "standard",
                                                 write = TRUE)

  if (mapsource != "google") googlemap <- FALSE else googlemap <- TRUE # in case user forgot to set both the same

  if (expandfactor != 0) { # grow bounds extents if requested
    xmid <- mean(myLocation[c(1,3)])
    ymid <- mean(myLocation[c(2,4)])
    xmax <- ((myLocation[3] - xmid) * expandfactor) + xmid #updated for sf/st
    xmin <- xmid - ((xmid - myLocation[1]) * expandfactor)
    ymax <- ((myLocation[4] - ymid) * expandfactor) + ymid
    ymin <- ymid - ((ymid - myLocation[2]) * expandfactor)
    myLocation <- c(xmin, ymin, xmax, ymax)
  }

  if (googlemap) myLocation <- c(mean(c(myLocation[1], myLocation[3])), mean(c(myLocation[2], myLocation[4]))) # googlemap needs a center lon lat

  # mapsource
  if (mapsource == "gbm.basemap") {
    shape = st_read(dsn = shape)

    # autoheight <- (6 / (attr(myMap, "bb")[[4]] - attr(myMap, "bb")[[2]])) * (attr(myMap, "bb")[[3]] - attr(myMap, "bb")[[1]]) * 1.2
    # if (googlemap) autoheight <- 6.4 # googlemap pulls tiles for a centre point hence will always be square. But needs a bit extra for title area.
    autoheight <- 6.4 # dummy
    #autoheight####
  } else { # mapsource ifelse
    myMap <- ggmap::get_map(
      location = myLocation, # -62.57564  28.64368  33.78889  63.68533 # stamen etc want a bounding box
      zoom = mapzoom, # 3 (continent) - 21 (building). Stamen: 0-18
      # scale = "auto", # default "auto", 1, 2, 4 all the same
      messaging = TRUE, #
      source = mapsource, # "google" # using stamen as fallback
      maptype = maptype, # "satellite"
      crop = TRUE # google maps crs = 4326
    )

    # Define a function to fix the bbox to be in EPSG:3857
    # https://stackoverflow.com/a/50844502/1736291 Fixes "error no lon value" in ggmap below
    ggmap_bbox <- function(map) {
      if (!inherits(map, "ggmap")) stop("map must be a ggmap object")
      # Extract the bounding box (in lat/lon) from the ggmap to a numeric vector,
      # and set the names to what sf::st_bbox expects:
      map_bbox <- setNames(unlist(attr(map, "bb")), c("ymin", "xmin", "ymax", "xmax"))
      # Convert the bbox to an sf polygon, transform it to 3857,
      # and convert back to a bbox (convoluted, but it works)
      bbox_3857 <- sf::st_bbox(sf::st_transform(sf::st_as_sfc(sf::st_bbox(map_bbox, crs = 4326)), 3857))
      # Overwrite the bbox of the ggmap object with the transformed coordinates
      attr(map, "bb")$ll.lat <- bbox_3857["ymin"]
      attr(map, "bb")$ll.lon <- bbox_3857["xmin"]
      attr(map, "bb")$ur.lat <- bbox_3857["ymax"]
      attr(map, "bb")$ur.lon <- bbox_3857["xmax"]
      map
    }
    myMap <- ggmap_bbox(myMap) # Use the function. Resulting map is CRS 3857

    # Automate width * height adjustments for different map extent / ratio
    # 6 (manually chosen width, below), divided by width range times by height range
    # Maintains ratio by scales height to width(6). Then *1.2 because it still wasn't perfect.
    # attr(myMap, "bb")[[4]] - attr(myMap, "bb")[[2]] # longitude, x, width, bind as 6
    # attr(myMap, "bb")[[3]] - attr(myMap, "bb")[[1]] # latitude, y, height
    autoheight <- (6 / (attr(myMap, "bb")[[4]] - attr(myMap, "bb")[[2]])) * (attr(myMap, "bb")[[3]] - attr(myMap, "bb")[[1]]) * 1.2
    if (googlemap) autoheight <- 6.4 # googlemap pulls tiles for a centre point hence will always be square. But needs a bit extra for title area.
  } # close mapsource ifelse

  # remove extra columns from predabund so they don't become additional layers in stars raster
  predabund <- data.frame(Latitude = predabund[, predabundlat],
                          Longitude = predabund[, predabundlon],
                          PredAbund = predabund[, predabundpreds])
  # create stars from predabund
  predabundstars <- stars::st_as_stars(predabund, coords = c(predabundlon, predabundlat)) |>
    sf::st_set_crs(4326) # one of (i) character: a string accepted by GDAL, (ii) integer, a valid EPSG value (numeric), or (iii) an object of class crs.

  if (trim) { # trim raster extent to data?
    is.na(predabundstars[[1]]) <- predabundstars[[1]] == 0 # replace char pattern (0) in whole df/tbl with NA
    is.na(predabundstars[[1]]) <- predabundstars[[1]] < (max(predabundstars[[1]], na.rm = TRUE) * 0.05) # replace anything < 95% contour with NA since it won't be drawn
  }
  predabundstars <- starsExtra::trim2(predabundstars) # remove NA columns, which were all zero columns. This changes the bbox accordingly
  if (scale100) predabundstars[[1]] <- (predabundstars[[1]] / max(predabundstars[[1]], na.rm = TRUE)) * 100 # convert from raw values to 0:100 scale so legend is 0:100%

  # set colourscale from colorscale
  if (!is.null(colorscale)) colourscale <- colorscale

  # ggplot() +   #plot lines by year
  #   ggspatial::layer_spatial(shape, col = "black") +
  #   stars::geom_stars(data = predabundstars |> sf::st_transform(3857), inherit.aes = FALSE)
  # # don't need to see the rest, geom_stars & afters work outside of {ifelse}, not within: stops after later_spatial
  # # if order flipped (ggmap then layerspatial: layerspatial works)


  # plot map

  # if (mapsource != "gbm.basemap") {
  #   ggmap::ggmap(myMap, # basemap CRS = 3857? Doesn't plot on its own, don't worry
  #                darken = c(darkenproportion, "black"))
  # } + # close mapsource if2

  # if (mapsource == "gbm.basemap") {
  #   ggplot() +   #plot lines by year
  #     ggspatial::layer_spatial(shape, col = "black")
  # } + # close mapsource if

  if (mapsource == "gbm.basemap") {
    p <- ggplot() +   #plot lines by year
      ggspatial::layer_spatial(shape, col = "black")
  } else {
    p <- ggmap::ggmap(myMap, # basemap CRS = 3857? Doesn't plot on its own, don't worry
                      darken = c(darkenproportion, "black"))
  } # close mapsource ifelse

  p <- p +
    # Prediction surface
    # need to convert points to raster using byx code from gbm.map
    stars::geom_stars(data = predabundstars |> sf::st_transform(3857), inherit.aes = FALSE) +

    # # receiver centrepoints
    {if (!is.null(receiverlats) & !is.null(receiverlons))
      ggplot2::geom_sf(data = receiver |>
                         sf::st_transform(3857), # Vector transform after st_contour
                       # already 3857 above so converting twice but it ain't broke
                       colour = recpointscol,
                       fill = recpointsfill,
                       alpha = recpointsalpha,
                       size = recpointssize,
                       shape = recpointsshape,
                       inherit.aes = FALSE,
      )
    } +

    # receiver buffer circles
    {if (!is.null(receiverlats) & !is.null(receiverlons) & !is.null(receiverrange))
      ggplot2::geom_sf(data = sf::st_buffer(receiver, dist = receiverrange)  |>
                         sf::st_transform(3857), # Vector transform after st_contour
                       # already 3857 above so converting twice but it ain't broke
                       colour = recbufcol,
                       fill = recbuffill,
                       alpha = recbufalpha,
                       inherit.aes = FALSE
      )
    } +

    # receiver labels
    {if (!is.null(receiverlats) & !is.null(receiverlons) & !is.null(receivernames))
      ggplot2::geom_sf_label(data = receiver  |>
                               sf::st_transform(3857), # Vector transform after st_contour
                             # already 3857 above so converting twice but it ain't broke
                             colour = reclabcol,
                             fill = reclabfill,
                             inherit.aes = FALSE,
                             nudge_x = reclabnudgex,
                             nudge_y = reclabnudgey,
                             label.padding = unit(reclabpad, "lines"), # 0.25
                             label.r = unit(reclabrad, "lines"),
                             label.size = reclabbord, # 0.25
                             ggplot2::aes(label = receivernames)
      )
    } +

    # CPUE scale colours
    {if (colourscale == "viridis")
      viridis::scale_fill_viridis(
        if (!is.null(colourscalelimits)) {limits = colourscalelimits}, # c(0, 1),
        if (!is.null(colourscalebreaks)) {breaks = colourscalebreaks},
        if (!is.null(colourscalelabels)) {labels = colourscalelabels},
        if (!is.null(colourscaleexpand)) {expand = colourscaleexpand}, # c(0, 0),
        alpha = 1, # 0:1
        begin = 0, # hue
        end = 1, # hue
        direction = 1, # colour order, 1 or -1
        discrete = FALSE, # false = continuous
        option = "D", # A magma B inferno C plasma D viridis E cividis F rocket G mako H turbo
        space = "Lab",
        na.value = "grey50",
        guide = "colourbar",
        aesthetics = "fill",
        # name = waiver(),
        name = legendtitle,
        # limits = NA,
        # position = "left"
        position = "right"
      )
    } +

    # CPUE scale colours
    {if (colourscale == "gradient")
      scale_fill_gradientn(
        # if (!is.null(colourscalelimits)) limits = colourscalelimits, # c(0, 1),
        # if (!is.null(colourscalebreaks)) breaks = colourscalebreaks,
        # if (!is.null(colourscalelabels)) labels = colourscalelabels,
        # if (!is.null(colourscaleexpand)) expand = colourscaleexpand, # c(0, 0),
        name = legendtitle,
        position = "right",
        colours = colorRampPalette(heatcolours)(colournumber), # Vector of colours to use for n-colour gradient.
        na.value = "grey50"
      )
    } +

    ggplot2::ggtitle(plottitle, subtitle = plotsubtitle) +
    ggplot2::labs(x = axisxlabel, y = axisylabel, caption = plotcaption) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = legendposition, #%dist (of middle? of legend box) from L to R, %dist from Bot to Top
      legend.spacing.x = ggplot2::unit(0, 'cm'), #compress spacing between legend items, this is min
      legend.spacing.y = ggplot2::unit(0, 'cm'), #compress spacing between legend items, this is min
      legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 8),
      legend.background = ggplot2::element_rect(fill = "white", colour = NA), # element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = "grey50"), # white background
      plot.background = ggplot2::element_rect(fill = "white", colour = "grey50"), # white background
      legend.key = ggplot2::element_blank(),
      text = ggplot2::element_text(size = fontsize,  family = fontfamily)
    ) # removed whitespace buffer around legend boxes which is nice

  ggplot2::ggsave(filename = filesavename,
                  plot = p, # ggplot2::last_plot()
                  device = "png",
                  path = savedir,
                  scale = 1,
                  #changes how big lines & legend items & axes & titles are relative to basemap. Smaller number = bigger items
                  width = 6, height = autoheight, units = "in", dpi = 600, limitsize = TRUE)

} # close function

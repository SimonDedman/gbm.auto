#' Data: Numbers of 4 adult female rays caught in 2137 Irish Sea trawls, 1994 to 2014
#'
#' 2137 capture events of adult female cuckoo, thornback, spotted and blonde
#' rays in the Irish Sea from 1994 to 2014 by the ICES IBTS, including
#' explanatory variables: Length Per Unit Effort in that area by the commercial
#' fishery, depth, temperature, distance to shore, and current speed at the
#' bottom. Note column 9 between Distance_to_shore and Current_Speed is blank.
#'
#' \itemize{
#'   \item Longitude. Decimal longitudes in the Irish Sea
#'   \item Latitude. Decimal latitudes in the Irish Sea
#'   \item Haul_Index. ICES IBTS area, survey, station, and year
#'   \item F_LPUE. Commercial fishery LPUE in Kg/Hr
#'   \item Depth. Metres, decimal.
#'   \item Temperature. Degrees, decimal.
#'   \item Salinity. PPM.
#'   \item Distance_to_Shore. Metres, decimal.
#'   \item Current_Speed. Metres per second at the seabed.
#'   \item Cuckoo. Numbers of cuckoo rays caught, standardised to 1 hour
#'   \item Thornback. Numbers of thornback rays caught, standardised to 1 hour
#'   \item Blonde. Numbers of blonde rays caught, standardised to 1 hour
#'   \item Spotted. Numbers of spotted rays caught, standardised to 1 hour
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Adult_Females
#' @usage data(Adult_Females)
#' @format A data frame with 2137 rows and 13 variables
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#' @references \url{http://datras.ices.dk}
"Adult_Females"

#' Data: Scaled abundance data for 2 subsets of 4 rays in the Irish Sea, by gbm.cons
#'
#' A dataset containing the output of the gbm.cons example run, conservation
#' priority areas within the Irish Sea for juvenile and adult female cuckoo,
#' blonde, thornback and spotted rays.
#'
#' \itemize{
#'   \item Longitude. Decimal longitudes in the Irish Sea
#'   \item Latitude. Decimal latitudes in the Irish Sea
#'   \item allscaled. Relative abundance. Each juvenile and adult female
#'   cuckoo, blonde, thornback and spotted ray scaled to 1 and added together.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name AllScaledData
#' @usage data(AllScaledData)
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#' @format A data frame with 378570 rows and 3 variables
"AllScaledData"

#' Data: Explanatory and response variables for 4 juvenile rays in the Irish Sea
#'
#' A dataset containing explanatory variables for environment, fishery and
#' predators of juvenile rays in the Irish Sea, and the response variables,
#' abundance CPUEs of cuckoo, thornback, blonde and spotted rays.
#'
#' \itemize{
#'   \item Survey_StNo_HaulNo_Year.
#'   \item Latitude. Decimal latitudes in the Irish Sea
#'   \item Longitude. Decimal longitudes in the Irish Sea
#'   \item Depth. Metres, decimal.
#'   \item Temperature. Degrees, decimal.
#'   \item Salinity. PPM.
#'   \item Current_Speed. Metres per second at the seabed.
#'   \item Distance_to_Shore. Metres, decimal.
#'   \item F_LPUE. Commercial fishery LPUE in Kg/Hr
#'   \item Scallop. Average KwH Scallop effort from logbooks, Marine Institute and MMO combined
#'   \item MI_Av_E_Hr. Average effort hours, Marine Institute Scallop VMS,
#'   0.03*0.02 rectangles, all Irish Sea, 2006-14
#'   \item MI_Av_LPUE. Average scallop CPUE, Marine Institute Scallop VMS,
#'   0.03*0.02 rectangles, all Irish Sea, 2006-14
#'   \item MI_Sum_Liv. Sum of live weight. Average scallop CPUE, Marine Institute Scallop VMS,
#'   0.03*0.02 rectangles, all Irish Sea, 2006-14
#'   \item Whelk. MMO Whelk LPUE 2009-12, pivot, polygons to points
#'   \item MmoAvScKwh. MMO Scallop Effort 2009-12, pivot, polygons to points. ICES rectangles.
#'   \item Cod_C. ICES IBTS CPUE of cod caught between 1994 - 2014 large enough to predate upon <= year 1 cuckoo rays.
#'   \item Cod_T. As Cod_C for yr1 thornback rays.
#'   \item Cod_B. As Cod_C for yr1 blonde rays.
#'   \item Cod_S. As Cod_C for yr1 spotted rays.
#'   \item Haddock_C. As Cod_C, haddock predating upon cuckoo rays
#'   \item Haddock_T. As Cod_C, haddock predating upon thornback rays
#'   \item Haddock_B. As Cod_C, haddock predating upon blonde rays
#'   \item Haddock_S. As Cod_C, haddock predating upon spotted rays
#'   \item Plaice_C. As Cod_C, plaice predating upon cuckoo rays
#'   \item Plaice_T. As Cod_C, plaice predating upon thornback rays
#'   \item Plaice_B. As Cod_C, plaice predating upon blonde rays
#'   \item Plaice_S. As Cod_C, plaice predating upon spotted rays
#'   \item Whiting_C. As Cod_C, whiting predating upon cuckoo rays
#'   \item Whiting_T. As Cod_C, whiting predating upon thornback rays
#'   \item Whiting_B. As Cod_C, whiting predating upon blonde rays
#'   \item Whiting_S. As Cod_C, whiting predating upon spotted rays
#'   \item ComSkt_C. As Cod_C, common skate predating upon cuckoo rays
#'   \item ComSkt_T. As Cod_C, common skate predating upon thornback rays
#'   \item ComSkt_B. As Cod_C, common skate predating upon blonde rays
#'   \item ComSkt_S. As Cod_C, common skate predating upon spotted rays
#'   \item Blonde_C. As Cod_C, blonde ray predating upon cuckoo rays
#'   \item Blonde_T. As Cod_C, blonde ray predating upon thornback rays
#'   \item Blonde_S. As Cod_C, blonde ray predating upon spotted rays
#'   \item C_Preds. All predator CPUEs combined for cuckoo rays.
#'   \item T_Preds. All predator CPUEs combined for thornback rays.
#'   \item B_Preds. All predator CPUEs combined for blonde rays.
#'   \item S_Preds. All predator CPUEs combined for spotted rays.
#'   \item Cuckoo. Numbers of juvenile cuckoo rays caught, standardised to 1 hour
#'   \item Thornback. Numbers of juvenile thornback rays caught, standardised to 1 hour
#'   \item Blonde. Numbers of juvenile blonde rays caught, standardised to 1 hour
#'   \item Spotted. Numbers of juvenile spotted rays caught, standardised to 1 hour
#' }
#'
#' @docType data
#' @keywords datasets
#' @name Juveniles
#' @usage data(Juveniles)
#' @format A data frame with 2136 rows and 46 variables
#' @author Simon Dedman, \email{simondedman@@gmail.com}
"Juveniles"

#' Data: Explanatory variables for rays in the Irish Sea
#'
#' A dataset containing explanatory variables for environment, fishery and
#' predators of rays including juveniles in the Irish Sea.
#'
#' \itemize{
#'   \item Longitude. Decimal longitudes in the Irish Sea
#'   \item Latitude. Decimal latitudes in the Irish Sea
#'   \item Depth. Metres, decimal.
#'   \item Temperature. Degrees, decimal.
#'   \item Salinity. PPM.
#'   \item Current_Speed. Metres per second at the seabed.
#'   \item Distance_to_Shore. Metres, decimal.
#'   \item F_LPUE. Commercial fishery LPUE in Kg/Hr
#'   \item Scallop. Average KwH Scallop effort from logbooks, Marine Institute and MMO combined
#'   \item MI_Av_E_Hr. Average effort hours, Marine Institute Scallop VMS,
#'   0.03*0.02 rectangles, all Irish Sea, 2006-14
#'   \item MI_Av_LPUE. Average scallop CPUE, Marine Institute Scallop VMS,
#'   0.03*0.02 rectangles, all Irish Sea, 2006-14
#'   \item MI_Sum_Liv. Sum of live weight. Average scallop CPUE, Marine Institute Scallop VMS,
#'   0.03*0.02 rectangles, all Irish Sea, 2006-14
#'   \item Whelk. MMO Whelk LPUE 2009-12, pivot, polygons to points
#'   \item MmoAvScKwh. MMO Scallop Effort 2009-12, pivot, polygons to points. ICES rectangles.
#'   \item HubDist. map calc, distance of grid point to nearest datras point representing it (for preds)
#'   \item Cod_C. ICES IBTS CPUE of cod caught between 1994 - 2014 large enough to predate upon <= year 1 cuckoo rays.
#'   \item Cod_T. As Cod_C for yr1 thornback rays.
#'   \item Cod_B. As Cod_C for yr1 blonde rays.
#'   \item Cod_S. As Cod_C for yr1 spotted rays.
#'   \item Haddock_C. As Cod_C, haddock predating upon cuckoo rays
#'   \item Haddock_T. As Cod_C, haddock predating upon thornback rays
#'   \item Haddock_B. As Cod_C, haddock predating upon blonde rays
#'   \item Haddock_S. As Cod_C, haddock predating upon spotted rays
#'   \item Plaice_C. As Cod_C, plaice predating upon cuckoo rays
#'   \item Plaice_T. As Cod_C, plaice predating upon thornback rays
#'   \item Plaice_B. As Cod_C, plaice predating upon blonde rays
#'   \item Plaice_S. As Cod_C, plaice predating upon spotted rays
#'   \item Whiting_C. As Cod_C, whiting predating upon cuckoo rays
#'   \item Whiting_T. As Cod_C, whiting predating upon thornback rays
#'   \item Whiting_B. As Cod_C, whiting predating upon blonde rays
#'   \item Whiting_S. As Cod_C, whiting predating upon spotted rays
#'   \item ComSkt_C. As Cod_C, common skate predating upon cuckoo rays
#'   \item ComSkt_T. As Cod_C, common skate predating upon thornback rays
#'   \item ComSkt_B. As Cod_C, common skate predating upon blonde rays
#'   \item ComSkt_S. As Cod_C, common skate predating upon spotted rays
#'   \item Blonde_C. As Cod_C, blonde ray predating upon cuckoo rays
#'   \item Blonde_T. As Cod_C, blonde ray predating upon thornback rays
#'   \item Blonde_S. As Cod_C, blonde ray predating upon spotted rays
#'   \item C_Preds. All predator CPUEs combined for cuckoo rays.
#'   \item T_Preds. All predator CPUEs combined for thornback rays.
#'   \item B_Preds. All predator CPUEs combined for blonde rays.
#'   \item S_Preds. All predator CPUEs combined for spotted rays.
#'   \item Effort.  Irish commercial beam trawler effort 2012
#' }
#'
#' @docType data
#' @keywords datasets
#' @name grids
#' @usage data(grids)
#' @format A data frame with 378570 rows and 43 variables
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#' @references \url{http://oar.marine.ie/handle/10793/958}
"grids"

#' Data: Numbers of 4 ray species caught in 2137 Irish Sea trawls, 1994 to 2014
#'
#' 2244 capture events of cuckoo, thornback, spotted and blonde rays in the Irish
#' Sea from 1994 to 2014 by the ICES IBTS, including explanatory variables:
#' Length Per Unit Effort in that area by the commercial fishery, fishing effort
#' by same, depth, temperature, distance to shore, and current speed at the
#' bottom. Note column 8 between Distance_to_shore and Current_Speed is blank.
#'
#' \itemize{
#'   \item Survey_StNo_HaulNo_Year.
#'   \item Latitude. Decimal latitudes in the Irish Sea
#'   \item Longitude. Decimal longitudes in the Irish Sea
#'   \item Depth. Metres, decimal.
#'   \item Temperature. Degrees, decimal.
#'   \item Salinity. PPM.
#'   \item Current_Speed. Metres per second at the seabed.
#'   \item Distance_to_Shore. Metres, decimal.
#'   \item F_LPUE. Commercial fishery LPUE in Kg/Hr
#'   \item Effort. Irish commercial beam trawler effort 2012
#'   \item Cuckoo. Numbers of juvenile cuckoo rays caught, standardised to 1 hour
#'   \item Thornback. Numbers of juvenile thornback rays caught, standardised to 1 hour
#'   \item Blonde. Numbers of juvenile blonde rays caught, standardised to 1 hour
#'   \item Spotted. Numbers of juvenile spotted rays caught, standardised to 1 hour
#' }
#'
#' @docType data
#' @keywords datasets
#' @name samples
#' @usage data(samples)
#' @format A data frame with 2244 rows and 14 variables
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#' @references \url{http://oar.marine.ie/handle/10793/958}
"samples"

#' Data: Predicted abundances of 4 ray species generated using gbm.auto
#'
#' Predicted abundances of 4 ray species generated using gbm.auto, and
#' Irish commercial beam trawler effort 2012.
#'
#' \itemize{
#'   \item Latitude. Decimal latitudes in the Irish Sea
#'   \item Longitude. Decimal longitudes in the Irish Sea
#'   \item Cuckoo. Predicted abundances of cuckoo rays in the Irish Sea, generated using gbm.auto
#'   \item Thornback. Predicted abundances of thornback rays in the Irish Sea, generated using gbm.auto
#'   \item Blonde. Predicted abundances of blonde rays in the Irish Sea, generated using gbm.auto
#'   \item Spotted. Predicted abundances of spotted rays in the Irish Sea, generated using gbm.auto
#'   \item Effort. Irish commercial beam trawler effort 2012
#' }
#'
#' @docType data
#' @keywords datasets
#' @name AllPreds_E
#' @usage data(AllPreds_E)
#' @format A data frame with 378570 rows and 7 variables
#' @author Simon Dedman, \email{simondedman@@gmail.com}
"AllPreds_E"

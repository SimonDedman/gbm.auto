% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{MySamples}
\alias{MySamples}
\title{Data: Numbers of 4 ray species caught in 2137 Irish Sea trawls, 1994 to 2014}
\format{
A data frame with 2244 rows and 14 variables:
\describe{
\item{Survey_StNo_HaulNo_Year}{Index column of combined Survey number, station number, haul number, and year}
\item{Latitude}{Decimal latitudes in the Irish Sea}
\item{Longitude}{Decimal longitudes in the Irish Sea}
\item{Depth}{ Metres, decimal}
\item{Temperature}{Degrees, decimal}
\item{Salinity}{PPM}
\item{Current_Speed}{Metres per second at the seabed}
\item{Distance_to_Shore}{Metres, decimal}
\item{F_LPUE}{Commercial fishery LPUE in Kg/Hr}
\item{Effort}{Irish commercial beam trawler effort 2012}
\item{Cuckoo}{Numbers of juvenile cuckoo rays caught, standardised to 1 hour}
\item{Thornback}{Numbers of juvenile thornback rays caught, standardised to 1 hour}
\item{Blonde}{Numbers of juvenile blonde rays caught, standardised to 1 hour}
\item{Spotted}{Numbers of juvenile spotted rays caught, standardised to 1 hour}
}
}
\source{
\url{http://oar.marine.ie/handle/10793/958}
}
\usage{
data(MySamples)
}
\description{
2244 capture events of cuckoo, thornback, spotted and blonde rays in the Irish
Sea from 1994 to 2014 by the ICES IBTS, including explanatory variables:
Length Per Unit Effort in that area by the commercial fishery, fishing effort
by same, depth, temperature, distance to shore, and current speed at the
bottom.
}
\author{
Simon Dedman, \email{simondedman@gmail.com}
}
\keyword{datasets}

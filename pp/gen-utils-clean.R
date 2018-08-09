library(pracma)
library('SDMTools')


#' Obtains the coordinates from string values stored in regions for JSONs
#' @param coordsStr "latitude,longitude"
#' @return data.frame(lat, lon)
getCoords <- function(coordsStr) {
  coords <- strsplit(coordsStr, ",")[[1]]
  return(data.frame(lon = as.numeric(coords[2]), lat = as.numeric(coords[1])))
}

#' It generates an intensity function that is the sum of more intensity
#' functions centered at certain locations.
#' Note: the distance used is geodetic WGS84 using Vicenty's formula
#'       hence centers are WGS84 coordinates (latitude, longitude)
#'
#' @param centers dataframe with lon and lat coordinates of centers
#' @param centers_intensity probability density function to place at each center
#' @param ... arguments for the probability density function
#'
#' @return intensity function that is the sum of the centered ones.
#'         function(latitude, longitude)
gen_manta <- function(centers, centers_intensity, ...) {
  # Defines the intensity functions for each center
  center_ints <- list()
  for (i in 1:dim(centers)[1]) {
    center_ints[[i]] <- local({
      origin <- centers[i,]
      
      function(longi, lati) {
        # Vicenty's distance
        r <- SDMTools::distance(lat1=lati, lon1=longi,
                 lat2=origin$lat, lon2=origin$lon)$distance
         
        return(centers_intensity(r, ...))
      }
      
    })
  }
  
  # Creates the function that is sum of intensities
  sum_int <- function(longi, lati) {
    sum_prob <- 0
    for (center_int in center_ints) {
      sum_prob <- sum_prob + center_int(longi, lati)
    }
    return(sum_prob)
  }
  
  return(sum_int)
}


#' Generates a dataframe evaluating an intensity measure function in a grid of
#' coordinates.
#'
#' @param latis vector (latitude_min, latitude_max)
#' @param longis vector (longitude_min, longitude_max)
#' @param latisLen sampling length in latitude axis
#' @param longisLen sampling length in longitude axis
#' @param intF intensity measure function used to generate the values
#'
#' @return a datafram having latitude, longitude and intensity
genIntFrame <- function(latis, longis, latisLen, longisLen, intF) {
  latiAxes <- seq(from = latis[1], to = latis[2], length = latisLen)
  longiAxes <- seq(from = longis[1], to = longis[2], length = longisLen)
  
  intFrame <- expand.grid(lon = longiAxes, lat = latiAxes)
  intFrame <- data.frame(lon = intFrame$lon, lat = intFrame$lat,
                          intensity = array(dim = nrow(intFrame)))
  for (row in 1:nrow(intFrame)) {
    intFrame$intensity[row] <- intF(intFrame$lon[row], intFrame$lat[row])
  }
  
  return(intFrame)
}


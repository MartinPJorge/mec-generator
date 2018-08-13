library(cluster)
library('SDMTools') # for Vicenty's distance

#' @description Obtains the medoids for a grid within which antennas fall in.
#' @param anntenasLoc matrix with two columns -latitude,longitude-
#' @param lonSteps array of longitudes for the grid
#' @param latSteps array of latitudes for the grid
#' @return a dataframe with lat and lon of the medoid antennas
getMedoids <- function(antennasLoc, lonSteps, latSteps) {
  medoidsLat <- c()
  medoidsLon <- c()
  
  # Iterate through each square
  for (lon in 2:length(lonSteps)) {
    for (lat in 2:length(latSteps)) {
      lonA <- lonSteps[lon - 1]
      lonB <- lonSteps[lon]
      latA <- latSteps[lat - 1]
      latB <- latSteps[lat]
      
      # Find first antenna fitting in
      for (antenna in nrow(antennasLoc)) {
        antenaLat <- antennasLoc$lat[antenna]
        antenaLon <- antennasLoc$lon[antenna]
        insideLat <- latA <= antennaLat & antennaLat <= latB
        insideLon <- lonA <= antennaLon & antennaLon <= lonB
        
        if (insideLat & insideLon) {
          medoidsLat <- c(medoidsLat, antennaLat)
          medoidsLon <- c(medoidsLon, antennaLon)
          break
        }
      }
    }
  }
  
  return(data.frame(lat = medoidsLat, lon = medoidsLon))
}


#' @description Performs the PAM clustering based on a given grid.
#' @param anntenasLoc matrix with two columns -latitude,longitude-
#' @param lonSteps array of longitudes for the grid
#' @param latSteps array of latitudes for the grid
#' @return TODO
pam_antennas <- function(antennasLoc, lonSteps, latSteps) {
  
}

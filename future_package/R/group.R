#############################################################
############# START FROM HERE ON THE ADAPTATION #############
#############################################################


library(fields)
library(geosphere)

#' @description It groups all elements with coordinates (lats, lons) within
#' groups of a maximum of groupN elements that are <maxDis from the group
#' center. It uses the agglomerative clustering with the complete method
#' @param lats elements' latitudes
#' @param lons elements' longitudes
#' @param groupN number of elements per group
#' @param maxDis maximum distance, in meters, between elements in the same group
#' @return a data frame with the group assignments for each element coordinate
#' data.frame(lat, lon, group)
groupElems <- function(lats, lons, groupN, maxDis) {
  
  # Complete agglom clustering to ensure intra distances are <maxDis
  coords = data.frame(lats, lons)
  dist <- rdist.earth(coords, miles = FALSE) 
  fit <- hclust(as.dist(dist), method = "complete")
  clusters <- cutree(fit, h = maxDis / 1000) # our distance matrix is in km
  
  groupLons <- c()
  groupLats <- c()
  groupIds <- c()
  currGroup <- 0
  
  # Create the groups using the clusters of elements <maxDis
  for (cluster in clusters) {
    elemsIdx <- which(clusters == cluster)
    while (length(elemsIdx) > 0) {
      poppedIdx <- head(elemsIdx, n = groupN)
      groupLons <- c(groupLons, lons[poppedIdx])
      groupLats <- c(groupLats, lats[poppedIdx])
      groupIds <- c(groupIds, rep(x = currGroup, times = length(poppedIdx)))
      currGroup <- currGroup + 1
      
      tailElems <- max(0, length(elemsIdx) - length(poppedIdx))
      elemsIdx <- tail(elemsIdx, n = tailElems)
    }
  }
  
  return(data.frame(lat = groupLats, lon = groupLons, group = groupIds))
}


#' @description It creates centers for grouped elements
#' @param lats elements' latitudes
#' @param lons elements' longitudes
#' @param groups elements' assigned group
#' @return the coordinates of each group center data.frame(lat, lon, group)
groupCenters <- function(lats, lons, groups) {
  centerLats <- c()
  centerLons <- c()
  centerGroups <- c()
  
  for (group in unique(groups)) {
    groupLats <- lats[which(groups == group)]
    groupLons <- lons[which(groups == group)]
    groupCenter <- centroid(matrix(c(groupLats, groupLons), ncol = 2,
                                   nrows = length(groupLats)))
    
    centerLats <- c(centerLats, groupCenter[1])
    centerLons <- c(centerLons, groupCenter[])
    centerGroups <- c(centerGroups, group)
  }
  
  return(data.frame(lat = centerLats, lon = centerLons, group = centerGroups))
}


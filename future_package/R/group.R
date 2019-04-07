#############################################################
############# START FROM HERE ON THE ADAPTATION #############
#############################################################


library(fields)
library(geosphere)
library(igraph)

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
  for (cluster in unique(clusters)) {
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
    if (length(groupLons) == 1) {
      centerLats <- c(centerLats, groupLats[1])
      centerLons <- c(centerLons, groupLons[1])
    } else if (length(groupLons) == 2) {
      centerLats <- c(centerLats, (groupLats[1] + groupLats[2]) / 2)
      centerLons <- c(centerLons, (groupLons[1] + groupLons[2]) / 2)
    } else {
      groupCenter <- centroid(matrix(c(groupLats, groupLons), ncol = 2,
                                     nrow = length(groupLats)))
      
      centerLats <- c(centerLats, groupCenter[1])
      centerLons <- c(centerLons, groupCenter[2])
    }
    centerGroups <- c(centerGroups, group)
  }
  
  return(data.frame(lat = centerLats, lon = centerLons, group = centerGroups))
}


#' @description Writes the 5G scenario as a graph. Each data.frame parameter
#' contains lat and lon coordinates and group columns to indicate associations 
#' among elements.
#' @param m1Assoc data.frame cell coordinates and the M1 ids
#' @param m1Coords data.frame M1 coordinates and their ids
#' @param m1AccAssocs data.frame M1 coordinates and the access ring ids
#' @param accCentCoords data.frame access ring center coordinates and their ids
#' @param m2Assocs data.frame access ring coords and associated M2 ids
#' @param m2Switches data.frame M2 switches coordinates and their ids
#' @param m2AggAssocs data.frame M2 coordinates and aggregation rings ids
#' @param aggCentCoords data.frame aggregation ring centers and their ids
#' @param m3Assocs data.frame aggregation ring centers and M3 ids
#' @param m3Switches data.frame with M3 coordinates and their ids
#' @param format output graph format
#' @param file path to the file where the graph is stored
write5GtoGraph <- function(m1Assoc, m1Coords, m1AccAssocs, accCentCoords,
                           m2Assocs, m2Switches, m2AggAssocs, aggCentCoords,
                           m3Assocs, m3Switches, format, file) {
  
  g <- make_empty_graph(directed = FALSE)
  
  ################## NODES CREATION ##################
  nodeIds <- c()
  nodeTypes <- c()
  nodeLons <- c()
  nodeLats <- c()
  
  # Create the cell nodes
  for (row in nrow(m1Assoc)) {
    nodeIds <- c(nodeIds, cat("cell", row, sep = ""))
    nodeTypes <- c(nodeTypes, "cell")
    nodeLons <- c(nodeLons, m1Assoc[row,]$lon)
    nodeLats <- c(nodeLats, m1Assoc[row,]$lat)
  }
  
  # Create the M1 nodes
  for (row in nrow(m1Coords)) {
    nodeIds <- c(nodeIds, cat("m1_", m1Coords[row,]$group, sep = ""))
    nodeTypes <- c(nodeTypes, "m1")
    nodeLons <- c(nodeLons, m1Coords[row,]$lon)
    nodeLats <- c(nodeLats, m1Coords[row,]$lat)
  }
  
  # Create the M2 nodes
  for (row in nrow(m2Switches)) {
    nodeIds <- c(nodeIds, cat("m2_", m2Switches[row,]$group, sep = ""))
    nodeTypes <- c(nodeTypes, "m2")
    nodeLons <- c(nodeLons, m2Switches[row,]$lon)
    nodeLats <- c(nodeLats, m2Switches[row,]$lat)
  }
  
  # Create the M3 nodes with the redundancy
  for (row in nrow(m3Switches)) {
    nodeIds <- c(nodeIds, cat("m3_", m3Switches[row,]$group, sep = ""))
    nodeTypes <- c(nodeTypes, "m3")
    nodeLons <- c(nodeLons, m3Switches[row,]$lon)
    nodeLats <- c(nodeLats, m3Switches[row,]$lat)
    
    nodeIds <- c(nodeIds, cat("m3_rep_", m3Switches[row,]$group, sep = ""))
    nodeTypes <- c(nodeTypes, "m3")
    nodeLons <- c(nodeLons, m3Switches[row,]$lon)
    nodeLats <- c(nodeLats, m3Switches[row,]$lat)
  }
  nodes <- cbind(id = nodeIds, type = nodeTypes, lon = nodeLons, lat = nodeLats)
  
  
  ################## LINKS CREATION ##################
  origins <- c()
  ends <- c()
  bandwidths <- c()
  bandwidthUnits <- c()
  distances <- c()
  distanceUnits <- c()
  
  # Link cells to M1 switches
  for (row in nrow(m1Assoc)) {
    origins <- c(origins, cat("cell", row, sep = ""))
    ends <- c(ends, cat("m1_", m1Assoc[row,]$group, sep = ""))
    bandwidths <- c(bandwidths, "10")
    bandwidthUnits <- c(bandWidthUnits, "Gb/s")
    
    m1Coords <- m1Coords[which(m1Coords$group == m1Assoc$group),][1,]
    dis <- SDMTools::distance(lat1 = m1Assoc[row,]$lat,
                              lon1 = m1Assoc[row,]$lon,
                              lat2 = m1Coords$lat, lon2 = m1Coords$lon)$distance
    distances <- c(distances, dis)
    distanceUnits <- c(distances, "meters")
  }
  
  # TODO - add remaining links
}


#' @description Builds the 5G scenario of Figure 1 in:
#' https://ieeexplore.ieee.org/mediastore_new/IEEE/content/media/8412472/8436592/8436847/paper_23-fig-1-source-large.gif
#' @param lats cell antenna latitudes
#' @param lons cell antenna longitudes
build5GScenario <- function(lats, lons) {
  # Create the M1 associations
  m1Assoc <- groupElems(lats = lats, lons = lons, groupN = 6,
                             maxDis = 10000)
  m1Coords <- groupCenters(lats = m1Assoc$lat, lons = m1Assoc$lon,
                           groups = m1Assoc$group)
  
  # Create the access rings with their center coordinates
  m1AccAssocs <- groupElems(lats = m1Coords$lat, lons = m1Coords$lon,
                            groupN = 6, maxDis = 20000)
  accCentCoords <- groupCenters(lats = m1AccAssocs$lat, lons = m1AccAssocs$lon,
                                groups = m1AccAssocs$group)
  
  # Create the M2 switches grouping by 4 the access rings centers
  m2Assocs <- groupElems(lats = accCentCoords$lat, lons = accCentCoords$lon,
                           groupN = 4, maxDis = 20000)
  m2Switches <- groupCenters(lats = m2Assocs$lat, lons = m2Assocs$lon,
                             groups = m2Assocs$group)
  
  # Create the aggregation rings with their center coordinates
  m2AggAssocs <- groupElems(lats = m2Switches$lat, lons = m2Switches$lon,
                            groupN = 6, maxDis = 40000)
  aggCentCoords <- groupCenters(lats = m2AggAssocs$lat, lons = m2AggAssocs$lon,
                                groups = m2AggAssocs$group)
  
  # Create the M3 switches grouping by 2 the aggregation rings
  m3Assocs <- groupElems(lats = aggCentCoords$lat, lons = aggCentCoords$lon,
                         groupN = 2, maxDis = 80000)
  m3Switches <- groupCenters(lats = m3Assocs$lat, lons = m3Assocs$lon,
                             groups = m3Assocs$group)
}



g <- make_empty_graph(n = 2, directed = FALSE)
add_edges(graph = g,  edges = c(1,2))



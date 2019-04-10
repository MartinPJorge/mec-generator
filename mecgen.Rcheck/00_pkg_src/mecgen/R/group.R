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
  if (length(lats) == 1) {
    return(data.frame(lat = lats[1], lon = lons[1], group = 0))
  }

  # Complete agglom clustering to ensure intra distances are <maxDis
  coords = data.frame(lats, lons)
  dist <- fields::rdist.earth(coords, miles = FALSE)
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
      groupCenter <- geosphere::centroid(matrix(c(groupLons, groupLats), ncol = 2,
                                     nrow = length(groupLats)))

      # For multipart polygons, just average coordinates
      if (groupCenter[2] < min(groupLats) | groupCenter[2] > max(groupLats) |
          groupCenter[1] < min(groupLons) | groupCenter[1] > max(groupLons)) {
        groupCenter[1] <- mean(groupLons)
        groupCenter[2] <- mean(groupLats)
      }

      centerLons <- c(centerLons, groupCenter[1])
      centerLats <- c(centerLats, groupCenter[2])
    }
    centerGroups <- c(centerGroups, group)
  }

  return(data.frame(lat = centerLats, lon = centerLons, group = centerGroups))
}


#' @description Create the access rings' edges list
#' @param m1Assoc data.frame cell coordinates and the M1 ids
#' @param m1Coords data.frame M1 coordinates and their ids
#' @param m1AccAssocs data.frame M1 coordinates and the access ring ids
#' @param accCentCoords data.frame access ring center coordinates and their ids
#' @param m2Assocs data.frame access ring coords and associated M2 ids
#' @return list(from, to, bandwidth, bandwidthUnits, distance,
#' distanceUnits)
linkAccessRings <- function(m1Assoc, m1Coords, m1AccAssocs, accCentCoords,
                            m2Assocs) {
  origins <- c()
  ends <- c()
  bandwidths <- c()
  bandwidthUnits <- c()
  distances <- c()
  distanceUnits <- c()

  for (accRing in unique(m1AccAssocs$group)) {
    # Get the access ring M1 ids
    assocM1s <- m1AccAssocs[which(m1AccAssocs$group == accRing),]
    ringM1s <- NULL

    ringM1s <- m1Coords[which(m1Coords$lon == assocM1s$lon),]
    ringM1s <- ringM1s[which(ringM1s$lat == assocM1s$lat),]

    # Get the M2 for the access ring
    accRingCenter <- accCentCoords[which(accCentCoords$group == accRing),]
    m2 <- m2Assocs[which(m2Assocs$lon == accRingCenter$lon),]
    m2 <- m2[which(m2$lat == accRingCenter$lat),]

    # Create the access ring in the graph
    if (nrow(ringM1s) > 1) {
      for (m1Row in 2:nrow(ringM1s)) {
        m1Src <- ringM1s[m1Row-1,]
        m1Dst <- ringM1s[m1Row,]
        origins <- c(origins, paste("m1_", m1Src$group, sep = ""))
        ends <- c(ends, paste("m1_", m1Dst$group, sep = ""))
        bandwidths <- c(bandwidths, "300")
        bandwidthUnits <- c(bandwidthUnits, "Gb/s")
        dis <- SDMTools::distance(lat1 = m1Src$lat, lon1 = m1Src$lon,
                                  lat2 = m1Dst$lat, lon2 = m1Dst$lon)$distance
        distances <- c(distances, dis)
        distanceUnits <- c(distanceUnits, "meters")
      }
    }

    # Link last M1 to M2
    lastM1 <- tail(ringM1s, n = 1)
    origins <- c(origins, paste("m1_", lastM1$group, sep = ""))
    ends <- c(ends, paste("m2_", m2$group, sep = ""))
    bandwidths <- c(bandwidths, "300")
    bandwidthUnits <- c(bandwidthUnits, "Gb/s")
    dis <- SDMTools::distance(lat1 = lastM1$lat, lon1 = lastM1$lon,
                              lat2 = m2$lat, lon2 = m2$lon)$distance
    distances <- c(distances, dis)
    distanceUnits <- c(distanceUnits, "meters")

    # Link M2 to first M1
    firstM1 <- head(ringM1s, n = 1)
    origins <- c(origins, paste("m2_", m2$group, sep = ""))
    ends <- c(ends, paste("m1_", firstM1$group, sep = ""))
    bandwidths <- c(bandwidths, 300)
    bandwidthUnits <- c(bandwidthUnits, "Gb/s")
    dis <- SDMTools::distance(lat1 = m2$lat, lon1 = m2$lon,
                              lat2 = firstM1$lat, lon2 = firstM1$lon)$distance
    distances <- c(distances, dis)
    distanceUnits <- c(distanceUnits, "meters")
  }


  return(list(from = origins, to = ends, bandwidth = bandwidths,
                    bandwidthUnits = bandwidthUnits, distance = distances,
                    distanceUnits = distanceUnits))
}


#' @description Create the aggregation ring edges list
#' @param m2Assocs data.frame access ring coords and associated M2 ids
#' @param m2Switches data.frame M2 switches coordinates and their ids
#' @param m2AggAssocs data.frame M2 coordinates and aggregation rings ids
#' @param aggCentCoords data.frame aggregation ring centers and their ids
#' @param m3Assocs data.frame aggregation ring centers and M3 ids
#' @param m3Switches data.frame with M3 coordinates and their ids
#' @return list(from, to, bandwidth, bandwidthUnits, distance,
#' distanceUnits)
linkAggRings <- function(m2Assocs, m2Switches, m2AggAssocs, aggCentCoords,
                         m3Assocs, m3Switches) {
  links <- NULL

  origins <- c()
  ends <- c()
  bandwidths <- c()
  bandwidthUnits <- c()
  distances <- c()
  distanceUnits <- c()

  # Iterate through each aggregation
  for (aggRing in unique(m2AggAssocs$group)) {
    # Find the M2 switches of the aggregation ring
    m2Coords <- m2AggAssocs[which(m2AggAssocs$group == aggRing),]
    m2s <- m2Switches[which(m2Switches$lon == m2Coords$lon),]
    m2s <- m2s[which(m2s$lat == m2Coords$lat),]

    # Find the M3 switch of the aggegation ring
    aggCenter <- aggCentCoords[which(aggCentCoords$group == aggRing),]
    m3Id <- m3Assocs[which(m3Assocs$lon == aggCenter$lon &&
                           m3Assocs$lat == aggCenter$lat),]$group
    m3 <- m3Switches[which(m3Switches$group == m3Id),]

    # Create M2 ring links
    if (nrow(m2s) > 1) {
      for (row in 2:nrow(m2s)) {
        dis <- SDMTools::distance(lat1 = m2s[row-1,]$lat,
                                  lon1 = m2s[row-1,]$lon,
                                  lat2 = m2s[row,]$lat,
                                  lon2 = m2s[row,]$lon)$distance
        origins <- c(origins, paste("m2_", m2s[row-1,]$group, sep = ""))
        ends <- c(ends, paste("m2_", m2s[row,]$group, sep = ""))
        bandwidths <- c(bandwidths, 6)
        bandwidthUnits <- c(bandwidthUnits, "Tb/s")
        distances <- c(distances, dis)
        distanceUnits <- c(distanceUnits, "meters")
      }
    }

    # Link last M2 to M3
    lastM2 <- tail(m2s, n = 1)
    dis <- SDMTools::distance(lat1 = lastM2$lat, lon1 = lastM2$lon,
                              lat2 = m3$lat, lon2 = m3$lon)$distance

    origins <- c(origins, paste("m2_", lastM2$group, sep = ""))
    ends <- c(ends, paste("m3_", m3$group, sep = ""))
    bandwidths <- c(bandwidths, 6)
    bandwidthUnits <- c(bandwidthUnits, "Tb/s")
    distances <- c(distances, dis)
    distanceUnits <- c(distanceUnits, "meters")

    # Link M3 to redundancy M3
    origins <- c(origins, paste("m3_", m3$group, sep = ""))
    ends <- c(ends, paste("m3_rep_", m3$group, sep = ""))
    bandwidths <- c(bandwidths, 6)
    bandwidthUnits <- c(bandwidthUnits, "Tb/s")
    distances <- c(distances, dis)
    distanceUnits <- c(distanceUnits, "meters")


    # Link redundancy M3 to first M2 of the ring
    firstM2 <- head(m2s, n = 1)
    dis <- SDMTools::distance(lat1 = firstM2$lat, lon1 = firstM2$lon,
                              lat2 = m3$lat, lon2 = m3$lon)$distance
    origins <- c(origins, paste("m3_rep_", m3$group, sep = ""))
    ends <- c(ends, paste("m2_", firstM2$group, sep = ""))
    bandwidths <- c(bandwidths, 6)
    bandwidthUnits <- c(bandwidthUnits, "Tb/s")
    distances <- c(distances, dis)
    distanceUnits <- c(distanceUnits, "meters")
  }

  return(list(from = origins, to = ends, bandwidth = bandwidths,
                    bandwidthUnits = bandwidthUnits, distance = distances,
                    distanceUnits = distanceUnits))
}


#' @export
#' @name graphFrames
#' @title graphFrames
#' @description creates the nodes and links data frames for the graph creation
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
#' @param servers data.frame with created servers
#' @param fogNodes data.frame with fog computing nodes
#' @param endpoints data.frame with endpoints
#' @return list(links=data.frame, nodes=data.frame)
graphFrames <- function(m1Assoc, m1Coords, m1AccAssocs, accCentCoords,
                           m2Assocs, m2Switches, m2AggAssocs, aggCentCoords,
                           m3Assocs, m3Switches, servers = NULL,
                        fogNodes = NULL, endpoints = NULL) {

  ################## NODES CREATION ##################
  nodes <- NULL
  nodeIds <- c()
  nodeTypes <- c()
  nodeLons <- c()
  nodeLats <- c()

  # Create the cell nodes
  for (row in 1:nrow(m1Assoc)) {
    nodeIds <- c(nodeIds, paste("cell", row, sep = ""))
    nodeTypes <- c(nodeTypes, "cell")
    nodeLons <- c(nodeLons, m1Assoc[row,]$lon)
    nodeLats <- c(nodeLats, m1Assoc[row,]$lat)
  }

  # Create the M1 nodes
  for (row in 1:nrow(m1Coords)) {
    nodeIds <- c(nodeIds, paste("m1_", m1Coords[row,]$group, sep = ""))
    nodeTypes <- c(nodeTypes, "m1")
    nodeLons <- c(nodeLons, m1Coords[row,]$lon)
    nodeLats <- c(nodeLats, m1Coords[row,]$lat)
  }

  # Create the M2 nodes
  for (row in 1:nrow(m2Switches)) {
    nodeIds <- c(nodeIds, paste("m2_", m2Switches[row,]$group, sep = ""))
    nodeTypes <- c(nodeTypes, "m2")
    nodeLons <- c(nodeLons, m2Switches[row,]$lon)
    nodeLats <- c(nodeLats, m2Switches[row,]$lat)
  }

  # Create the M3 nodes with the redundancy
  for (row in 1:nrow(m3Switches)) {
    nodeIds <- c(nodeIds, paste("m3_", m3Switches[row,]$group, sep = ""))
    nodeTypes <- c(nodeTypes, "m3")
    nodeLons <- c(nodeLons, m3Switches[row,]$lon)
    nodeLats <- c(nodeLats, m3Switches[row,]$lat)

    nodeIds <- c(nodeIds, paste("m3_rep_", m3Switches[row,]$group, sep = ""))
    nodeTypes <- c(nodeTypes, "m3")
    nodeLons <- c(nodeLons, m3Switches[row,]$lon)
    nodeLats <- c(nodeLats, m3Switches[row,]$lat)
  }

  nodes <- data.frame(id = nodeIds, type = nodeTypes, lon = nodeLons,
                      lat = nodeLats)


  ################## LINKS CREATION ##################
  links = NULL
  origins <- c()
  ends <- c()
  bandwidths <- c()
  bandwidthUnits <- c()
  distances <- c()
  distanceUnits <- c()

  # Link cells to M1 switches
  for (row in 1:nrow(m1Assoc)) {
    m1Coord <- m1Coords[which(m1Coords$group == m1Assoc[row,]$group),][1,]
    dis <- SDMTools::distance(lat1 = m1Assoc[row,]$lat,
                              lon1 = m1Assoc[row,]$lon,
                              lat2 = m1Coord$lat, lon2 = m1Coord$lon)$distance
    origins <- c(origins, paste("cell", row, sep = ""))
    ends <- c(ends, paste("m1_", m1Assoc[row,]$group, sep = ""))
    bandwidths <- c(bandwidths, 10)
    bandwidthUnits <- c(bandwidthUnits, "Gb/s")
    distances <- c(distances, dis)
    distanceUnits <- c(distanceUnits, "meters")
  }

  # Link the access ring elements
  accessRingLinks <- linkAccessRings(m1Assoc, m1Coords, m1AccAssocs,
                                     accCentCoords, m2Assocs)
  origins <- c(origins, accessRingLinks$from)
  ends <- c(ends, accessRingLinks$to)
  bandwidths <- c(bandwidths, accessRingLinks$bandwidth)
  bandwidthUnits <- c(bandwidthUnits, accessRingLinks$bandwidthUnits)
  distances <- c(distances, accessRingLinks$distance)
  distanceUnits <- c(distanceUnits, accessRingLinks$distanceUnits)

  # Link the aggregation ring elements
  aggRingLinks <- linkAggRings(m2Assocs, m2Switches, m2AggAssocs, aggCentCoords,
                               m3Assocs, m3Switches)
  origins <- c(origins, aggRingLinks$from)
  ends <- c(ends, aggRingLinks$to)
  bandwidths <- c(bandwidths, aggRingLinks$bandwidth)
  bandwidthUnits <- c(bandwidthUnits, aggRingLinks$bandwidthUnits)
  distances <- c(distances, aggRingLinks$distance)
  distanceUnits <- c(distanceUnits, aggRingLinks$distanceUnits)


  # Link all the M3 switches in a ring to ensure reachability
  if (nrow(m3Switches) > 1) {
    for (row in 2:nrow(m3Switches)) {
      dis <- SDMTools::distance(lat1 = m3Switches[row-1,]$lat,
                                lon1 = m3Switches[row-1,]$lon,
                                lat2 = m3Switches[row,]$lat,
                                lon2 = m3Switches[row,]$lon)$distance
      origins <- c(origins, paste("m3_", m3Switches[row-1,]$group, sep = ""))
      ends <- c(ends, paste("m3_rep_", m3Switches[row,]$group, sep = ""))
      bandwidths <- c(bandwidth, 36)
      bandwidthUnits <- c(bandwidthUnits, "Tb/s")
      distances <- c(distances, dis)
      distanceUnits <- c(distanceUnits, "meters")
    }

    # Connect last M3 to first replica M3
    lastM3 <- tail(m3Switches, n = 1)
    firstM3 <- head(m3Switches, n = 1)
    dis <- SDMTools::distance(lat1 = lastM3$lat, lon1 = lastM3$lon,
                              lat2 = firstM3$lat, lon2 = firstM3$lon)$distance
    origins <- c(origins, paste("m3_", lastM3$group, sep = ""))
    ends <- c(ends, paste("m3_rep_", firstM3$group, sep = ""))
    bandwidths <- c(bandwidths, 36)
    bandwidthUnits <- c(bandwidthUnits, "Tb/s")
    distances <- c(distances, dis)
    distanceUnits <- c(distanceUnits, "meters")
  }


  links <- data.frame(from = origins, to = ends, bandwidth = bandwidths,
                      bandwidthUnits = bandwidthUnits, distance = distances,
                      distanceUnits = distanceUnits)

  return(list(links = links, nodes = nodes))
  #### # Build and write the graph
  #### g = igraph::graph_from_data_frame(links, vertices = nodes, directed = FALSE)
  #### igraph::write_graph(graph = g, file = file, format = format)
}


#' @export
#' @name build5GScenario
#' @title build5GScenario
#' @description Builds the 5G scenario of Figure 1 in:
#' https://ieeexplore.ieee.org/mediastore_new/IEEE/content/media/8412472/8436592/8436847/paper_23-fig-1-source-large.gif
#' @param lats cell antenna latitudes
#' @param lons cell antenna longitudes
#' @return list of data.frames with associations
#' [[1]]  m1Assoc data.frame cell coordinates and the M1 ids
#' [[2]]  m1Coords data.frame M1 coordinates and their ids
#' [[3]]  m1AccAssocs data.frame M1 coordinates and the access ring ids
#' [[4]]  accCentCoords data.frame access ring center coordinates and their ids
#' [[5]]  m2Assocs data.frame access ring coords and associated M2 ids
#' [[6]]  m2Switches data.frame M2 switches coordinates and their ids
#' [[7]]  m2AggAssocs data.frame M2 coordinates and aggregation rings ids
#' [[8]]  aggCentCoords data.frame aggregation ring centers and their ids
#' [[9]]  m3Assocs data.frame aggregation ring centers and M3 ids
#' [[10]] m3Switches data.frame with M3 coordinates and their ids
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

  return(list(m1Assoc, m1Coords, m1AccAssocs, accCentCoords, m2Assocs,
              m2Switches, m2AggAssocs, aggCentCoords, m3Assocs, m3Switches))
}




######### TEST
# coboCells <- mecgen::cobo
# assocs <- build5GScenario(lats = coboCells$lat, lons = coboCells$lon)
#
# m1Assoc <- assocs[[1]]
# m1Coords <- assocs[[2]]
# m1AccAssocs <- assocs[[3]]
# accCentCoords <- assocs[[4]]
# m2Assocs <- assocs[[5]]
# m2Switches <- assocs[[6]]
# m2AggAssocs <- assocs[[7]]
# aggCentCoords <- assocs[[8]]
# m3Assocs <- assocs[[9]]
# m3Switches <- assocs[[10]]
#
#
# frames <- graphFrames(m1Assoc, m1Coords, m1AccAssocs, accCentCoords,
#                            m2Assocs, m2Switches, m2AggAssocs, aggCentCoords,
#                            m3Assocs, m3Switches)

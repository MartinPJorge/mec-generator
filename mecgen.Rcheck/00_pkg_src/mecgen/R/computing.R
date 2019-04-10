#' @name stackAndFill
#' @title stackAndFill
#' @description Stack two dataframes and fill with NULL values not present
#' @param dfA data.frame A
#' @param dfB data.frame B
#' @return stacked data.frame
stackAndFill <- function(dfA, dfB) {
  props <- list()

  # Stack common columns
  for (prop in intersect(names(dfA), names(dfB))) {
    props[[prop]] <- c(as.vector(dfA[[prop]]), as.vector(dfB[[prop]]))
  }

  # Fill in dfB those elements of dfA that it doesn't have
  for (prop in setdiff(names(dfA), names(dfB))) {
    props[[prop]] <- c(as.vector(dfA[[prop]]), rep(x = NULL, times = nrow(dfB)))
  }

  # Fill in dfA those elements of dfB that it doesn't have
  for (prop in setdiff(names(dfB), names(dfA))) {
    props[[prop]] <- c(rep(x = NULL, times = nrow(dfA)), as.vector(dfB[[prop]]))
  }

  stacked <- data.frame()
  for (prop in props) {
    stacked[[prop]] <- props[[prop]]
  }

  return(stacked)
}


#' @export
#' @name attachServers
#' @title attachServers
#' @description Attaches numServers to the given infrastructure switches
#' @param nodes data.frame(id, type, lon, lat)
#' @param link data.frame(from, to, bandwidth, bandwidthUnits, distance,
#' distanceUnits)
#' @param numServers number of servers to be attached
#' @param delay delay for link between each server and the attached switch
#' @param bandwidth bandwidth for link between each server and the switch
#' @param bandwidthUnits string for bandwidth units
#' @param distance distance beween servers and switches
#' @param distanceUnits string for distance units
#' @param switchType string to determine the type of switch to attach
#' @param properties list with each server properties, e.g., list(cpu=,disk=)
#' @param idPrefix string to be prefixed before the server id
#' @return list(nodes = data.frame, links = data.frame)
#' @note the function iterates through all the switches to assign the servers
attachServers <- function(nodes, links, numServers, delay, bandwidth,
                          bandwidthUnits, distance, distanceUnits, switchType,
                          properties, idPrefix) {
  switches <- nodes[which(nodes$type == switchType),]
  if (nrow(switches) == 0) {
    warning(paste("Trying to attach server to non existing switches of type ",
                  switchType, sep = ""))
    return(NULL)
  }

  if (nrow(subset(nodes, grepl(idPrefix, id))) > 0) {
    warning(paste("There are already node ids containing prefix ", idPrefix,
                  sep = ""))
    return(NULL)
  }

  # Assign servers to switches
  servIds <- c()
  switchIds <- c()
  servLons <- c()
  servLats <- c()
  for (servNum in 0:(numServers - 1)) {
    switchIdx <- servNum %% nrow(switches) + 1
    servIds <- c(servIds, paste(idPrefix, "_server_", servNum, sep = ""))
    switchIds <- c(switchIds, switches$id[switchIdx])
    servLons <- c(servLons, switches$lon[switchIdx])
    servLats <- c(servLats, switches$lat[switchIdx])
  }

  # Create the server nodes data.frame
  serverNodes <- data.frame(id = servIds, lon = servLons, lat = servLats,
                      type = rep(x = "server", times = numServers))
  for (att in attributes(properties)$names) {
    nodes[,att] <- properties[[att]]
  }

  # Create the server links data.frame
  serverLinks <- data.frame(from = servIds, to = switchIds,
                      bandwidth = rep(x = bandwidth, times = numServers),
                      bandwidthUnits = rep(x = bandwidthUnits,
                                           times = numServers),
                      distance = rep(x = distance, times = numServers),
                      distanceUnits = rep(x = distanceUnits,
                                          times = numServers))

  return(list(nodes = stackAndFill(nodes, serverNodes),
              links = stackAndFill(links, serverLinks)))
}


#' @export
#' @name genFogNodes
#' @title genFogNodes
#' @description Generates fog computing nodes within a certain region
#' @param latB region bottom latitude limit
#' @param latT region top latitude limit
#' @param lonL region left longitude limit
#' @param lonR region right longitude limit
#' @param cells data.frame(lon, lat) cell antennas' coordinates
#' @param numNodes number of fog compute nodes to generate
#' @param properties list with compute node resource properties
#' @return data.frame(lon, lat, cellLon, cellLat, cellPos, fogId,
#' prop1, ..., propN)
genFogNodes <- function(latB, latT, lonL, lonR, cells, numNodes, properties) {
  fogLats <- runif(n = numNodes, min = latB, max = latT)
  fogLons <- runif(n = numNodes, min = lonL, max = lonR)
  fogIds <- c()
  assignedLon <- c()
  assignedLat <- c()
  assignedPos <- c()

  for (fog in 1:numNodes) {
    fogLon <- fogLons[fog]
    fogLat <- fogLats[fog]
    fogIds <- c(fogIds, paste("fogNode_", fog, sep = ""))

    # Find nearest cell to be assigned
    assigned <- 1
    minDis <- Inf
    for (row in 1:nrow(cells)) {
      cell <- cells[row,]
      currDis <- SDMTools::distance(lat1 = fogLat, lon1 = fogLon,
                                    lat2 = cell$lat, lon2 = cell$lon)$distance
      if (currDis < minDis) {
        assigned <- row
        minDis <- currDis
      }
    }

    assignedPos <- c(assignedPos, assigned)
    assignedLon <- c(assignedLon, cells[assigned,]$lon)
    assignedLat <- c(assignedLat, cells[assigned,]$lat)
  }

  fogNodes <- data.frame(lon = fogLons, lat = fogLats, cellLon = assignedLon,
                    cellLat = assignedLat, cellPos = assignedPos,
                    fogId = fogIds)

  for (att in attributes(properties)$names) {
    fogNodes[,att] <- rep(x = properties[[att]], n = numNodes)
  }

  return(fogNodes)
}


#' @export
#' @name genEndpoints
#' @title genEndpoints
#' @description Generates endpoints within the selected region, asking for
#' certain bandwidth resources
#' @param latB region bottom latitude limit
#' @param latT region top latitude limit
#' @param lonL region left longitude limit
#' @param lonR region right longitude limit
#' @param cells data.frame(lon, lat) cell antennas' coordinates
#' @param numEndpoints number of endpoints to be generated
#' @param bandwidth bandwidth resources available for each user
#' @return data.frame(lon, lat, cellLon, cellLat, cellPos, enpointId, bandwidth)
genEndpoints <- function(latB, latT, lonL, lonR, cells, numEndpoints,
                         bandwidth) {
  endsLons <- runif(n = numEndpoints, min = lonL, max = lonR)
  endsLats <- runif(n = numEndpoints, min = latB, max = latT)
  endpointIds <- c()
  assignedLon <- c()
  assignedLat <- c()
  assignedPos <- c()

  for (endpoint in 1:numEndpoints) {
    endpointLon <- endsLons[endpoint]
    endpointLat <- endsLats[endpoint]
    endpointIds <- c(endpointIds, paste("endpoint_", endpoint, sep = ""))

    # Find nearest cell to be assigned
    assigned <- 1
    minDis <- Inf
    for (row in 1:nrow(cells)) {
      cell <- cells[row,]
      currDis <- SDMTools::distance(lat1 = endpointLat, lon1 = endpointLon,
                                    lat2 = cell$lat, lon2 = cell$lon)$distance
      if (currDis < minDis) {
        assigned <- row
        minDis <- currDis
      }
    }

    assignedPos <- c(assignedPos, assigned)
    assignedLon <- c(assignedLon, cells[assigned,]$lon)
    assignedLat <- c(assignedLat, cells[assigned,]$lat)
  }

  return(data.frame(lon = endsLons, lat = endsLats, cellLon = assignedLon,
                    cellLat = assignedLat, cellPos = assignedPos,
                    endpointId = endpointIds,
                    bandwidth = rep(x = bandwidth, times = numEndpoints)))
}



# nodes <- mecgen::coboNodes
# links <- mecgen::coboLinks
#
# nodes[["strangeVal"]] <- rep(x = 23, times = nrow(nodes))
# delay <- 20
# bandwidth <- 1
# bandwidthUnits <- "Gb/s"
# distance <- 2
# distanceUnits <- "meters"
# properties <- list(cpu = 2, disk = 100, mem = 16)
# numServers <- 5
#
# attachFrames <- attachServers(nodes = nodes, links = links,
#                               numServers = numServers, delay = delay,
#                               bandwidth = bandwidth,
#                               bandwidthUnits = bandwidthUnits,
#                               distance = distance,
#                               distanceUnits = distanceUnits, switchType = "m2",
#                               properties = properties, idPrefix = "dell_")
# newNodes <- attachFrames$nodes
# newLinks <- attachFrames$links
# servNodes <- tail(newNodes, n = numServers)
# servLinks <- tail(newNodes, n = numServers)
#

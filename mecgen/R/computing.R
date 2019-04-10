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
    if (is.numeric(dfA[[prop]][1])) {
      props[[prop]] <- c(as.vector(dfA[[prop]]), rep(x = 0, times = nrow(dfB)))
    } else {
      props[[prop]] <- c(as.vector(dfA[[prop]]), rep(x = "", times = nrow(dfB)))
    }
  }

  # Fill in dfA those elements of dfB that it doesn't have
  for (prop in setdiff(names(dfB), names(dfA))) {
    if (is.numeric(dfB[[prop]][1])) {
      props[[prop]] <- c(rep(x = 0, times = nrow(dfA)), as.vector(dfB[[prop]]))
    } else {
      props[[prop]] <- c(rep(x = "", times = nrow(dfA)), as.vector(dfB[[prop]]))
    }
  }

  return(as.data.frame(props))
}


#' @export
#' @name attachServers
#' @title attachServers
#' @description Attaches numServers to the given infrastructure switches
#' @param nodes data.frame(id, type, lon, lat)
#' @param link data.frame(from, to, bandwidth, bandwidthUnits, distance,
#' distanceUnits)
#' @param numServers number of servers to be attached
#' @param bandwidth bandwidth for link between each server and the switch
#' @param bandwidthUnits string for bandwidth units
#' @param distance distance beween servers and switches
#' @param distanceUnits string for distance units
#' @param switchType string to determine the type of switch to attach
#' @param properties list with each server properties, e.g., list(cpu=,disk=)
#' @param idPrefix string to be prefixed before the server id
#' @return list(nodes = data.frame, links = data.frame)
#' @note the function iterates through all the switches to assign the servers
attachServers <- function(nodes, links, numServers, bandwidth, bandwidthUnits,
                          distance, distanceUnits, switchType, properties,
                          idPrefix) {
  switches <- nodes[which(nodes$type == switchType),]
  if (nrow(switches) == 0) {
    warning(paste("Trying to attach server to non existing switches of type ",
                  switchType, sep = ""))
    return(NULL)
  }

  if (nrow(subset(switches, grepl(idPrefix, id))) > 0) {
    warning(paste("There are already server ids containing prefix ", idPrefix,
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
    switchIds <- c(switchIds, as.vector(switches$id)[switchIdx])
    servLons <- c(servLons, switches$lon[switchIdx])
    servLats <- c(servLats, switches$lat[switchIdx])
  }

  # Create the server nodes data.frame
  serverNodes <- data.frame(id = servIds, lon = servLons, lat = servLats,
                      type = rep(x = "server", times = numServers))
  for (att in attributes(properties)$names) {
    serverNodes[,att] <- properties[[att]]
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
#' @name attachFogNodes
#' @title attachFogNodes
#' @description Generates fog computing nodes within a certain region, and
#' attaches them to the existing infrastructure
#' @param nodes data.frame(id, type, lon, lat)
#' @param link data.frame(from, to, bandwidth, bandwidthUnits, distance,
#' distanceUnits)
#' @param latB region bottom latitude limit
#' @param latT region top latitude limit
#' @param lonL region left longitude limit
#' @param lonR region right longitude limit
#' @param numNodes number of fog compute nodes to generate
#' @param properties list with compute node resource properties
#' @param bandwidth bandwidth for link between each server and the switch
#' @param bandwidthUnits string for bandwidth units
#' @param distanceUnits string for distance units
#' @param idPrefix string to be prefixed before the fog node id
#' @return list(nodes = data.frame, links = data.frame)
attachFogNodes <- function(nodes, links, latB, latT, lonL, lonR, numNodes,
                        properties, bandwidth, bandwidthUnits, idPrefix) {
  cells <- nodes[which(nodes$type == "cell"),]
  if (nrow(cells) == 0) {
    warning("No cells within the nodes data.frame")
    return(NULL)
  }
  if (nrow(subset(cells, grepl(idPrefix, id))) > 0) {
    warning(paste("There are already fog node ids containing prefix ", idPrefix,
                  sep = ""))
    return(NULL)
  }

  fogLats <- runif(n = numNodes, min = latB, max = latT)
  fogLons <- runif(n = numNodes, min = lonL, max = lonR)
  fogIds <- c()
  assignedCell <- c()
  distances <- c()

  for (fog in 1:numNodes) {
    fogLon <- fogLons[fog]
    fogLat <- fogLats[fog]
    fogIds <- c(fogIds, paste(idPrefix, "_fogNode_", fog, sep = ""))

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

    assignedCell <- c(assignedCell, cells[assigned,]$id)
    distances <- c(distances, minDis)
  }

  # Create the nodes data.frame
  fogNodes <- data.frame(id = fogIds, type = "fogNode", lon = fogLons,
                      lat = fogLats)
  for (att in attributes(properties)$names) {
    fogNodes[,att] <- rep(x = properties[[att]], n = numNodes)
  }

  # Create the links data.frame
  fogLinks <- data.frame(from = fogIds, to = assignedCell,
                      bandwidth = rep(x = bandwidth, times = numNodes),
                      bandwidthUnits = rep(x = bandwidthUnits,
                                           times = numNodes),
                      distance = distances,
                      distanceUnits = rep(x = "meters", times = numNodes))

  return(list(nodes = stackAndFill(nodes, fogNodes),
              links = stackAndFill(links, fogLinks)))
}


#' @export
#' @name attachEndpoints
#' @title attachEndpoints
#' @description Generates endpoints within the selected region, asking for
#' certain bandwidth resources. And attaches them to the nearest cells.
#' @param nodes data.frame(id, type, lon, lat)
#' @param link data.frame(from, to, bandwidth, bandwidthUnits, distance,
#' distanceUnits)
#' @param latB region bottom latitude limit
#' @param latT region top latitude limit
#' @param lonL region left longitude limit
#' @param lonR region right longitude limit
#' @param numEndpoints number of endpoints to be generated
#' @param bandwidth bandwidth resources available for each user
#' @param bandwidthUnits string for bandwidth units
#' @param idPrefix string to be prefixed before the fog node id
#' @return list(nodes = data.frame, links = data.frame)
attachEndpoints <- function(nodes, links, latB, latT, lonL, lonR, numEndpoints,
                         bandwidth, bandwidthUnits, idPrefix) {
  cells <- nodes[which(nodes$type == "cell"),]
  if (nrow(cells) == 0) {
    warning("No cells within the nodes data.frame")
    return(NULL)
  }
  if (nrow(subset(cells, grepl(idPrefix, id))) > 0) {
    warning(paste("There are already endpoint ids containing prefix ", idPrefix,
                  sep = ""))
    return(NULL)
  }

  endsLons <- runif(n = numEndpoints, min = lonL, max = lonR)
  endsLats <- runif(n = numEndpoints, min = latB, max = latT)
  endpointIds <- c()
  assignedCells <- c()
  distances <- c()

  for (endpoint in 1:numEndpoints) {
    endpointLon <- endsLons[endpoint]
    endpointLat <- endsLats[endpoint]
    endpointIds <- c(endpointIds, paste(idPrefix, "_endpoint_", endpoint,
                                        sep = ""))

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

    distances <- c(distances, minDis)
    assignedCells <- c(assignedCells, cells[assigned,]$id)
  }

  # Create the nodes data.frame
  endpointNodes <- data.frame(id = endpointIds,
                          type = rep(x = "endpoint", times = numEndpoints),
                          lon = endsLons, lat = endsLats)

  # Create the links data.frame
  endpointLinks <- data.frame(from = endpointIds, to = assignedCells,
                      bandwidth = rep(x = bandwidth, times = numEndpoints),
                      bandwidthUnits = rep(x = bandwidthUnits,
                                           times = numEndpoints),
                      distance = distances,
                      distanceUnits = rep(x = "meters", times = numEndpoints))

  return(list(nodes = stackAndFill(nodes, endpointNodes),
              links = stackAndFill(links, endpointLinks)))
}



#  nodes <- mecgen::coboNodes
#  links <- mecgen::coboLinks
#
#  nodes[["strangeVal"]] <- rep(x = 23, times = nrow(nodes))
#  bandwidth <- 1
#  bandwidthUnits <- "Gb/s"
#  distance <- 2
#  distanceUnits <- "meters"
#  properties <- list(cpu = 2, disk = 100, mem = 16, rare = 20)
#  numServers <- 5
#
#  attachFrames <- attachServers(nodes = nodes, links = links,
#                                numServers = numServers,
#                                bandwidth = bandwidth,
#                                bandwidthUnits = bandwidthUnits,
#                                distance = distance,
#                                distanceUnits = distanceUnits,
#                                switchType = "m2",
#                                properties = properties, idPrefix = "dell_")
#  newNodes <- attachFrames$nodes
#  newLinks <- attachFrames$links
#  servNodes <- tail(newNodes, n = numServers)
#  servLinks <- tail(newLinks, n = numServers)

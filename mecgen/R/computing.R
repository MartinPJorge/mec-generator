#' @export
#' @name attachServers
#' @title attachServers
#' @description Attaches numServers to the given infrastructure switches
#' @param switches data.frame(lon, lat, group) with switches coordinates and ids
#' @param numServers number of servers to be attached
#' @param delay delay for link between each server and the attached switch
#' @param bandwidth bandwidth for link between each server and the switch
#' @param properties list with each server properties, e.g., list(cpu=,disk=)
#' @return data.frame(properties.a, ..., properties.b, delay, bandwidth,
#' serverId, switchGroup)
#' @note the function iterates through all the switches to assign the servers
attachServers <- function(switches, numServers, delay, bandwidth, properties) {
  # Assign servers to switches
  servIds <- c()
  switchIds <- c()
  for (servNum in 0:(numServers - 1)) {
    servIds <- c(servIds, paste("server_", servNum, sep = ""))
    switchIds <- c(switchIds,
                   switches$group[servNum %% nrow(switches) + 1])
  }

  attachments <- data.frame(serverId = servIds, switchGroup = switchIds,
                            delay = rep(x = delay, times = numServers),
                            bandwidth = rep(x = bandwidth, times = numServers))
  for (att in attributes(properties)$names) {
    attachments[,att] <- properties[[att]]
  }

  return(attachments)
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


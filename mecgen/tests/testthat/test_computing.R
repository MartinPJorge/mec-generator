context("Computing assignments")
library(mecgen)

test_that("Server assignments work", {
  nodes <- data(coboNodes)
  links <- data(coboLinks)

  nodes[["strangeVal"]] <- rep(x = 23, times = nrow(nodes))
  bandwidth <- 1
  bandwidthUnits <- "Gb/s"
  distance <- 2
  distanceUnits <- "meters"
  properties <- list(cpu = 2, disk = 100, mem = 16, rare = 20)
  numServers <- 5
  switchType <- "m2"

  attachFrames <- attachServers(nodes = nodes, links = links,
                                numServers = numServers,
                                bandwidth = bandwidth,
                                bandwidthUnits = bandwidthUnits,
                                distance = distance,
                                distanceUnits = distanceUnits,
                                switchType = switchType,
                                properties = properties, idPrefix = "dell")
  newNodes <- attachFrames$nodes
  newLinks <- attachFrames$links
  servNodes <- tail(newNodes, n = numServers)
  servLinks <- tail(newLinks, n = numServers)
  serverIds <- c()
  for (id in 0:(numServers - 1)) {
    serverIds <- c(serverIds, paste("dell_server_", id, sep = ""))
  }


  expect_equal(as.vector(servLinks$bandwidth),
               rep(x = bandwidth, times = numServers))
  expect_equal(as.vector(servNodes$strangeVal),
               rep(x = 0, times = numServers))
  expect_equal(as.vector(servNodes$cpu),
               rep(x = properties$cpu, times = numServers))
  expect_equal(as.vector(servNodes$disk),
               rep(x = properties$disk, times = numServers))
  expect_equal(as.vector(servNodes$mem),
               rep(x = properties$mem, times = numServers))
  expect_equal(as.vector(servNodes$rare),
               rep(x = properties$rare, times = numServers))
  expect_equal(head(newNodes, n = nrow(newNodes) - numServers)$rare,
               rep(x = 0, times = nrow(newNodes) - numServers))
  expect_equal(as.vector(servNodes$type),
               rep(x = "server", times = numServers))
  expect_equal(as.vector(servNodes$id), serverIds)
})


test_that("Fog nodes generation works", {
  nodes <- data(coboNodes)
  links <- data(coboLinks)

  nodes[["strangeVal"]] <- rep(x = 23, times = nrow(nodes))
  coboLonL <- -3.775409
  coboLonR <- -3.737324
  coboLatB <- 40.253541
  coboLatT <- 40.276686
  bandwidth <- 1
  bandwidthUnits <- "Gb/s"
  properties <- list(cpu = 2, disk = 100, mem = 16, rare = 20)
  numFogNodes <- 3

  attachFrames <- attachFogNodes(nodes = nodes, links = links,
                                 latB = coboLatB, latT = coboLatT,
                                 lonL = coboLonL, lonR = coboLonR,
                                 numNodes = numFogNodes,
                                 properties = properties, bandwidth = bandwidth,
                                 bandwidthUnits = bandwidthUnits,
                                 idPrefix = "test")

  newNodes <- attachFrames$nodes
  newLinks <- attachFrames$links
  fogNodes <- tail(newNodes, n = numFogNodes)
  fogLinks <- tail(newLinks, n = numFogNodes)
  fogIds <- c()
  for (id in 1:numFogNodes) {
    fogIds <- c(fogIds, paste("test_fogNode_", id, sep = ""))
  }


  expect_equal(as.vector(fogLinks$bandwidth),
               rep(x = bandwidth, times = numFogNodes))
  expect_equal(as.vector(fogNodes$strangeVal),
               rep(x = 0, times = numFogNodes))
  expect_equal(as.vector(fogNodes$cpu),
               rep(x = properties$cpu, times = numFogNodes))
  expect_equal(as.vector(fogNodes$disk),
               rep(x = properties$disk, times = numFogNodes))
  expect_equal(as.vector(fogNodes$mem),
               rep(x = properties$mem, times = numFogNodes))
  expect_equal(as.vector(fogNodes$rare),
               rep(x = properties$rare, times = numFogNodes))
  expect_equal(head(newNodes, n = nrow(newNodes) - numFogNodes)$rare,
               rep(x = 0, times = nrow(newNodes) - numFogNodes))
  expect_equal(as.vector(fogNodes$type),
               rep(x = "fogNode", times = numFogNodes))
  expect_equal(as.vector(fogNodes$id), fogIds)
})


test_that("Fog nodes generation works with multiple cell assignments", {

  # Build up the nodes and links for just 6 cells
  cells <- 6
  coboCells <- data(cobo)
  regions <-   data(regions)
  coboCells <- head(coboCells, n = cells) # just 6 AAUs
  assocs <- build5GScenario(lats = coboCells$lat, lons = coboCells$lon)
  m1Assoc <- assocs[[1]]
  m1Coords <- assocs[[2]]
  m1AccAssocs <- assocs[[3]]
  accCentCoords <- assocs[[4]]
  m2Assocs <- assocs[[5]]
  m2Switches <- assocs[[6]]
  m2AggAssocs <- assocs[[7]]
  aggCentCoords <- assocs[[8]]
  m3Assocs <- assocs[[9]]
  m3Switches <- assocs[[10]]
  frames <- graphFrames(m1Assoc, m1Coords, m1AccAssocs, accCentCoords,
                        m2Assocs, m2Switches, m2AggAssocs, aggCentCoords,
                        m3Assocs, m3Switches)

  nodes <- frames$nodes
  links <- frames$links

  nodes[["strangeVal"]] <- rep(x = 23, times = nrow(nodes))
  coboLonL <- -3.775409
  coboLonR <- -3.737324
  coboLatB <- 40.253541
  coboLatT <- 40.276686
  bandwidth <- 1
  bandwidthUnits <- "Gb/s"
  properties <- list(cpu = 2, disk = 100, mem = 16, rare = 20)
  numFogNodes <- 3
  dis <- 100000 # 100 km of distance expressed in meters -> all cells assigned

  attachFrames <- attachFogNodes(nodes = nodes, links = links,
                                 latB = coboLatB, latT = coboLatT,
                                 lonL = coboLonL, lonR = coboLonR,
                                 numNodes = numFogNodes,
                                 properties = properties, bandwidth = bandwidth,
                                 bandwidthUnits = bandwidthUnits,
                                 idPrefix = "test", dis = dis)

  newNodes <- attachFrames$nodes
  newLinks <- attachFrames$links
  fogNodes <- tail(newNodes, n = numFogNodes)
  fogLinks <- tail(newLinks, n = numFogNodes * cells)
  fogIds <- c()
  for (id in 1:numFogNodes) {
    fogIds <- c(fogIds, paste("test_fogNode_", id, sep = ""))
  }


  expect_equal(as.vector(fogLinks$bandwidth),
               rep(x = bandwidth, times = numFogNodes * cells))
  expect_equal(as.vector(fogNodes$strangeVal),
               rep(x = 0, times = numFogNodes))
  expect_equal(as.vector(fogNodes$cpu),
               rep(x = properties$cpu, times = numFogNodes))
  expect_equal(as.vector(fogNodes$disk),
               rep(x = properties$disk, times = numFogNodes))
  expect_equal(as.vector(fogNodes$mem),
               rep(x = properties$mem, times = numFogNodes))
  expect_equal(as.vector(fogNodes$rare),
               rep(x = properties$rare, times = numFogNodes))
  expect_equal(head(newNodes, n = nrow(newNodes) - numFogNodes)$rare,
               rep(x = 0, times = nrow(newNodes) - numFogNodes))
  expect_equal(as.vector(fogNodes$type),
               rep(x = "fogNode", times = numFogNodes))
  expect_equal(as.vector(fogNodes$id), fogIds)
})


test_that("Fog node endpoints are generated correctly", {

  # Build up the nodes and links for just 6 cells
  cells <- 6
  coboCells <- data(cobo)
  regions <-   data(regions)
  coboCells <- head(coboCells, n = cells) # just 6 AAUs
  assocs <- build5GScenario(lats = coboCells$lat, lons = coboCells$lon)
  m1Assoc <- assocs[[1]]
  m1Coords <- assocs[[2]]
  m1AccAssocs <- assocs[[3]]
  accCentCoords <- assocs[[4]]
  m2Assocs <- assocs[[5]]
  m2Switches <- assocs[[6]]
  m2AggAssocs <- assocs[[7]]
  aggCentCoords <- assocs[[8]]
  m3Assocs <- assocs[[9]]
  m3Switches <- assocs[[10]]
  frames <- graphFrames(m1Assoc, m1Coords, m1AccAssocs, accCentCoords,
                        m2Assocs, m2Switches, m2AggAssocs, aggCentCoords,
                        m3Assocs, m3Switches)

  nodes <- frames$nodes
  links <- frames$links

  nodes[["strangeVal"]] <- rep(x = 23, times = nrow(nodes))
  coboLonL <- -3.775409
  coboLonR <- -3.737324
  coboLatB <- 40.253541
  coboLatT <- 40.276686
  bandwidth <- 1
  bandwidthUnits <- "Gb/s"
  properties <- list(cpu = 2, disk = 100, mem = 16, rare = 20)
  numFogNodes <- 3
  dis <- 100000 # 100 km of distance expressed in meters -> all cells assigned

  attachFrames <- attachFogEndpoints(nodes = nodes, links = links,
                                     latB = coboLatB, latT = coboLatT,
                                     lonL = coboLonL, lonR = coboLonR,
                                     numNodes = numFogNodes,
                                     properties = properties,
                                     bandwidth = bandwidth,
                                     bandwidthUnits = bandwidthUnits,
                                     idPrefix = "test", dis = dis)

  newNodes <- attachFrames$nodes
  newLinks <- attachFrames$links
  fogNodes <- tail(newNodes, n = 2*numFogNodes)
  fogLinks <- tail(newLinks, n = numFogNodes * (1 + cells))
  fogIds <- c()
  for (id in 1:numFogNodes) {
    fogIds <- c(fogIds, paste("fogEndpoint_test_fogNode_", id, sep = ""))
  }
  for (id in 1:numFogNodes) {
    fogIds <- c(fogIds, paste("fogEndpoint_test_endpoint_", id, sep = ""))
  }

  expect_equal(as.vector(fogNodes$id), fogIds)
  expect_equal(as.vector(tail(fogLinks, n = numFogNodes)$distance),
               rep(x = 0, times = numFogNodes))
  expect_equal(as.vector(head(fogNodes, n = numFogNodes)$type),
               rep(x = "fogNode", times = numFogNodes))
  expect_equal(as.vector(tail(fogNodes, n = numFogNodes)$type),
                rep(x = "endpoint", times = numFogNodes))
  expect_equal(as.vector(tail(fogNodes, n = numFogNodes)$lon),
               as.vector(head(fogNodes, n = numFogNodes)$lon))
  expect_equal(as.vector(tail(fogNodes, n = numFogNodes)$lat),
               as.vector(head(fogNodes, n = numFogNodes)$lat))
})


test_that("Endpoints generation works", {
  nodes <- data(coboNodes)
  links <- data(coboLinks)

  nodes[["strangeVal"]] <- rep(x = 23, times = nrow(nodes))
  coboLonL <- -3.775409
  coboLonR <- -3.737324
  coboLatB <- 40.253541
  coboLatT <- 40.276686
  bandwidth <- 1
  bandwidthUnits <- "Gb/s"
  numEndpoints <- 3
  endpointIds <- c()
  for (endpoint in 1:numEndpoints) {
    endpointIds <- c(endpointIds, paste("person_endpoint_", endpoint, sep = ""))
  }

  attachFrames <- attachEndpoints(nodes = nodes, links = links, latB = coboLatB,
                                  latT = coboLatT, lonL = coboLonL,
                                  lonR = coboLonR, numEndpoints = numEndpoints,
                                  bandwidth = bandwidth,
                                  bandwidthUnits = bandwidthUnits,
                                  idPrefix = "person")
  newNodes <- attachFrames$nodes
  newLinks <- attachFrames$links
  endpointNodes <- tail(newNodes, n = numEndpoints)
  endpointLinks <- tail(newLinks, n = numEndpoints)


  expect_equal(as.vector(endpointLinks$bandwidth),
               rep(x = bandwidth, times = numEndpoints))
  expect_equal(as.vector(endpointNodes$strangeVal),
               rep(x = 0, times = numEndpoints))
  expect_equal(as.vector(endpointNodes$type),
               rep(x = "endpoint", times = numEndpoints))
  expect_equal(as.vector(endpointNodes$id), endpointIds)
})

context("Computing assignments")
library(mecgen)

test_that("Server assignments work", {
  nodes <- mecgen::coboNodes
  links <- mecgen::coboLinks

  nodes["strangeVal",] <- rep(x = 23, times = nrow(nodes))
  delay <- 20
  bandwidth <- 1
  bandwidthUnits <- "Gb/s"
  distance <- 2
  distanceUnits <- "meters"
  properties <- list(cpu = 2, disk = 100, mem = 16)
  numServers <- 5

  attachFrames <- attachServers(nodes = nodes, links = links,
                          numServers = rep(x = numServers, times = numServers),
                          delay = rep(x = delay, times = numServers),
                          bandwidth = rep(x = bandwidth, times = numServers),
                          bandwidthUnits = rep(x = bandwidthUnits,
                                               times = numServers),
                          distance = rep(x = distance, times = numServers),
                          distanceUnits = rep(x = distanceUnits,
                                              times = numServers),
                          switchType = "m2", properties = properties,
                          idPrefix = "dell_")
  newNodes <- attachFrames$nodes
  newLinks <- attachFrames$links
  servNodes <- tail(newNodes, n = numServers)
  servLinks <- tail(newNodes, n = numServers)

  expect_equal(newNodes$strangeVal, rep(x = NULL, times = numServers))
  expect_equal(newNodes$delay, rep(x = delay, times = numServers))
  expect_equal(newNodes$bandwidth, rep(x = bandwidth, times = numServers))
  expect_equal(newNodes$cpu, rep(x = properties$cpu, times = numServers))
  expect_equal(newNodes$disk, rep(x = properties$disk, times = numServers))
  expect_equal(newNodes$mem, rep(x = properties$mem, times = numServers))
  expect_equal(as.vector(attachFrames$switchType),
               rep(x = switchType, times = numServers))
  expect_equal(as.vector(attachFrames$serverId), serverIds)
  expect_equal(attachFrames$switchGroup, c(1,2,3,1,2))
})


test_that("Fog nodes generation works", {
  data(cobo)
  tenCells <- head(cobo, n = 10)
  coboLonL <- -3.775409
  coboLonR <- -3.737324
  coboLatB <- 40.253541
  coboLatT <- 40.276686
  numNodes <- 20
  properties <- list(cpu = 2, mem = 16, disk = 100)
  fogIds <- c()
  for (fog in 1:numNodes) {
    fogIds <- c(fogIds, paste("fogNode_", fog, sep = ""))
  }

  fogNodes <- genFogNodes(latB = coboLatB, latT = coboLatT, lonL = coboLonL,
                          lonR = coboLonR, cells = tenCells,
                          numNodes = numNodes, properties = properties)

  expect_equal(as.vector(fogNodes$fogId), fogIds)
  expect_equal(as.vector(fogNodes$cpu), rep(x = properties$cpu,
                                            times = numNodes))
  expect_equal(as.vector(fogNodes$mem), rep(x = properties$mem,
                                            times = numNodes))
  expect_equal(as.vector(fogNodes$disk), rep(x = properties$disk,
                                             times = numNodes))
})



test_that("Endpoints generation works", {
  data(cobo)
  tenCells <- head(cobo, n = 10)
  coboLonL <- -3.775409
  coboLonR <- -3.737324
  coboLatB <- 40.253541
  coboLatT <- 40.276686
  numEndpoints <- 20
  bandwidth <- 100
  endpointIds <- c()
  for (endpoint in 1:numEndpoints) {
    endpointIds <- c(endpointIds, paste("endpoint_", endpoint, sep = ""))
  }

  endpointNodes <- genEndpoints(latB = coboLatB, latT = coboLatT,
                                    lonL = coboLonL, lonR = coboLonR,
                                    cells = tenCells,
                                    numEndpoints = numEndpoints,
                                    bandwidth = bandwidth)

  expect_equal(as.vector(endpointNodes$endpointId), endpointIds)
  expect_equal(as.vector(endpointNodes$bandwidth), rep(x = bandwidth,
                                                       times = numEndpoints))
})

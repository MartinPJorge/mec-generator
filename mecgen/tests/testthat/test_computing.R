context("Computing assignments")
library(mecgen)

test_that("Server assignments work", {
  switches <- data.frame(group = c(1,2,3))
  delay <- 20
  bandwidth <- 1
  properties <- list(cpu = 2, disk = 100, mem = 16)
  numServers <- 5
  attachedServs <- attachServers(switches = switches, numServers = numServers,
                              delay = delay, bandwidth = bandwidth,
                              properties = properties)
  serverIds <- c()
  for (servNum in 0:(numServers - 1)) {
    serverIds <- c(serverIds, paste("server_", servNum, sep = ""))
  }

  expect_equal(attachedServs$delay, rep(x = delay, times = numServers))
  expect_equal(attachedServs$bandwidth, rep(x = bandwidth, times = numServers))
  expect_equal(attachedServs$cpu, rep(x = properties$cpu, times = numServers))
  expect_equal(attachedServs$disk, rep(x = properties$disk, times = numServers))
  expect_equal(attachedServs$mem, rep(x = properties$mem, times = numServers))
  expect_equal(as.vector(attachedServs$serverId), serverIds)
  expect_equal(attachedServs$switchGroup, c(1,2,3,1,2))
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

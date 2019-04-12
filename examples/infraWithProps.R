library(mecgen)

coboCells <- mecgen::cobo
regions <- mecgen::regions
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

attachFrames <- attachServers(nodes = frames$nodes, links = frames$links,
                              numServers = 3,
                              bandwidth = 12,
                              bandwidthUnits = "Mbps",
                              distance = 0,
                              distanceUnits = "meter",
                              switchType = "m2",
                              properties = list(cpu=2, mem=20, disk=100), idPrefix = "dell")

# Retrieve the Cobo Calleja region limiting coordinates
coboRegion <- regions$regions[[2]]
coboLonL <- as.numeric(strsplit(coboRegion$bl, split = ",")[[1]][2])
coboLonR <- as.numeric(strsplit(coboRegion$tr, split = ",")[[1]][2])
coboLatB <- as.numeric(strsplit(coboRegion$bl, split = ",")[[1]][1])
coboLatT <- as.numeric(strsplit(coboRegion$tr, split = ",")[[1]][1])

attachFrames <- attachFogNodes(nodes = attachFrames$nodes,
                               links = attachFrames$links,
                               latB = coboLatB, latT = coboLatT,
                               lonL = coboLonL, lonR = coboLonR,
                               numNodes = 10,
                               properties = list(cpu = 2, mem = 20), bandwidth = 20,
                               bandwidthUnits = "Mpbs",
                               idPrefix = "test")


attachFrames <- attachEndpoints(nodes = attachFrames$nodes,
                                links = attachFrames$links, latB = coboLatB,
                                latT = coboLatT, lonL = coboLonL,
                                lonR = coboLonR, numEndpoints = 10,
                                bandwidth = 10,
                                bandwidthUnits = "Mbps",
                                idPrefix = "person")

# Attach links properties
froms = as.vector(attachFrames$links$from)
tos = as.vector(attachFrames$links$to)
properties = list(reliability = rep(x = 1, times = length(froms)),
                  trafficCost = rep(x = 10, times = length(froms)))
newLinks <- addLinkProps(links = as.data.frame(attachFrames$links),
                         from_ = froms, to_ = tos, properties = properties)


# Attach nodes properties
reliab <- rep(x = 1, times = nrow(attachFrames$nodes))
newNodes <- addNodeProps(nodes = as.data.frame(attachFrames$nodes),
                         id_ = as.vector(attachFrames$nodes$id),
                         properties = list(reliability = reliab))

cellIdxs <- which(attachFrames$nodes$type == "cell")
cellIdxs <- as.vector(attachFrames$nodes$id)[cellIdxs]
radioTechs <- rep(x = "LTE", times = length(cellIdxs))
radioCost <- rep(x = 20, times = length(radioTechs))

serverIdxs <- which(attachFrames$nodes$type == "server")
serverIdxs <- as.vector(attachFrames$nodes$id)[serverIdxs]
serversCost <- rep(x = 100, times = length(serverIdxs))

newNodes <- addNodeProps(nodes = newNodes,
                         id_ = cellIdxs, properties = list(radio = radioTechs,
                                           radioCost = radioCost))
newNodes <- addNodeProps(node = as.data.frame(newNodes), id_ = serverIdxs,
                         properties = list(resCost = serversCost))

# Write the graphcellsIdx

g = igraph::graph_from_data_frame(newLinks, vertices = newNodes,
                                  directed = FALSE)
igraph::write_graph(graph = g, file = "/tmp/5g-mec.gml", format = "gml")

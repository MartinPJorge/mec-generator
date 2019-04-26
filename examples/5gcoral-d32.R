library(mecgen)

#' @description generates the 5g coral d3.2 test scenario for the placement
#'              algorithm
#' @param cells number of antennas of mecgen::cobo (Cobo Calleja) dataset
#' @param m1_servers number of servers attached to m1 nodes
#' @param m2_servers number of servers attached to m2 nodes
#' @param fog_nodes number of fog nodes
#' @param m2_cpu ammount of CPU for the m2 servers
#' @param m2_mem ammount of memory for the m2 servers
#' @param m2_disk ammount of disk for the m2 servers
#' @param m1_cpu ammount of CPU for the m1 servers
#' @param m1_mem ammount of memory for the m1 servers
#' @param m1_disk ammount of disk for the m1 servers
#' @param fog_cpu ammount of CPU for the fog servers
#' @param fog_mem ammount of memory for the fog servers
#' @param fog_disk ammount of disk for the fog servers
#' @param m1_vol_mean central volatility of m1 servers
#' @param m1_vol_dif percentual difference from the m1_vol_mean
#' @param m2_vol_mean central volatility of m2 servers
#' @param m2_vol_dif percentual difference from the m2_vol_mean
#' @param out path of the GML file where the scenario is generated
#' @references 
#'  [RPi3B+] https://www.raspberrypi.org/products/raspberry-pi-3-model-b-plus/
#'  [edgeDBox] https://azure.microsoft.com/mediahandler/files/resourcefiles/azure-data-box-edge-datasheet/Azure%20Data%20Box%20Edge%20Datasheet.pdf
#'  [truncExp] https://math.stackexchange.com/questions/28004/random-exponential-like-distribution
coralD32 <- function(cells=36, m1_servers=6, m2_servers=1, fog_nodes=128,
                     m2_cpu=200, m2_mem=1280, m2_disk=120000,
                     m1_cpu=20, m1_mem=128, m1_disk=12000,
                     fog_cpu=1, fog_mem=1, fog_disk=2,
                     m1_vol_mean=0.01, m1_vol_dif=0.1,
                     fog_vol_mean=0.1, fog_vol_dif=0.1,
                     out="/tmp/5gcorald32.gml") {
    
  coboCells <- mecgen::cobo
  regions <- mecgen::regions
  coboCells <- head(coboCells, n = cells)
  
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
  
  # Attach the M1 and M2 server
  edgeIdPrefix <- "azure_data_box_edge_m1"
  attachFrames <- attachServers(nodes = frames$nodes, links = frames$links,
                                numServers = m1_servers,
                                bandwidth = 1024,
                                bandwidthUnits = "Mbps",
                                distance = 0,
                                distanceUnits = "meter",
                                switchType = "m1",
                                properties = list(cpu=m1_cpu, mem=m1_mem,
                                                  disk=m1_disk),
                                idPrefix = edgeIdPrefix)
  
  # x10 times an azure data box capabilities
  cloudIdPrefix <- "dell_m2"
  attachFrames <- attachServers(nodes = attachFrames$nodes, links = attachFrames$links,
                                numServers = m2_servers,
                                bandwidth = 1024,
                                bandwidthUnits = "Mbps",
                                distance = 0,
                                distanceUnits = "meter",
                                switchType = "m2",
                                properties = list(cpu=m2_cpu, mem=m2_mem,
                                                  disk=m2_disk),
                                idPrefix = cloudIdPrefix)
  cloud.servers <- subset(attachFrames$nodes, grepl(cloudIdPrefix, id))
  cloud.servers.idx <- which(attachFrames$nodes$id %in% cloud.servers$id)
  cloud.servers.ids <- as.vector(attachFrames$nodes[cloud.servers.idx,]$id)
  attachFrames$nodes <- addNodeProps(nodes=attachFrames$nodes,
                                     id_=cloud.servers.ids,
                                     properties=list(
                                       radio=rep(x="cloud", times=m2_servers)))
  
  # Retrieve the Cobo Calleja region limiting coordinates
  coboRegion <- regions$regions[[2]]
  coboLonL <- as.numeric(strsplit(coboRegion$bl, split = ",")[[1]][2])
  coboLonR <- as.numeric(strsplit(coboRegion$tr, split = ",")[[1]][2])
  coboLatB <- as.numeric(strsplit(coboRegion$bl, split = ",")[[1]][1])
  coboLatT <- as.numeric(strsplit(coboRegion$tr, split = ",")[[1]][1])
  
  # Attach 2 fog endpoints, and 1 fog compute node
  fogIdPrefix <- "rpi"
  attachFrames <- attachFogNodes(nodes = attachFrames$nodes,
                                 links = attachFrames$links,
                                 latB = coboLatB, latT = coboLatT,
                                 lonL = coboLonL, lonR = coboLonR,
                                 numNodes = fog_nodes,
                                 properties = list(cpu = 1, mem = 1, disk=2),
                                 bandwidth = 11, bandwidthUnits = "Mpbs",
                                 idPrefix = fogIdPrefix)
  
  
  
  #' @description generates exponentialy distributed values inside the interval
  #'              (center - center*diff, center + center*diff)
  #' @param center center of the interval
  #' @param diff percentual difference from the center to the limits of interval
  #' @param n number of values to generate
  #' @return vector with n numbers exponentially distributed in the specified
  #'         interval
  #' @references https://math.stackexchange.com/questions/28004/random-exponential-like-distribution
  intervalExp <- function(center, diff, n) {
    left <- center * (1 - diff)
    right <- center * (1 + diff)
    values <- rexp(n = n, rate = 1 / center)
    
    
    for (i in 1:length(values)) {
      if (values[i] > 0) {
        values[i] <- values[i] - floor(values[i])
      }
      values[i] <- values[i] * (right - left) + left
    }
    
    return(values)
  }
  
  
  ######## Fog nodes volatility ##########
  fogVolatility.values <- intervalExp(center = fog_vol_mean,
                                      diff = fog_vol_dif, n = fog_nodes)
  
  fog.nodes <- subset(attachFrames$nodes, grepl(fogIdPrefix, id))
  fog.nodes.idx <- which(attachFrames$nodes$id %in% fog.nodes$id)
  fog.nodes.ids <- as.vector(attachFrames$nodes[fog.nodes.idx,]$id)
  attachFrames$nodes <- addNodeProps(nodes=attachFrames$nodes, id_=fog.nodes.ids,
                                     properties=list(
                                       volatility=fogVolatility.values,
                                       radio=rep(x="fog",
                                                 times=length(fog.nodes.ids))))
  
  ######## Edge nodes volatility ########
  edgeVolatility.values <- intervalExp(center = m1_vol_mean,
                                       diff = m1_vol_dif,
                                       n = m1_servers)
  edge.servers <- subset(attachFrames$nodes, grepl(edgeIdPrefix, id))
  edge.servers.idx <- which(attachFrames$nodes$id %in% edge.servers$id)
  edge.servers.ids <- as.vector(attachFrames$nodes[edge.servers.idx,]$id)
  attachFrames$nodes <- addNodeProps(nodes=attachFrames$nodes,
                                     id_=edge.servers.ids,
                                     properties=list(
                                       volatility=edgeVolatility.values,
                                       radio=rep(x="fog,edge,cloud",
                                                 times=length(edge.servers.ids))))
  
  
  ###### Set link reliability as (1-volatility) of the attached nodes ######
  froms <- c()
  tos <- c()
  linkReliabs <- c()
  for (row in 1:nrow(attachFrames$nodes)) {
    linkFrom <- as.vector(attachFrames$links[row,]$from)
    linkTo <- as.vector(attachFrames$links[row,]$to)
    
    # Link for a fog node
    fromFog <- length(grep(fogIdPrefix, linkFrom)) > 0
    toFog <- length(grep(fogIdPrefix, linkTo)) > 0
    fromEdge <- length(grep(edgeIdPrefix, linkFrom)) > 0
    toEdge <- length(grep(edgeIdPrefix, linkTo)) > 0
    if (fromFog || toFog) {
      froms <- c(froms, linkFrom)
      tos <- c(tos, linkTo)
      fogId <- ifelse(fromFog, yes = linkFrom, no = linkTo)
      linkReliabs <- c(linkReliabs,
                       1 - as.vector(attachFrames$nodes[
                              which(attachFrames$nodes$id == fogId),]$volatility))
    } else if (fromEdge || fromFog) { # Link for an edge node
      froms <- c(froms, linkFrom)
      tos <- c(tos, linkTo)
      edgeId <- ifelse(fromEdge, yes = linkFrom, no = linkTo)
      linkReliabs <- c(linkReliabs,
                       1 - as.vector(attachFrames$nodes[
                             which(attachFrames$nodes$id == edgeId),]$volatility))
    }
  }
  attachFrames$links <- addLinkProps(links = as.data.frame(attachFrames$links),
                                     from_ = froms, to_ = tos,
                                     properties = list(
                                       reliability = linkReliabs))
  
  
  # Write the graphcellsIdx
  g = igraph::graph_from_data_frame(attachFrames$links,
                                    vertices=attachFrames$nodes,
                                    directed = FALSE)
  igraph::write_graph(graph = g, file = out, format = "gml")
}


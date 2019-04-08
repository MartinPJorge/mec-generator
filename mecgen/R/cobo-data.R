#' Generated antennas data for Cobo Calleja
#'
#' Antenna cells generated within the industrial area of Cobo Calleja in Madrid
#'
#' @docType data
#'
#' @usage data(cobo)
#'
#' @format An object of class \code{"list"};  TODO - specify better the format
#'
#' @keywords datasets
#'
#' @source \href{https://ieeexplore.ieee.org/document/8667076}{Region 5g MEC case study}
#'
#' @examples
#' data(cobo)
#' assocs <- build5GScenario(lats = cobo$lat, lons = cobo$lon)
#' 
#' m1Assoc <- assocs[[1]]
#' m1Coords <- assocs[[2]]   
#' m1AccAssocs <- assocs[[3]]   
#' accCentCoords <- assocs[[4]]   
#' m2Assocs <- assocs[[5]]   
#' m2Switches <- assocs[[6]]   
#' m2AggAssocs <- assocs[[7]]   
#' aggCentCoords <- assocs[[8]]   
#' m3Assocs <- assocs[[9]]   
#' m3Switches <- assocs[[10]]  
#' 
#' 
#' write5GtoGraph(m1Assoc, m1Coords, m1AccAssocs, accCentCoords,
#'                            m2Assocs, m2Switches, m2AggAssocs, aggCentCoords,
#'                            m3Assocs, m3Switches, format = "gml",
#'                file = "/tmp/5g-mec.gml")
"cobo"

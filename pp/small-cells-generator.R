source("gen-utils-clean.R")
library(spatstat)
library(rjson)


REGION <- "Cobo-Calleja"
REGIONS <- "../data/regions.json"


# Parse arguments if existing to change default global variables
args <- commandArgs(trailingOnly=TRUE)
if(length(args) > 0)
  if(length(args) != 1) {
    stop(cat("Arguments to receive are: ",
      "regionId\n"))
  } else {
    REGION <- args[1]
  }

region <- selectRegion(regionsFile = REGIONS, regionName = REGION)
regionWin <- owin(xrange = c(region$bl$lon, region$tr$lon),
                   yrange = c(region$bl$lat, region$tr$lat))
antennasCSV <- read.csv(paste("../data/antennas/", REGION, "/", REGION, ".csv",
                           sep = ""))
lteAntennas <- subset(antennasCSV, antennasCSV$radio == "LTE")


# Generate small cells around LTE antennas
extendedAntennas <- lteAntennas
for (row in 1:nrow(lteAntennas)) {
  popu <- findAssocPopulation(regionList = region, lon = lteAntennas[row,]$lon,
                      lat = lteAntennas[row,]$lat)
  smallCellInt <- genSmallCellManta(lteLon = lteAntennas[row,]$lon,
                             lteLat = lteAntennas[row,]$lat,
                             b = popu$smallCells$b, c = popu$smallCells$c,
                             avgSmallCells = popu$smallCells$avgSmallCells)
  
  popou
  popuRect <- outerRect(centLon = popuCent$lon, centLat = popuCent$lat,
                        radius = popu$radius)
  popuWin <- owin(xrange = c(popuRect$lonL, popuRect$lonR),
                  yrange = c(popuRect$latB, popuRect$latT))
  smallCells <- jorgeMaternIImapI(lambda = smallCellInt, win = popuWin,
                                  r = popu$smallCells$repulsion)
  
  extendedAntennas <- appendSmallCells(antennasDf = extendedAntennas,
                                       smallCellsPP = smallCells)
}




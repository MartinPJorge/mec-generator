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


# Generate small-cell antennas arround LTE ones
extendedAntennas <- lteAntennas
genSmallCells <- c()
for (row in 1:nrow(lteAntennas)) {
  print(row)
  # Find assoc population and generate intensity function for the antenna
  popu <- findAssocPopulation(regionList = region, lon = lteAntennas[row,]$lon,
                      lat = lteAntennas[row,]$lat)
  smallCellsManta <- genMatIImapIsmallCellManta(
    lteCenters = data.frame(lon = lteAntennas[row,]$lon,
                            lat = lteAntennas[row,]$lat),
    bs = c(popu$smallCells$b), cs = c(popu$smallCells$c),
    rs = c(popu$smallCells$repulsion),
    avgSmallCells = c(popu$smallCells$avgSmallCells), approxM = "approx3")
  
  # Generate the small-cells
  stepFrad <- popu$smallCells$b*2 + popu$smallCells$c
  antennaRect <- outerRect(centLon = lteAntennas[row,]$lon,
                           centLat = lteAntennas[row,]$lat,
                           radius = stepFrad)
  antennaWin <- owin(xrange = c(antennaRect$lonL, antennaRect$lonR),
                  yrange = c(antennaRect$latB, antennaRect$latT))
  smallCells <- jorgeMaternIImapI(lambda = smallCellsManta, win = antennaWin,
                                  r = popu$smallCells$repulsion)
  genSmallCells <- c(genSmallCells, smallCells$n)
  
  # Append the generated small-cell antennas
  extendedAntennas <- appendSmallCells(antennasDf = extendedAntennas,
                                       smallCellsPP = smallCells)
}




